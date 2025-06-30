/*******************************************************************************
 * Copyright (c) 2024 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 * 
 * SPDX-License-Identifier: EPL-2.0 
 *******************************************************************************/
package org.eclipse.elk.alg.layered.intermediate.unzipping;

import java.util.ArrayList;
import java.util.List;
import java.util.ListIterator;

import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LNode.NodeType;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.intermediate.LongEdgeSplitter;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.core.alg.ILayoutProcessor;
import org.eclipse.elk.core.options.PortConstraints;
import org.eclipse.elk.core.util.IElkProgressMonitor;
import org.eclipse.elk.core.util.Pair;

import com.google.common.collect.Lists;

/**
 * Divides nodes up between layers to create a more compact final layout.
 * Reads the property of each layer to determine how many sub-layers it 
 * should be split into.
 * 
 * <dl>
 * <dt>Preconditions:</dt>
 *   <dd>A layered graph whose node order has been decided.</dd>
 * <dt>Postconditions:</dt>
 *   <dd>Layers are split up into multiple layers with the nodes alternating between them. For example, if layerSplit
 *   is set to 3 and there are 5 nodes in a layer, then node 1 is placed in sublayer 1, node 2 in sublayer 2, node 3 in
 *   sublayer 3, node 4 in sublayer 1 and node 5 in sublayer 2.</dd>
 * <dt>Slots:</dt>
 *   <dd>Before phase 4.</dd>
 * <dt>Same-slot dependencies:</dt>
 *   <dd>None</dd>
 * </dl>
 *
 */
public class AlternatingLayerUnzipper implements ILayoutProcessor<LGraph> {
    
    @Override
    public void process (LGraph graph, IElkProgressMonitor progressMonitor) {
        
        int insertionLayerOffset = 1;
        List<Pair<Layer, Integer>> newLayers = new ArrayList<>();
        for (int i = 0; i < graph.getLayers().size(); i++) {
            
            int N = getLayerSplitProperty(graph.getLayers().get(i));
            boolean resetOnLongEdges = getResetOnLongEdgesProperty(graph.getLayers().get(i));
            boolean minimizeEdgeLength = getminimizeEdgeLengthProperty(graph.getLayers().get(i));
            
            if (minimizeEdgeLength) {
                // get maximum node width and average node height
                double maxWidth = 0.0;
                double averageHeight = 0.0;
                
                for (LNode node : graph.getLayers().get(i)) {
                    maxWidth = Math.max(maxWidth, node.getSize().x);
                    averageHeight += node.getSize().y;
                }
                
                averageHeight /= graph.getLayers().get(i).getNodes().size();
                
                // add spacings for heuristic to get closer to zero spacing situation
                maxWidth += Math.max(2 * graph.getProperty(LayeredOptions.SPACING_EDGE_NODE_BETWEEN_LAYERS),
                        Math.max(graph.getLayers().get(i).getNodes().size() 
                                * graph.getProperty(LayeredOptions.SPACING_EDGE_EDGE_BETWEEN_LAYERS),
                                graph.getProperty(LayeredOptions.SPACING_NODE_NODE_BETWEEN_LAYERS)));
                
                averageHeight += Math.max(
                        graph.getProperty(LayeredOptions.SPACING_NODE_NODE),
                        graph.getProperty(LayeredOptions.SPACING_EDGE_NODE));
                
                // apply heuristic
                // TODO: this heuristic is for the decision between N=1 and N=2, it shouldn't be used with N > 2
                if (maxWidth / averageHeight >= graph.getLayers().get(i).getNodes().size() / 4) {
                    // skip this layer and do not split it
                    continue;
                }
                
            }
            
            // only split if there are more nodes than the resulting sub-layers
            // an alternative would be to reduce N for this layer, this may or may
            // not be desirable
            if (graph.getLayers().get(i).getNodes().size() > N) {
                
                List<Layer> subLayers = new ArrayList<>();
                // add current layer as first sub-layer
                subLayers.add(graph.getLayers().get(i));
                for (int j = 0; j < N - 1; j++) {
                    Layer newLayer = new Layer(graph);
                    newLayers.add(new Pair<>(newLayer, i+j+insertionLayerOffset));
                    subLayers.add(newLayer);
                }
                insertionLayerOffset += N - 1;
                
                int nodesInLayer = subLayers.get(0).getNodes().size();
                for (int j = 0, nodeIndex = 0, targetLayer = 0; j < nodesInLayer; j++, nodeIndex++, targetLayer++) {
                    LNode node = subLayers.get(0).getNodes().get(nodeIndex);
                    if (node.getType() != NodeType.NONSHIFTING_PLACEHOLDER) {
                        nodeIndex += shiftNode(graph, subLayers, targetLayer % N, nodeIndex);
                    } else {
                        j -= 1;
                        targetLayer -= 1;
                    }
                    if (resetOnLongEdges && node.getType() == NodeType.LONG_EDGE) {
                        // reset next iterations target layer to 0
                        targetLayer = -1;
                    }
                    
                }
            }
        }
        for (Pair<Layer, Integer> newLayer : newLayers) {
            graph.getLayers().add(newLayer.getSecond(), newLayer.getFirst());
        }
        
        // remove unconnected placeholder nodes
        for (Layer layer : graph.getLayers()) {
            ListIterator<LNode> nodeIterator = layer.getNodes().listIterator();
            while (nodeIterator.hasNext()) {
                LNode node = nodeIterator.next();
                if (node.getType() == NodeType.PLACEHOLDER || node.getType() == NodeType.NONSHIFTING_PLACEHOLDER) {
                    nodeIterator.remove();
                }
            }
        }
        
        
    }

    /**
     * Checks all nodes of a layer for the layerSplit property and returns the lowest set value.
     * 
     * @param layer The layer to determine the layerSplit property for
     * @return the layerSplit value
     */
    private int getLayerSplitProperty(Layer layer) {
        int layerSplit = Integer.MAX_VALUE;
        boolean propertyUnset = true;
        for (int i = 0; i < layer.getNodes().size(); i++) {
            if (layer.getNodes().get(i).hasProperty(LayeredOptions.LAYER_UNZIPPING_LAYER_SPLIT)) {
                propertyUnset = false;
                int nodeValue = layer.getNodes().get(i).getProperty(LayeredOptions.LAYER_UNZIPPING_LAYER_SPLIT);
                layerSplit = layerSplit < nodeValue ? layerSplit : nodeValue;
            }
        }
        if (propertyUnset) {
            layerSplit = LayeredOptions.LAYER_UNZIPPING_LAYER_SPLIT.getDefault();
        }
        return layerSplit;
    }
    
    /**
     * Checks all nodes of a layer for the resetOnLongEdges property and if any sets the value to false, returns false.
     * 
     * @param layer The layer to determine the resetOnLongEdges property for.
     * @return the resetOnLongEdges value
     */
    private boolean getResetOnLongEdgesProperty(Layer layer) {
        for (int i = 0; i < layer.getNodes().size(); i++) {
            if (layer.getNodes().get(i).hasProperty(LayeredOptions.LAYER_UNZIPPING_RESET_ON_LONG_EDGES)) {
                if (!layer.getNodes().get(i).getProperty(LayeredOptions.LAYER_UNZIPPING_RESET_ON_LONG_EDGES)) {
                    return false;
                }
            }
        }
        return true;
    }
    
    /**
     * Checks all nodes of a layer for the minimizeEdgeLength property and if any sets the value to true, returns true.
     * 
     * @param layer The layer to determine the minimizeEdgeLength property for.
     * @return the minimizeEdgeLength value
     */
    private boolean getminimizeEdgeLengthProperty(Layer layer) {
        for (int i = 0; i < layer.getNodes().size(); i++) {
            if (layer.getNodes().get(i).hasProperty(LayeredOptions.LAYER_UNZIPPING_MINIMIZE_EDGE_LENGTH)) {
                if (layer.getNodes().get(i).getProperty(LayeredOptions.LAYER_UNZIPPING_MINIMIZE_EDGE_LENGTH)) {
                    return true;
                }
            }
        }
        return false;
    }

    /**
     * Shifts a node from one layer to another and adds dummy nodes for the long edges this introduces.
     * @param graph
     * @param subLayers
     * @param targetLayer
     * @param nodeIndex
     * 
     * @return the number new nodes in the original layer
     */
    private int shiftNode(LGraph graph, List<Layer> subLayers, int targetLayer, int nodeIndex) {
        LNode node = subLayers.get(0).getNodes().get(nodeIndex);
        if (targetLayer > 0){
            node.setLayer(subLayers.get(targetLayer));
        }
        // handle incoming edges and preceding layers
        int edgeCount = 0;
        // If there are no incoming edges, the nodeindex will have to be decreased by one
        boolean noIncomingEdges = true; 
        List<LEdge> reversedIncomingEdges = Lists.reverse(Lists.newArrayList(node.getIncomingEdges()));
        for (LEdge incomingEdge : reversedIncomingEdges) {
            noIncomingEdges = false;
            LEdge nextEdgeToSplit = incomingEdge;
            for (int layerIndex = 0; layerIndex < targetLayer; layerIndex++) {
                LNode dummyNode = createDummyNode(graph, nextEdgeToSplit);
                if (nodeIndex + edgeCount > subLayers.get(layerIndex).getNodes().size()) {
                    dummyNode.setLayer(subLayers.get(layerIndex));
                } else {
                    dummyNode.setLayer(nodeIndex + edgeCount, subLayers.get(layerIndex)); 
                }
                nextEdgeToSplit = LongEdgeSplitter.splitEdge(nextEdgeToSplit, dummyNode);
            }
            if (targetLayer > 0) {
                edgeCount += 1;
            }
        }
        
        // create unconnected dummy nodes to fill the layers if there are no incoming edges
        if (noIncomingEdges) {
            for (int layerIndex = 0; layerIndex < targetLayer; layerIndex++) {
                LNode dummyNode = new LNode(graph);
                dummyNode.setType(NodeType.PLACEHOLDER);
                if (nodeIndex + edgeCount > subLayers.get(layerIndex).getNodes().size()) {
                    dummyNode.setLayer(subLayers.get(layerIndex));
                } else {
                    dummyNode.setLayer(nodeIndex + edgeCount, subLayers.get(layerIndex));
                }
            }
            if (targetLayer > 0) {
                edgeCount += 1;
            }
        }
        
        // handle outgoing edges and following layers
        boolean extraEdge = false;
        for (LEdge outgoingEdge : node.getOutgoingEdges()) {
            LEdge nextEdgeToSplit = outgoingEdge;
            for (int layerIndex = targetLayer + 1; layerIndex < subLayers.size(); layerIndex++) {
                LNode dummyNode = createDummyNode(graph, nextEdgeToSplit);
                dummyNode.setLayer(subLayers.get(layerIndex));
                nextEdgeToSplit = LongEdgeSplitter.splitEdge(nextEdgeToSplit, dummyNode);
            }

            for (int layerIndex = 0; layerIndex <= targetLayer; layerIndex++) {
                if (extraEdge) {
                    // add a placeholder beneath node's old position so that later
                    LNode placeholder = new LNode(graph);
                    placeholder.setType(NodeType.NONSHIFTING_PLACEHOLDER);
                    
                    if (nodeIndex + 1 > subLayers.get(layerIndex).getNodes().size()) {
                        placeholder.setLayer(subLayers.get(layerIndex));
                    } else {
                        placeholder.setLayer(nodeIndex + 1, subLayers.get(layerIndex));
                    }
                }
            }

            if (extraEdge) {
                edgeCount += 1;
            }
            
            extraEdge = true;
        }

        if (edgeCount > 0) {
            return edgeCount - 1;
        } else {
            return 0;
        }
    }

    /**
     * Creates a dummy node for an edge that should be split into a long edge.
     * @param graph
     * @param nextEdgeToSplit
     * @return
     */
    private LNode createDummyNode(LGraph graph, LEdge nextEdgeToSplit) {
        LNode dummyNode = new LNode(graph);
        dummyNode.setType(NodeType.LONG_EDGE);
        dummyNode.setProperty(InternalProperties.ORIGIN, nextEdgeToSplit);
        dummyNode.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_POS);
        return dummyNode;
    }

}
