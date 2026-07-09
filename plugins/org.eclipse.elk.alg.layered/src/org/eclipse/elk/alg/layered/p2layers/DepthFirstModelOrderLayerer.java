/*******************************************************************************
 * Copyright (c) 2023 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.

 * This Source Code may also be made available under the following Secondary
 * Licenses when the conditions for such availability set forth in the Eclipse
 * Public License v. 2.0 are satisfied: GPL-3.0 which is available at
 * https://www.gnu.org/licenses/gpl-3.0-standalone.html.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.layered.p2layers;

import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;

import org.eclipse.elk.alg.layered.LayeredPhases;
import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LNode.NodeType;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.intermediate.IntermediateProcessorStrategy;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.core.UnsupportedGraphException;
import org.eclipse.elk.core.alg.ILayoutPhase;
import org.eclipse.elk.core.alg.LayoutProcessorConfiguration;
import org.eclipse.elk.core.util.IElkProgressMonitor;

/**
 * A model order layer assigner that assumes a depth first node ordering.
 * Places nodes depth first by model order in strips as seen below.
 * 1 -> 2 -> 3
 *        \
 *      4 -> 5 -> 6
 * 
 * <dl>
 *   <dt>Precondition:</dt>
 *      <dd>The graph has no cycles.</dd>
 *      <dd>Each node is a normal node and has a model order set or is a label dummy node.</dd>
 *   <dt>Postcondition:</dt>
 *      <dd>All nodes have been assigned a layer such that
 *        edges connect only nodes from layers with increasing indices.</dd>
 * </dl>
 */
public final class DepthFirstModelOrderLayerer implements ILayoutPhase<LayeredPhases, LGraph> {
    
    /** intermediate processing configuration. */
    private static final LayoutProcessorConfiguration<LayeredPhases, LGraph> BASELINE_PROCESSING_CONFIGURATION =
        LayoutProcessorConfiguration.<LayeredPhases, LGraph>create()
            .addBefore(LayeredPhases.P1_CYCLE_BREAKING,
                    IntermediateProcessorStrategy.EDGE_AND_LAYER_CONSTRAINT_EDGE_REVERSER)
            .addBefore(LayeredPhases.P2_LAYERING, IntermediateProcessorStrategy.LAYER_CONSTRAINT_PREPROCESSOR)
            .addBefore(LayeredPhases.P3_NODE_ORDERING, IntermediateProcessorStrategy.LAYER_CONSTRAINT_POSTPROCESSOR);
    
    /** The layered graph to which layers are added. */
    private LGraph layeredGraph;
    
    /** The current layer without the nodesToPlace. */
    private Layer currentLayer;

    /** The current layer id. */
    private int currentLayerId;

    /** The current dummy layer without the nodesToPlace. */
    private Layer currentDummyLayer;

    /**
     * Holds nodes that will be placed once a new strip is connected to the rest of the graph.
     * The node id holds the currently desired layer.
     */
    private List<LNode> nodesToPlace;
    
    /**
     * Saves the maximum layer the nodes that still have to be assigned to a layer will be placed in.
     * This is necessary to determine the actual offset of the nodes to place.
     */
    private int maxToPlace = 0;
    
    @Override
    public LayoutProcessorConfiguration<LayeredPhases, LGraph> getLayoutProcessorConfiguration(final LGraph graph) {
        return BASELINE_PROCESSING_CONFIGURATION;
    }
    
    @Override
    public void process(final LGraph thelayeredGraph, final IElkProgressMonitor monitor) {
        monitor.begin("Depth first model order layering", 1);
        
        layeredGraph = thelayeredGraph;
        // Add first layer.
        List<LNode> realNodes = new ArrayList<LNode>();
        for (LNode node : layeredGraph.getLayerlessNodes()) {
            if (node.getType() == NodeType.NORMAL) {
                realNodes.add(node);
            }
        }
        // Sort real nodes by model order.
        Collections.sort(realNodes, (n1, n2) -> {
            if (n1.hasProperty(InternalProperties.MODEL_ORDER) && n2.hasProperty(InternalProperties.MODEL_ORDER)) {
                return Integer.compare(n1.getProperty(InternalProperties.MODEL_ORDER),
                        n2.getProperty(InternalProperties.MODEL_ORDER));
            }
            throw new UnsupportedGraphException("The DF model order layer assigner requires all real nodes to have"
                    + " a model order.");
        });
        boolean firstNode = true;
        currentLayer = new Layer(layeredGraph);
        currentDummyLayer = null;
        layeredGraph.getLayers().add(currentLayer);
        currentLayer.id = 0;
        currentLayerId = 0;
        nodesToPlace = new LinkedList<>();
        for (LNode node : realNodes) {
            if (firstNode) {
                // Set the first node in the first layer.
                node.setLayer(currentLayer);
                firstNode = false;
            } else {
                // Check whether the current node connects to the current layer.
                // If this is the case, place it in the layer after the current one (and add a dummy layer).
                // If not, begin a new "strip" based on the incoming connection of the node.
                if (isConnectedToCurrentLayer(node)) {
                    // Check whether there are incoming edges to layers higher than currendLayer.id + 2
                    int maxLayer = currentLayerId;
                    maxLayer = getMaxConnectedLayer(maxLayer, node);
                    int desiredLayer = maxLayer + 2;
                    int layerDiff = maxLayer - currentLayerId;
                    // Move all nodes that where not already placed by the layer difference.
                    if (!nodesToPlace.isEmpty()) {
                        if (layerDiff > 0) {
                            // Case some dependency to the existing graph was found add all nodes to their layers.
                            for (LNode toPlace : nodesToPlace) {
                                toPlace.id += maxLayer - maxToPlace;
                            }
                            placeNodesToPlace();
                            // Remove all nodes since they are now placed correctly.
                            nodesToPlace.clear();
                            addNodeToLayer(desiredLayer, node);
                        } else {
                            // Case nodes cannot be placed without the need to move them later.
                            nodesToPlace.add(node);
                            node.id = desiredLayer;
                            maxToPlace = Math.max(maxToPlace, desiredLayer);
                            // Add dummy nodes and give them their desired layer id.
                            for (LEdge edge : node.getIncomingEdges()) {
                                if (edge.getSource().getNode().getLayer() == null
                                        && edge.getSource().getNode().getType() == NodeType.LABEL) {
                                    nodesToPlace.add(edge.getSource().getNode());
                                    edge.getSource().getNode().id = desiredLayer - 1;
                                }
                            }
                            // Update the current layer id.
                            currentLayerId = desiredLayer;
                        }
                    } else {
                        addNodeToLayer(desiredLayer, node);
                    }
                } else {
                    // Case a new strip has to begin.
                    placeNodesToPlace();
                    nodesToPlace.clear();
                    
                    // Find the layer for the first element of the strip.
                    // If it has no incoming connections it is placed in the first layer.
                    // Once a node with a connection to the previous strip is added previous nodes move to be moved,
                    // therefore, we need to safe the dependent nodes.
                    // If it has an incoming connection, it has to be to an already placed node and hence the new node
                    // Should be placed in the current layer after the connected node.
                    int desiredLayer = 0;
                    if (!node.getIncomingEdges().iterator().hasNext()) {
                        // Case no incoming connections, save the node to be placed.
                        nodesToPlace.add(node);
                        node.id = 0; // Save the layer the node will be placed in.
                        maxToPlace = Math.max(maxToPlace, 0);
                        currentLayer = layeredGraph.getLayers().get(0);
                        currentLayerId = 0;
                        
                    } else {
                        // Find the last layer the node is connected to.
                        int maxLayer = 0;
                        maxLayer = getMaxConnectedLayer(maxLayer, node);
                        desiredLayer = maxLayer + 2;
                        // Add node to desired layer.
                        addNodeToLayer(desiredLayer, node);
                    }
                }
            }
        }
        // Finish the last 
        if (!nodesToPlace.isEmpty()) {
            placeNodesToPlace();
        }
        // Empty the list of unlayered nodes.
        layeredGraph.getLayerlessNodes().clear();
        // Delete all empty dummy layers.
        List<Layer> toDelete = new ArrayList<>();
        for (Layer layer : layeredGraph.getLayers()) {
            if (layer.getNodes().isEmpty()) {
                toDelete.add(layer);
            }
        }
        layeredGraph.getLayers().removeAll(toDelete);
        // Adjust layer ids
        int layerId = 0;
        for (Layer layer : layeredGraph.getLayers()) {
            layer.id = layerId;
            layerId++;
        }
        monitor.done();
    }
    
    /**
     * Returns true if a node is connected (directly or by label dummy) to the current layer.
     * @param node The node
     * @return True if the node is connected to the current layer.
     */
    private boolean isConnectedToCurrentLayer(final LNode node) {
        for (LEdge edge : node.getIncomingEdges()) {
            boolean directlyConnected;
            boolean connectedViaLabelDummy;
            if (nodesToPlace.isEmpty()) {
                // Case the node to check already has a layer.
                // This may case an NPE if getLayer is null;
                directlyConnected = edge.getSource().getNode().getType() == NodeType.NORMAL
                        && edge.getSource().getNode().getLayer() != null
                        && edge.getSource().getNode().getLayer().id == currentLayerId;
                if (edge.getSource().getNode().getIncomingEdges().iterator().hasNext()) {
                    Layer connectedLayerViaDummy = edge.getSource().getNode().getIncomingEdges().iterator().next()
                            .getSource().getNode().getLayer();
                    connectedViaLabelDummy = edge.getSource().getNode().getType() == NodeType.LABEL
                            && connectedLayerViaDummy != null
                            && connectedLayerViaDummy.id == currentLayerId;
                } else {
                    connectedViaLabelDummy = false;
                }
            } else {
                // Case the node to check might not be placed yet and the node id holds its desired layer.
                directlyConnected = edge.getSource().getNode().getType() == NodeType.NORMAL
                        && edge.getSource().getNode().id == currentLayerId;
                connectedViaLabelDummy = edge.getSource().getNode().getType() == NodeType.LABEL
                        && edge.getSource().getNode().getIncomingEdges().iterator().next()
                        .getSource().getNode().id == currentLayerId;
            }
            if (directlyConnected || connectedViaLabelDummy) {
                return true;
            }
        }
        return false;
    }
    
    /**
     * Adds a node and its label dummies to the desired layer.
     * @param layerId The id of the desired layer.
     * @param node The node.
     */
    private void addNodeToLayer(final int layerId, final LNode node) {
        if (layerId < layeredGraph.getLayers().size()) {
            // Case the layers already exist.
            currentLayer = layeredGraph.getLayers().get(layerId);
            currentDummyLayer = layeredGraph.getLayers().get(layerId - 1);
            currentLayerId = layerId;
        } else {
            // Case the layers need to be added.
            currentDummyLayer = new Layer(layeredGraph);
            currentDummyLayer.id = layerId - 1;
            layeredGraph.getLayers().add(currentDummyLayer);
            currentLayer = new Layer(layeredGraph);
            currentLayer.id = layerId;
            layeredGraph.getLayers().add(currentLayer);
            currentLayerId = layerId;
        }
        // Add node
        node.setLayer(currentLayer);
        // Add dummy node
        for (LEdge edge : node.getIncomingEdges()) {
            if (edge.getSource().getNode().getLayer() == null
                    && edge.getSource().getNode().getType() == NodeType.LABEL) {
                edge.getSource().getNode().setLayer(currentDummyLayer);
            }
        }
    }
    
    /**
     * Returns the highest layer id a node is connected to.
     * @param maxLayer The currently highest layer id.
     * @param node The node.
     * @return The highest layer id the node is connected to.
     */
    private int getMaxConnectedLayer(final int layerId, final LNode node) {
        int maxLayer = layerId;
        for (LEdge edge : node.getIncomingEdges()) {
            if (edge.getSource().getNode().getLayer() != null) {
                maxLayer = Math.max(maxLayer, edge.getSource().getNode().getLayer().id);
            }
        }
        return maxLayer;
    }
    
    /**
     * Places the nodes to place in their desired layer.
     */
    public void placeNodesToPlace() {
        maxToPlace = 0;
        for (LNode nodeToPlace : nodesToPlace) {
            if (nodeToPlace.id >= layeredGraph.getLayers().size()) {
                // Add a normal and a dummy layer.
                Layer dummyLayer = new Layer(layeredGraph);
                dummyLayer.id = nodeToPlace.id - 1;
                layeredGraph.getLayers().add(dummyLayer);
                Layer newLayer = new Layer(layeredGraph);
                newLayer.id = nodeToPlace.id;
                layeredGraph.getLayers().add(newLayer);
            }
            nodeToPlace.setLayer(layeredGraph.getLayers().get(nodeToPlace.id));
        }
    }

}
