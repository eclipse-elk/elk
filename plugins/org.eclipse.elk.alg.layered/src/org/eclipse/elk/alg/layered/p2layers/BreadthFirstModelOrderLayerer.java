/*******************************************************************************
 * Copyright (c) 2023 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.alg.layered.p2layers;

import java.util.ArrayList;
import java.util.Collections;
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
 * A model order layer assigner that assumes a breadth first node ordering.
 * 
 * <dl>
 *   <dt>Precondition:</dt>
 *      <dd>The graph has no cycles.</dd>
 *      <dd>Each node is a normal node and has a model order set or is a label dummy node.</dd>
 *   <dt>Postcondition:</dt>
 *      <dd>All nodes have been assigned a layer such that
 *        edges connect only nodes from layers with increasing indices.</dd>
 *      <dd>No node with a higher model order is in a layer before a node with a lower model order.</dd>
 * </dl>
 */
public final class BreadthFirstModelOrderLayerer implements ILayoutPhase<LayeredPhases, LGraph> {
    
    /** intermediate processing configuration. */
    private static final LayoutProcessorConfiguration<LayeredPhases, LGraph> BASELINE_PROCESSING_CONFIGURATION =
        LayoutProcessorConfiguration.<LayeredPhases, LGraph>create()
            .addBefore(LayeredPhases.P1_CYCLE_BREAKING,
                    IntermediateProcessorStrategy.EDGE_AND_LAYER_CONSTRAINT_EDGE_REVERSER)
            .addBefore(LayeredPhases.P2_LAYERING, IntermediateProcessorStrategy.LAYER_CONSTRAINT_PREPROCESSOR)
            .addBefore(LayeredPhases.P3_NODE_ORDERING, IntermediateProcessorStrategy.LAYER_CONSTRAINT_POSTPROCESSOR);
    
    /** the layered graph to which layers are added. */
    private LGraph layeredGraph;
    
    @Override
    public LayoutProcessorConfiguration<LayeredPhases, LGraph> getLayoutProcessorConfiguration(final LGraph graph) {
        return BASELINE_PROCESSING_CONFIGURATION;
    }
    
    @Override
    public void process(final LGraph thelayeredGraph, final IElkProgressMonitor monitor) {
        monitor.begin("Breadth first model order layering", 1);
        
        layeredGraph = thelayeredGraph;
        // Add first layer.
        List<LNode> realNodes = new ArrayList<LNode>();
        for (LNode node : layeredGraph.getLayerlessNodes()) {
            if (node.getType() == NodeType.NORMAL) {
                realNodes.add(node);
            }
        }
        // Sort real nodes by model order.
        // FIXME group model order?
        Collections.sort(realNodes, (n1, n2) -> {
            if (n1.hasProperty(InternalProperties.MODEL_ORDER) && n2.hasProperty(InternalProperties.MODEL_ORDER)) {
                return Integer.compare(n1.getProperty(InternalProperties.MODEL_ORDER),
                        n2.getProperty(InternalProperties.MODEL_ORDER));
            }
            throw new UnsupportedGraphException("The BF model order layer assigner requires all real nodes to have"
                    + " a model order.");
        });
        boolean firstNode = true;
        Layer currentLayer = new Layer(layeredGraph);
        Layer currentDummyLayer = null;
        layeredGraph.getLayers().add(currentLayer);
        for (LNode node : realNodes) {
            if (firstNode) {
                // Set the first node in the first layer.
                node.setLayer(currentLayer);
                firstNode = false;
            } else {
                // Check for each incoming edge whether it connects to the current layer.
                // If this is the case, add a new dummy and normal layer.
                for (LEdge edge : node.getIncomingEdges()) {
                    if (edge.getSource().getNode().getType() == NodeType.NORMAL // Case nodes are directly connected
                            && edge.getSource().getNode().getLayer() == currentLayer
                        || edge.getSource().getNode().getType() == NodeType.LABEL // Case dummy label in-between.
                            && edge.getSource().getNode().getIncomingEdges().iterator().next()
                                .getSource().getNode().getLayer() == currentLayer) {
                        currentDummyLayer = new Layer(layeredGraph);
                        layeredGraph.getLayers().add(currentDummyLayer);
                        currentLayer = new Layer(layeredGraph);
                        layeredGraph.getLayers().add(currentLayer);
                    }
                }
                // Add all label dummies to the dummy layer between the last layer and the current one.
                for (LEdge edge : node.getIncomingEdges()) {
                    if (edge.getSource().getNode().getType() == NodeType.LABEL
                            && edge.getSource().getNode().getLayer() == null) {
                        edge.getSource().getNode().setLayer(currentDummyLayer);
                    }
                }
                // Add the node to the (potentially new) current layer.
                node.setLayer(currentLayer);
            }
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

}
