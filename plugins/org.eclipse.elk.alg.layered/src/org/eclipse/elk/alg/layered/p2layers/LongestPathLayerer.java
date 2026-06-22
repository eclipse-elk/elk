/*******************************************************************************
 * Copyright (c) 2010, 2020 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.alg.layered.p2layers;

import java.util.List;
import java.util.Stack;

import org.eclipse.elk.alg.layered.LayeredPhases;
import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.intermediate.IntermediateProcessorStrategy;
import org.eclipse.elk.core.alg.ILayoutPhase;
import org.eclipse.elk.core.alg.LayoutProcessorConfiguration;
import org.eclipse.elk.core.util.IElkProgressMonitor;

/**
 * The most basic layering algorithm, which assign layers according to the
 * longest path to a sink.
 * 
 * <dl>
 *   <dt>Precondition:</dt><dd>the graph has no cycles</dd>
 *   <dt>Postcondition:</dt><dd>all nodes have been assigned a layer such that
 *     edges connect only nodes from layers with increasing indices</dd>
 * </dl>
 */
public final class LongestPathLayerer implements ILayoutPhase<LayeredPhases, LGraph> {
    
    /** intermediate processing configuration. */
    private static final LayoutProcessorConfiguration<LayeredPhases, LGraph> BASELINE_PROCESSING_CONFIGURATION =
        LayoutProcessorConfiguration.<LayeredPhases, LGraph>create()
            .addBefore(LayeredPhases.P1_CYCLE_BREAKING,
                    IntermediateProcessorStrategy.EDGE_AND_LAYER_CONSTRAINT_EDGE_REVERSER)
            .addBefore(LayeredPhases.P2_LAYERING, IntermediateProcessorStrategy.LAYER_CONSTRAINT_PREPROCESSOR)
            .addBefore(LayeredPhases.P3_NODE_ORDERING, IntermediateProcessorStrategy.LAYER_CONSTRAINT_POSTPROCESSOR);
    
    /** the layered graph to which layers are added. */
    private LGraph layeredGraph;
    /** map of nodes to their height in the layering. */
    private int[] nodeHeights;
    
    @Override
    public LayoutProcessorConfiguration<LayeredPhases, LGraph> getLayoutProcessorConfiguration(final LGraph graph) {
        return BASELINE_PROCESSING_CONFIGURATION;
    }
    
    @Override
    public void process(final LGraph thelayeredGraph, final IElkProgressMonitor monitor) {
        monitor.begin("Longest path layering", 1);
        
        layeredGraph = thelayeredGraph;
        List<LNode> nodes = layeredGraph.getLayerlessNodes();
        
        // initialize values required for the computation
        nodeHeights = new int[nodes.size()];
        int index = 0;
        for (LNode node : nodes) {
            // the node id is used as index for the nodeHeights array
            node.id = index;
            nodeHeights[index] = -1;
            index++;
        }
        
        // process all nodes
        for (LNode node : nodes) {
            visit(node);
        }
        // empty the list of unlayered nodes
        nodes.clear();
        
        // release the created resources
        this.layeredGraph = null;
        this.nodeHeights = null;
        
        monitor.done();
    }
    

    /**
     * Visit a node: if not already visited, find the longest path to a sink.
     * 
     * @param node node to visit
     */
    private void visit(final LNode node) {
        
        int [] internalNodeHeights = new int [nodeHeights.length];
        for (int i = 0; i < internalNodeHeights.length; i++) {
            internalNodeHeights[i] = nodeHeights[i];
        }
        Stack<LNode> stack = new Stack<>();
        stack.push(node);
        while (stack.size() > 0) { 
            LNode top = stack.pop();
            boolean outgoingEdges = false;
            for (LPort port : top.getPorts()) {
                for (LEdge edge : port.getOutgoingEdges()) {
                    LNode targetNode = edge.getTarget().getNode();
                    outgoingEdges = true;
                    if (node != targetNode && top != targetNode && internalNodeHeights[targetNode.id] < 0) {
                        // push unsolved problem back onto stack
                        stack.push(top);
                        stack.push(targetNode);
                    } else {
                        internalNodeHeights[top.id] = 
                                Math.max(internalNodeHeights[top.id], internalNodeHeights[targetNode.id] + 1);
                    }
                }
            }
            if (!outgoingEdges) {
                internalNodeHeights[top.id] = 1;
            }
        }
        // I think this is slower than the recursive function, because we have to 
        // call the outer visit method for each node and are not storing intermediate results
        // TODO: this should be able to be improved
        // it's not trivial though, because we need to keep track of the longest path somehow
        putNode(node, internalNodeHeights[node.id]);
    }
    
    /**
     * Puts the given node into the layered graph, adding new layers as necessary.
     * 
     * @param node a node
     * @param height height of the layer where the node shall be added
     *          (height = number of layers - layer index)
     */
    private void putNode(final LNode node, final int height) {
        List<Layer> layers = layeredGraph.getLayers();
        
        // add layers so as to guarantee that number of layers >= height
        for (int i = layers.size(); i < height; i++) {
            layers.add(0, new Layer(layeredGraph));
        }
        
        // layer index = number of layers - height
        node.setLayer(layers.get(layers.size() - height));
        nodeHeights[node.id] = height;
    }

}
