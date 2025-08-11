/*******************************************************************************
 * Copyright (c) 2025 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 * 
 * SPDX-License-Identifier: EPL-2.0 
 *******************************************************************************/
package org.eclipse.elk.alg.layered.p1cycles;

import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import org.eclipse.elk.alg.layered.LayeredPhases;
import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.Tarjan;
import org.eclipse.elk.alg.layered.intermediate.IntermediateProcessorStrategy;
import org.eclipse.elk.alg.layered.options.GroupOrderStrategy;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.core.alg.ILayoutPhase;
import org.eclipse.elk.core.alg.LayoutProcessorConfiguration;
import org.eclipse.elk.core.util.IElkProgressMonitor;

import com.google.common.collect.Lists;

/**
 * This Cycle Breaking Strategy relies on Tarjan's algorithm to find strongly connected components.
 * It than selects the node with the maximum model order in the strongly connected components and reverses its out-going
 * edges to nodes in the strongly connected component.
 * 
 * <dl>
 *   <dt>Precondition:</dt>
 *      <dd>no self loops</dd>
 *   <dt>Postcondition:</dt>
 *      <dd>the graph has no cycles</dd>
 * </dl>
 *
 */
public abstract class SCCModelOrderCycleBreaker implements ILayoutPhase<LayeredPhases, LGraph> {
    
    /**
     * List of strongly connected components calculated by tarjan.
     */
    protected List<Set<LNode>> stronglyConnectedComponents = new LinkedList<Set<LNode>>();
    
    /**
     * Maps node to id of its strongly connected component.
     */
    protected HashMap <LNode,Integer> nodeToSCCID = new HashMap<>();
    
    /**
     * The edges to reverse.
     */
    protected List<LEdge> revEdges = Lists.newArrayList();
    
    /**
     * The graph.
     */
    protected LGraph graph;

    /** intermediate processing configuration. */
    private static final LayoutProcessorConfiguration<LayeredPhases, LGraph> INTERMEDIATE_PROCESSING_CONFIGURATION =
        LayoutProcessorConfiguration.<LayeredPhases, LGraph>create()
            .addAfter(LayeredPhases.P5_EDGE_ROUTING, IntermediateProcessorStrategy.REVERSED_EDGE_RESTORER);

    @Override
    public LayoutProcessorConfiguration<LayeredPhases, LGraph> getLayoutProcessorConfiguration(final LGraph graph) {
        return INTERMEDIATE_PROCESSING_CONFIGURATION;
    }


    @Override
    public void process(final LGraph layeredGraph, final IElkProgressMonitor monitor) {
        monitor.begin("Model order cycle breaking", 1);
        
        this.graph = layeredGraph;

        // gather edges that point to the wrong direction
        revEdges = Lists.newArrayList();

        // One needs an offset to make sure that the model order of nodes with port constraints is
        // always lower/higher than that of other nodes.
        // E.g. A node with the LAST constraint needs to have a model order m = modelOrder + offset
        // such that m > m(n) with m(n) being the model order of a normal node n (without constraints).
        // Such that the highest model order has to be used as an offset
        int offset = Math.max(layeredGraph.getLayerlessNodes().size(), layeredGraph.getProperty(InternalProperties.MAX_MODEL_ORDER_NODES));
        
        while (true) {
            Tarjan t = new Tarjan(revEdges, stronglyConnectedComponents, nodeToSCCID);
            t.resetTarjan(layeredGraph);
            t.tarjan(layeredGraph);

            // If no Strongly connected components remain, the graph is acyclic.
            if (stronglyConnectedComponents.size() == 0) {
                break;
            }

            // highest model order only incoming
            findNodes(offset, offset * layeredGraph.getProperty(InternalProperties.CB_NUM_MODEL_ORDER_GROUPS));

            // reverse the gathered edges
            for (LEdge edge : revEdges) {
                edge.reverse(layeredGraph, false);
                edge.getSource().getNode().setProperty(LayeredOptions.LAYERING_LAYER_ID,
                        edge.getSource().getNode().getProperty(LayeredOptions.LAYERING_LAYER_ID) + 1);
                layeredGraph.setProperty(InternalProperties.CYCLIC, true);
            }

            stronglyConnectedComponents.clear();
            nodeToSCCID.clear();
            revEdges.clear();
        }

        monitor.done();
        monitor.log("Execution Time: " + monitor.getExecutionTime());
    }

    /**
     * Find the nodes with the highest model order or group model order to reverse all its outgoing edges.
     * @param offset Helper value to calculate constraint partitions. This has to be higher than model order such that
     *          "offset * primary criterion + secondary criterion" works.
     */
    public void findNodes(int offset, int bigOffset) {
        // All strongly connected components have one maximum element for which we can reverse all outgoing edges.
        for (int i = 0; i < stronglyConnectedComponents.size(); i++) {
            LNode max = null;
            GroupModelOrderCalculator calculator = new GroupModelOrderCalculator();
            int maxModelOrder = Integer.MIN_VALUE;
            for (LNode n : stronglyConnectedComponents.get(i)) {
                // Check whether model order or the group model order is the primary criterion.
                // If it is group model order, we need to handle this differently.
                boolean enforceGroupModelOrder = this.graph.getProperty(
                        LayeredOptions.CONSIDER_MODEL_ORDER_GROUP_MODEL_ORDER_CB_GROUP_ORDER_STRATEGY) == GroupOrderStrategy.ENFORCED;
                if (max == null) {
                    // Case first element
                    max = n;
                    maxModelOrder = enforceGroupModelOrder
                            ? calculator.computeConstraintGroupModelOrder(n, bigOffset, offset)
                            : calculator.computeConstraintModelOrder(n, offset);
                } else {
                    // Find a new maximum if possible.
                    int modelOrderCurrent = enforceGroupModelOrder
                            ? calculator.computeConstraintGroupModelOrder(n, bigOffset, offset)
                            : calculator.computeConstraintModelOrder(n, offset);
                    if (maxModelOrder < modelOrderCurrent) {
                        max = n;
                        maxModelOrder = modelOrderCurrent;
                    }
                }
            }
            for (LEdge edge : max.getOutgoingEdges()) {
                // Reverse all edges to the same strongly connected component.
                if (stronglyConnectedComponents.get(i).contains(edge.getTarget().getNode())) {
                    revEdges.add(edge);
                }
            }
        }
    }
}