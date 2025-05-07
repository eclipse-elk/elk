/*******************************************************************************
 * Copyright (c) 2023 mwr and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 * 
 * SPDX-License-Identifier: EPL-2.0 
 *******************************************************************************/
package org.eclipse.elk.alg.layered.p1cycles;

import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;
import java.util.Stack;

import org.eclipse.elk.alg.layered.LayeredPhases;
import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.intermediate.IntermediateProcessorStrategy;
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
 * @author mwr
 *
 */
public class SCCModelOrderCycleBreaker implements ILayoutPhase<LayeredPhases, LGraph> {

    private int firstSeparateModelOrder;
    private int lastSeparateModelOrder;
    private int index;
    protected List<Set<LNode>> stronglyConnectedComponents = new LinkedList<Set<LNode>>();
    private Stack<LNode> stack = new Stack<LNode>();
    protected List<LEdge> revEdges = Lists.newArrayList();

    /** intermediate processing configuration. */
    private static final LayoutProcessorConfiguration<LayeredPhases, LGraph> INTERMEDIATE_PROCESSING_CONFIGURATION =
        LayoutProcessorConfiguration.<LayeredPhases, LGraph>create()
            .addAfter(LayeredPhases.P5_EDGE_ROUTING, IntermediateProcessorStrategy.REVERSED_EDGE_RESTORER);

    @Override
    public LayoutProcessorConfiguration<LayeredPhases, LGraph> getLayoutProcessorConfiguration(final LGraph graph) {
        // TODO Auto-generated method stub
        return INTERMEDIATE_PROCESSING_CONFIGURATION;
    }
    
    
    @Override
    public void process(final LGraph layeredGraph, final IElkProgressMonitor monitor) {
        monitor.begin("Model order cycle breaking", 1);
        
        // Reset FIRST_SEPARATE and LAST_SEPARATE counters.
        firstSeparateModelOrder = 0;
        lastSeparateModelOrder = 0;
        
        // gather edges that point to the wrong direction
        revEdges = Lists.newArrayList();
        
        // One needs an offset to make sure that the model order of nodes with port constraints is
        // always lower/higher than that of other nodes.
        // E.g. A node with the LAST constraint needs to have a model order m = modelOrder + offset
        // such that m > m(n) with m(n) being the model order of a normal node n (without constraints).
        // Such that the highest model order has to be used as an offset
        int offset = layeredGraph.getLayerlessNodes().size();
        for (LNode node : layeredGraph.getLayerlessNodes()) {
            if (node.hasProperty(InternalProperties.MODEL_ORDER)) {
                offset = Math.max(offset, node.getProperty(InternalProperties.MODEL_ORDER) + 1);
            }
        }
        
        while (true) {
            resetTARJAN(layeredGraph);
            TARJAN(layeredGraph);
            
            // If no Strongly connected components remain, the graph is acyclic.
            if (stronglyConnectedComponents.size() == 0) {
                break;
            }
            
            // highest model order only incoming
            findNodes(offset);
            
            // reverse the gathered edges
            for (LEdge edge : revEdges) {
                edge.reverse(layeredGraph, false);
                edge.getSource().getNode().setProperty(LayeredOptions.LAYERING_LAYER_ID,
                        edge.getSource().getNode().getProperty(LayeredOptions.LAYERING_LAYER_ID) + 1);
                layeredGraph.setProperty(InternalProperties.CYCLIC, true);
            }
            
            stronglyConnectedComponents.clear();
            revEdges.clear();
        }
        
        monitor.done();
        monitor.log("Execution Time: " + monitor.getExecutionTime());
    }
    
    public void findNodes(int offset) {
        for (int i = 0; i < stronglyConnectedComponents.size(); i++) {
            LNode max = null;
            int maxModelOrder = Integer.MIN_VALUE;
            for (LNode n : stronglyConnectedComponents.get(i)) {
                List<Integer> groupmask = new LinkedList<Integer>();
//                FIXME this may not always work and should not be hard coded.
//                int groupID = n.getProperty(LayeredOptions.CONSIDER_MODEL_ORDER_GROUP_I_D);
//                if (!groupmask.contains(groupID)) {
//                    continue;
//                }
                if (max == null) {
                    max = n;
                    maxModelOrder = computeConstraintModelOrder(n, offset);
                    continue;
                }
                int modelOrderCurrent = computeConstraintModelOrder(n, offset);
                if (maxModelOrder < modelOrderCurrent ) {
                    max = n;
                    maxModelOrder = modelOrderCurrent;
                }
            }
            for (LEdge edge : max.getOutgoingEdges()) {
                 if (stronglyConnectedComponents.get(i).contains(edge.getTarget().getNode())) {
                     revEdges.add(edge);
                 }
            
            }
        }
    }
    
    
    
    /**
     * Set model order to a value such that the constraint is respected and the ordering between nodes with
     * the same constraint is preserved.
     * The order should be FIRST_SEPARATE < FIRST < NORMAL < LAST < LAST_SEPARATE. The offset is used to make sure the 
     * all nodes have unique model orders.
     * @param node The LNode
     * @param offset The offset between FIRST, FIRST_SEPARATE, NORMAL, LAST_SEPARATE, and LAST nodes for unique order
     * @return A unique model order
     */
    protected int computeConstraintModelOrder(final LNode node, final int offset) {
        int modelOrder = 0;
        switch (node.getProperty(LayeredOptions.LAYERING_LAYER_CONSTRAINT)) {
        case FIRST_SEPARATE:
            modelOrder = 2 * -offset + firstSeparateModelOrder;
            firstSeparateModelOrder++;
            break;
        case FIRST:
            modelOrder = -offset;
            break;
        case LAST:
            modelOrder = offset;
            break;
        case LAST_SEPARATE:
            modelOrder = 2 * offset + lastSeparateModelOrder;
            lastSeparateModelOrder++;
            break;
        default:
            break;
        }
        if (node.hasProperty(InternalProperties.MODEL_ORDER)) {
            modelOrder += node.getProperty(InternalProperties.MODEL_ORDER);
        }
        return modelOrder;
    }

protected int computeConstraintGroupModelOrder(final LNode node, final int offset) {
        int modelOrder = 0;
        switch (node.getProperty(LayeredOptions.LAYERING_LAYER_CONSTRAINT)) {
        case FIRST_SEPARATE:
            modelOrder = 2 * -offset + firstSeparateModelOrder;
            firstSeparateModelOrder++;
            break;
        case FIRST:
            modelOrder = -offset;
            break;
        case LAST:
            modelOrder = offset;
            break;
        case LAST_SEPARATE:
            modelOrder = 2 * offset + lastSeparateModelOrder;
            lastSeparateModelOrder++;
            break;
        default:
            break;
        }
        if (node.hasProperty(LayeredOptions.CONSIDER_MODEL_ORDER_GROUP_I_D)) {
            modelOrder += node.getProperty(LayeredOptions.CONSIDER_MODEL_ORDER_GROUP_I_D);
        }
        return modelOrder;
    }

    private void TARJAN(final LGraph graph) {
        index = 0;
        stack = new Stack<LNode>();
        for (LNode node : graph.getLayerlessNodes()) {
            if (node.getProperty(InternalProperties.TARJAN_ID) == -1) {
                stronglyConnected(node);
                stack.clear();
            }
        }
    }
    
    private void stronglyConnected(final LNode v) {
        v.setProperty(InternalProperties.TARJAN_ID, index);
        v.setProperty(InternalProperties.TARJAN_LOWLINK, index);
        index++;
        stack.push(v);
        v.setProperty(InternalProperties.TARJAN_ON_STACK, true);
        for (LEdge edge : v.getConnectedEdges()) {
            if (edge.getSource().getNode() != v && !revEdges.contains(edge)) {
                continue;
            }
            if (edge.getSource().getNode() == v && revEdges.contains(edge)) {
                continue;
            }
            LNode target = null;
            if (edge.getTarget().getNode() == v) {
                target = edge.getSource().getNode();
            } else {
                target = edge.getTarget().getNode();
            }
            if (target.getProperty(InternalProperties.TARJAN_ID) == -1) {
                stronglyConnected(target);
                v.setProperty(InternalProperties.TARJAN_LOWLINK, 
                        Math.min(v.getProperty(InternalProperties.TARJAN_LOWLINK),
                        target.getProperty(InternalProperties.TARJAN_LOWLINK)));
            } else if (target.getProperty(InternalProperties.TARJAN_ON_STACK)) {
                v.setProperty(InternalProperties.TARJAN_LOWLINK, 
                        Math.min(v.getProperty(InternalProperties.TARJAN_LOWLINK),
                                target.getProperty(InternalProperties.TARJAN_ID)));
            }
        }
        if (v.getProperty(InternalProperties.TARJAN_LOWLINK) == v.getProperty(InternalProperties.TARJAN_ID)) {
            Set<LNode> sCC = new HashSet<LNode>();
            LNode n = null;
            do {
                n = stack.pop();
                n.setProperty(InternalProperties.TARJAN_ON_STACK, false);
                sCC.add(n);
            } while (v != n);
            if (sCC.size() > 1) {
                stronglyConnectedComponents.add(sCC);
            }
        }
    }
    
    private void resetTARJAN(final LGraph graph) {
        for (LNode n : graph.getLayerlessNodes()) {
            n.setProperty(InternalProperties.TARJAN_ON_STACK, false);
            n.setProperty(InternalProperties.TARJAN_LOWLINK, -1);
            n.setProperty(InternalProperties.TARJAN_ID, -1);
            stack.clear();
            for (LEdge e : n.getConnectedEdges()) {
                e.setProperty(InternalProperties.IS_PART_OF_CYCLE, false);
            }
        }
    }
}