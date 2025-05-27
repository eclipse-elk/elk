/*******************************************************************************
 * Copyright (c) 2023 Sasuk and others.
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
import java.util.Queue;
import java.util.Set;
import java.util.Stack;
import java.util.concurrent.ConcurrentLinkedQueue;

import org.eclipse.elk.alg.layered.LayeredPhases;
import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.intermediate.IntermediateProcessorStrategy;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.layered.options.PortType;
import org.eclipse.elk.core.alg.ILayoutPhase;
import org.eclipse.elk.core.alg.LayoutProcessorConfiguration;
import org.eclipse.elk.core.util.IElkProgressMonitor;

import com.google.common.collect.Lists;

/**
 * Cycle Breaking Strategy that looks further down the graph for nodes of the same type to reverse edges.
 * Only works with Group Model Order.
 * @author Mwr
 *
 */
public class ModelOrderLookAHeadCycleBreaker implements ILayoutPhase<LayeredPhases, LGraph>{
    private int firstSeparateModelOrder;
    private int lastSeparateModelOrder;
    private Queue<LNode> nextNodes;
    private int index;
    private LinkedList<Set<LNode>> stronglyConnectedComponents = new LinkedList<Set<LNode>>();
    private Stack<LNode> stack = new Stack<LNode>();
    private List<LEdge> revEdges = Lists.newArrayList();

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
        
        // gather edges that point to the wrong direction.
        revEdges = Lists.newArrayList();
        
        int currentTopMO = 0;
        for (LNode node : layeredGraph.getLayerlessNodes()) {
            if (node.hasProperty(InternalProperties.MODEL_ORDER)) {
                currentTopMO = Math.max(currentTopMO, node.getProperty(InternalProperties.MODEL_ORDER) + 1);
            }
        }
        //add fake "Model Order" so nodes without one can be dealt with.
        for (LNode node : layeredGraph.getLayerlessNodes()) {
            if (node.getProperty(InternalProperties.MODEL_ORDER) == null 
                    && !node.toString().contains("external_port")) {
                node.setProperty(InternalProperties.MODEL_ORDER, currentTopMO);
                currentTopMO++;
            }
        }
        
        // create strongly connected components.
        TARJAN(layeredGraph);
        
        // mark edges that are part of a cycle (within a SCC).
        for (int i = 0; i < stronglyConnectedComponents.size(); i++) {
            for (LNode n : stronglyConnectedComponents.get(i)) {
             for (LEdge edge : n.getOutgoingEdges()) {
                 if (stronglyConnectedComponents.get(i).contains(edge.getTarget().getNode())) {
                     edge.setProperty(InternalProperties.IS_PART_OF_CYCLE, true);
                 }
             }
            }
        }
        stronglyConnectedComponents.clear();
        
        
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
        lookAhead(layeredGraph, offset);
        
        
        
        // reverse the gathered edges
        for (LEdge edge : revEdges) {
            edge.reverse(layeredGraph, false);
            layeredGraph.setProperty(InternalProperties.CYCLIC, true);
        }
        revEdges.clear();
        monitor.done();
    }
    
    /**
     * Looks Further down the graph and uses the skip sequential edges option with fallback edges.
     * @param layeredGraph the graph
     * @param offset the offset for working with Model Order
     */
    private void lookAhead(final LGraph layeredGraph,final int offset) {
        int modelOrderMask[] = {1,4};
        for (int i: modelOrderMask) {
            for (LNode source : layeredGraph.getLayerlessNodes()) {
                int groupID = source.getProperty(LayeredOptions.CONSIDER_MODEL_ORDER_GROUP_I_D);
                if (groupID != i) {    
                    continue;
                }
                int modelOrderSource = computeConstraintModelOrder(source, offset);
                
                
                for (LPort port : source.getPorts(PortType.OUTPUT)) {
                    for (LEdge edge : port.getOutgoingEdges()) {
                        if (revEdges.contains(edge) || !edge.getProperty(InternalProperties.IS_PART_OF_CYCLE)) {
                            continue;
                        }
                        LNode target = edge.getTarget().getNode();
                        int groupIDTarget = target.getProperty(LayeredOptions.CONSIDER_MODEL_ORDER_GROUP_I_D);
                        if (groupID == -1 || groupIDTarget == -1) {
                            continue;
                        }
                        if (groupIDTarget != groupID) {
                            nextNodes = new ConcurrentLinkedQueue<LNode>();
                            HashSet<LNode> visited = new HashSet<LNode>();
                            HashSet<LEdge> toRev = new HashSet<LEdge>();
                            visited.add(target);
                            toRev = getNextNodeWithGroupPriority(target, groupID, 
                                    offset, modelOrderSource, visited, source, edge, toRev);
                            revEdges.addAll(toRev);
                            
                        } else if (computeConstraintModelOrder(target, offset) <= modelOrderSource) {
                            revEdges.add(edge);
                        } else {
                            edge.setProperty(InternalProperties.IS_FIXED, true);
                        }
                        
                    }
                }
            }
        }
    }
    
    /**
     * Search for the next Node with the same groupID of the source Node, using BFS.
     * @param currentSource the current source node
     * @param groupPriority the groupID of the source
     * @param offset current offset of the modelOrder
     * @param modelOrderOriginalSource the original source node
     * @return true if the edge should be reversed
     */
    private HashSet<LEdge> getNextNodeWithGroupPriority(final LNode currentSource, 
            final int groupPriority, final int offset, 
            final int modelOrderOriginalSource, final HashSet<LNode> checked, final LNode originalSource, 
            final LEdge fallBackEdge, final HashSet<LEdge> toRev) {
        for (LPort port : currentSource.getPorts(PortType.OUTPUT)) {
            for (LEdge edge : port.getOutgoingEdges()) {
                if  (revEdges.contains(edge) || toRev.contains(edge) 
                        || !edge.getProperty(InternalProperties.IS_PART_OF_CYCLE)) {
                    continue;
                }
                LNode target = edge.getTarget().getNode();
                int groupPriorityTarget = target.getProperty(LayeredOptions.CONSIDER_MODEL_ORDER_GROUP_I_D);
                if (groupPriorityTarget == groupPriority) {
                    if (computeConstraintModelOrder(target, offset) <= modelOrderOriginalSource) {
                        if (edge.getProperty(InternalProperties.IS_FIXED)) {
                            toRev.clear();
                            toRev.add(fallBackEdge);
                            return toRev;
                        }
                        toRev.add(edge);
                    } else {
                        edge.setProperty(InternalProperties.IS_FIXED, true);
                    }
                } else if (!checked.contains(target)) {
                    nextNodes.add(target);
                    checked.add(target);
                }
                
            }
        }   
        if (!nextNodes.isEmpty()) {
            getNextNodeWithGroupPriority(nextNodes.poll(), groupPriority, offset, modelOrderOriginalSource, checked,
                    originalSource, fallBackEdge, toRev);
        }
        return toRev;
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
    private int computeConstraintModelOrder(final LNode node, final int offset) {
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
        for (LEdge edge : v.getOutgoingEdges()) {
            LNode target = edge.getTarget().getNode();
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
            stronglyConnectedComponents.add(sCC);
        }
    }
}