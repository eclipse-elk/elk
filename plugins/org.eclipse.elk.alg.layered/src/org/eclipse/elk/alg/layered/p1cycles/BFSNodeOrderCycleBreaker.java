/*******************************************************************************
 * Copyright (c) 2024 Sasuk and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 * 
 * SPDX-License-Identifier: EPL-2.0 
 *******************************************************************************/
package org.eclipse.elk.alg.layered.p1cycles;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;
import java.util.Set;
import java.util.SortedSet;
import java.util.Stack;
import java.util.TreeSet;

import org.eclipse.elk.alg.layered.DebugUtil;
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

import com.google.common.collect.Iterables;

/**
 * Uses the Breadth-First-Search to traverse the graph and reverses edges if the node is already explored.
 * This implementation uses Tarjan's Algorithm to reduce unnecessary reversals.
 * @author Mwr
 *
 */
public class BFSNodeOrderCycleBreaker implements ILayoutPhase<LayeredPhases, LGraph> {
    
    /** Intermediate processing configuration. */
    private static final LayoutProcessorConfiguration<LayeredPhases, LGraph> INTERMEDIATE_PROCESSING_CONFIGURATION =
        LayoutProcessorConfiguration.<LayeredPhases, LGraph>create()
            .addAfter(LayeredPhases.P5_EDGE_ROUTING, IntermediateProcessorStrategy.REVERSED_EDGE_RESTORER);

    /** Set of source nodes. */
    private HashSet<LNode> sources;
    
    /** Set of sink nodes. */
    private HashSet<LNode> sinks;
    
    /** Indicates whether a given node was already visited during BFS. */
    private boolean[] visited;
    
    /**
     * Queues the nodes for BFS.
     */
    private Queue<LNode> bfsQueue;
    
    /** The list of edges to be reversed at the end of our little algorithmic adventure. */
    private List<LEdge> edgesToBeReversed;
    
    
    private int index;
    protected List<Set<LNode>> stronglyConnectedComponents;
    private Stack<LNode> stack;
    private HashMap <LNode,Integer> nodeToSCCID;

    
    @Override
    public LayoutProcessorConfiguration<LayeredPhases, LGraph> getLayoutProcessorConfiguration(final LGraph graph) {
        return INTERMEDIATE_PROCESSING_CONFIGURATION;
    }
    
    @Override
    public void process(final LGraph graph, final IElkProgressMonitor monitor) {
        monitor.begin("Breadth-first cycle removal", 1);
        
        List<LNode> nodes = graph.getLayerlessNodes();

        // initialize values for the algorithm 
        bfsQueue = new LinkedList<LNode>();
        sources = new HashSet<>();
        sinks = new HashSet<>();
        visited = new boolean[nodes.size()];
        edgesToBeReversed = new ArrayList<>();
        stronglyConnectedComponents = new LinkedList<Set<LNode>>();
        stack = new Stack<LNode>();
        nodeToSCCID = new HashMap<>();
        
        TARJAN(graph);
        if (stronglyConnectedComponents.size() == 0) {
            // Cleanup
            this.sources = null;
            this.visited = null;
            this.bfsQueue = null;
            this.edgesToBeReversed = null;
            
            monitor.done();
            return;
        }
        int index = 0;
        for (LNode node : nodes) {
            // The node id is used as index into our arrays
            node.id = index;
            if (Iterables.isEmpty(node.getIncomingEdges())) {
                sources.add(node);
            }
            if (Iterables.isEmpty(node.getOutgoingEdges())) {
                sinks.add(node);
            }
            index++;
        }
        
        // Start BFS Search starting at each source sequentially.
        for (LNode source : sources) {
            
            //sequential bfs
            bfsQueue.add(source);
            bfsLoop();
            
            
            //Simultaneous starting bfs
            //bfs(source);
        }
        
        bfsLoop();
        
        // Start more BFS runs from the first node that has not been visited yet. This must be part of a cycle since it
        // is not a source nodes
        boolean changed = true;
        while(changed) {
            changed = false;
            for (int i = 0; i < nodes.size(); i++) {
                if (!visited[i]) {
                    LNode n = nodes.get(i);
                    assert n.id == i;
                    bfsQueue.add(n);
                    changed = true;
                    break;
                }
            }
            bfsLoop();
        }
        
        
        // Reverse "back edges"
        for (LEdge edge : edgesToBeReversed) {
            edge.reverse(graph, true);
            graph.setProperty(InternalProperties.CYCLIC, true);
        }
        
        // Cleanup
        this.sources = null;
        this.visited = null;
        this.bfsQueue = null;
        this.edgesToBeReversed = null;
        
        monitor.done();
    }
    
    private void bfsLoop() {
        while(!bfsQueue.isEmpty()) {
            bfs(bfsQueue.poll());
        }
    }
    
    private void bfs(final LNode n) {
        if (visited[n.id]) {
            return;
        }
        this.visited[n.id] = true;
        
        HashMap<Integer,HashSet<LEdge>> modelOrderMap = new HashMap<Integer,HashSet<LEdge>>();
        
        //Create a map of edges and the model order of the node they lead to
        n.getOutgoingEdges().forEach(e -> {
            if (e.getTarget().getNode().getProperty(InternalProperties.MODEL_ORDER) == null) {
                modelOrderMap.put(Integer.MAX_VALUE - modelOrderMap.size(),new HashSet<LEdge>(){{add(e);}});
            }
            else {
                int targetModelOrder = e.getTarget().getNode().getProperty(InternalProperties.MODEL_ORDER);
                if (modelOrderMap.containsKey(targetModelOrder)){
                    modelOrderMap.get(targetModelOrder).add(e);
                }
                else {
                    modelOrderMap.put(targetModelOrder, new HashSet<LEdge>(){{add(e);}});
                }
            }
        });
        SortedSet<Integer> modelOrderSet = new TreeSet<>(modelOrderMap.keySet());
        
        for (int key: modelOrderSet) {
            LEdge out = modelOrderMap.get(key).iterator().next();
         // for (LEdge out : n.getOutgoingEdges()) {
            // Get the SCC of the source, or -1 if not part of a SCC
            int source_SCC = nodeToSCCID.get(n) != null ?  nodeToSCCID.get(n) : -1;
            if(out.isSelfLoop()) {
                continue;
            }
            LNode target = out.getTarget().getNode();
            // Get the SCC of the target, or -2 (force difference between source and target) if not part of a SCC
            int target_SCC = nodeToSCCID.get(target) != null ? nodeToSCCID.get(target) : -2;
            if (this.visited[target.id] && !sources.contains(n) && !sinks.contains(target)
                    //&& source_SCC == target_SCC
                    ) {
                edgesToBeReversed.addAll(modelOrderMap.get(key));
            }else {
                bfsQueue.add(target);
            }
        }
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
            if (edge.getSource().getNode() != v && !edgesToBeReversed.contains(edge)) {
                continue;
            }
            if (edge.getSource().getNode() == v && edgesToBeReversed.contains(edge)) {
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
            if (sCC.size() >1) {
                int index = stronglyConnectedComponents.size();
                stronglyConnectedComponents.add(sCC);
                for (LNode node : sCC) {
                    nodeToSCCID.put(node, index);
                }
            }
        }
    }
}

