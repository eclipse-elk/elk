/*******************************************************************************
 * Copyright (c) 2025 KIel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 * 
 * SPDX-License-Identifier: EPL-2.0 
 *******************************************************************************/
package org.eclipse.elk.alg.layered.p1cycles;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;
import java.util.TreeSet;

import org.eclipse.elk.alg.layered.LayeredPhases;
import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.intermediate.IntermediateProcessorStrategy;
import org.eclipse.elk.alg.layered.options.GroupOrderStrategy;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.core.alg.ILayoutPhase;
import org.eclipse.elk.core.alg.LayoutProcessorConfiguration;
import org.eclipse.elk.core.util.IElkProgressMonitor;

import com.google.common.collect.Iterables;

/**
 * Uses the Breadth-First-Search to traverse the graph and reverses edges if the node is already explored.
 * 
 * <p>This cycle breaker does not support the {@link LayeredOptions#PRIORITY_DIRECTION} option 
 * that can be set on edges. Neither does it support layer constraints out of the box. 
 * If layer constraints should be observed,
 * {@link org.eclipse.elk.alg.layered.intermediate.EdgeAndLayerConstraintEdgeReverser} and
 * {@link org.eclipse.elk.alg.layered.intermediate.LayerConstraintProcessor} should
 * be used.</p>
 * 
 * <dl>
 *   <dt>Precondition:</dt><dd>an unlayered graph</dd>
 *   <dt>Postcondition:</dt><dd>the graph has no cycles</dd>
 * </dl>
 * 
 * @see org.eclipse.elk.alg.layered.intermediate.EdgeAndLayerConstraintEdgeReverser
 * @see org.eclipse.elk.alg.layered.intermediate.LayerConstraintProcessor
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
    
    private LGraph graph;


    @Override
    public LayoutProcessorConfiguration<LayeredPhases, LGraph> getLayoutProcessorConfiguration(final LGraph graph) {
        return INTERMEDIATE_PROCESSING_CONFIGURATION;
    }

    @Override
    public void process(final LGraph graph, final IElkProgressMonitor monitor) {
        monitor.begin("Breadth-first cycle removal", 1);

        this.graph = graph;
        List<LNode> nodes = graph.getLayerlessNodes();

        // initialize values for the algorithm 
        bfsQueue = new LinkedList<LNode>();
        sources = new HashSet<>();
        sinks = new HashSet<>();
        visited = new boolean[nodes.size()];
        edgesToBeReversed = new ArrayList<>();
        
        // Find all sources and sinks in the graph.
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
        // This means each source may add their connections to the queue such that we search breadth-first.
        for (LNode source : sources) {

            //sequential bfs
            bfsQueue.add(source);
            bfsLoop();
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

    /**
     * Visits a node and adds its connections to the BF-queue.
     * @param n the node to visit
     */
    private void bfs(final LNode n) {
        // Return if the node was already visited.
        if (visited[n.id]) {
            return;
        }
        this.visited[n.id] = true;

        // Map to save the node model order of each edge connection.
        HashMap<Integer, HashSet<LEdge>> modelOrderMap = new HashMap<Integer, HashSet<LEdge>>();
        boolean groupModelOrder = this.graph.getProperty(LayeredOptions.CONSIDER_MODEL_ORDER_GROUP_MODEL_ORDER_CB_GROUP_ORDER_STRATEGY) == GroupOrderStrategy.ENFORCED;

        // Create a map of edges and the model order of the node they lead to
        for (LEdge e : n.getOutgoingEdges()) {
            if (!e.getTarget().getNode().hasProperty(InternalProperties.MODEL_ORDER)) {
                // Handle edges that connect to nodes without model order.
                // They get a high unique value such that the first such node is the last one and the second the second last.
                modelOrderMap.put(Integer.MAX_VALUE - modelOrderMap.size(), new HashSet<LEdge>(Arrays.asList(e)));
            } else {
                int targetModelOrder = 0;
                LNode target = e.getTarget().getNode();
                // Find out whether the model order group id or the model order is more important.
                if (groupModelOrder) {
                    // Get the biggest cycle breaking model order group. Now scale all groups such that
                    // maxModelOrderGroupSize * <model order group id> + model order creates a total ordering on all nodes.
                    // This orders all nodes without a group model order at the top.
                    // I leave this for know and find out whether this is desired. Maybe all need a group model order to begin with.
                    int maxModelOrderGroupSize = this.graph.getProperty(InternalProperties.MAX_MODEL_ORDER_NODES);
                    targetModelOrder = maxModelOrderGroupSize * target.getProperty(LayeredOptions.CONSIDER_MODEL_ORDER_GROUP_MODEL_ORDER_CYCLE_BREAKING_ID)
                            + target.getProperty(InternalProperties.MODEL_ORDER);
                } else {
                    targetModelOrder = e.getTarget().getNode().getProperty(InternalProperties.MODEL_ORDER);
                }
                // If the long edge target node has a model order, add it to the map.
                if (modelOrderMap.containsKey(targetModelOrder)){
                    modelOrderMap.get(targetModelOrder).add(e);
                } else {
                    modelOrderMap.put(targetModelOrder, new HashSet<LEdge>(Arrays.asList(e)));
                }
            }
        }
        // This holds all model orders of nodes connected to the current node sorted by model order.
        // Basically this orders all different model orders (or group model orders) by priority.
        TreeSet<Integer> modelOrderSet = new TreeSet<>(modelOrderMap.keySet());

        // Since the model order determines the iteration order of e
        for (int key : modelOrderSet) {
            LEdge out = modelOrderMap.get(key).iterator().next();
            // Do not visit self loops
            if(out.isSelfLoop()) {
                continue;
            }
            // If the target was already visited, reverse the edge to it.
            LNode target = out.getTarget().getNode();
            // 
            if (this.visited[target.id] && !sources.contains(n) && !sinks.contains(target)) {
                edgesToBeReversed.addAll(modelOrderMap.get(key));
            } else {
                bfsQueue.add(target);
            }
        }
    }
}