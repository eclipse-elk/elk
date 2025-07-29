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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
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
 * Cycle breaker implementation that uses a depth-first traversal of the graph. Described in
 * <ul>
 * <li>Emden R. Gansner, Eleftherios Koutsofios, Stephen C. North, Kiem-Phong Vo, A technique for drawing directed
 * graphs. <i>Software Engineering</i> 19(3), pp. 214-230, 1993.</li>
 * </ul>
 * While the {@link org.eclipse.elk.alg.layered.p1cycles.DepthFirstCycleBreaker} uses the edge order as the visiting order
 * this cycle breaker uses the node (model order) for this.
 * 
 * <p>
 * This cycle breaker does not support the {@link LayeredOptions#PRIORITY_DIRECTION} option that can be set on edges.
 * Neither does it support layer constraints out of the box. If layer constraints should be observed,
 * {@link org.eclipse.elk.alg.layered.intermediate.EdgeAndLayerConstraintEdgeReverser} and
 * {@link org.eclipse.elk.alg.layered.intermediate.LayerConstraintProcessor} should be used.
 * </p>
 * 
 * <dl>
 * <dt>Precondition:</dt>
 * <dd>an unlayered graph</dd>
 * <dt>Postcondition:</dt>
 * <dd>the graph has no cycles</dd>
 * </dl>
 * 
 * @see org.eclipse.elk.alg.layered.intermediate.EdgeAndLayerConstraintEdgeReverser
 * @see org.eclipse.elk.alg.layered.intermediate.LayerConstraintProcessor
 */
public class DFSNodeOrderCycleBreaker implements ILayoutPhase<LayeredPhases, LGraph> {

    /** Intermediate processing configuration. */
    private static final LayoutProcessorConfiguration<LayeredPhases, LGraph> INTERMEDIATE_PROCESSING_CONFIGURATION =
            LayoutProcessorConfiguration.<LayeredPhases, LGraph> create().addAfter(LayeredPhases.P5_EDGE_ROUTING,
                    IntermediateProcessorStrategy.REVERSED_EDGE_RESTORER);

    /** List of source nodes. */
    private List<LNode> sources;
    /** Indicates whether a given node was already visited during DFS. */
    private boolean[] visited;
    /**
     * A node is active during DFS if it is on our current DFS path. Any edge that leads back to an active node induces
     * a cycle and needs to be reversed.
     */
    private boolean[] active;
    /** The list of edges to be reversed at the end of our little algorithmic adventure. */
    private List<LEdge> edgesToBeReversed;
    
    private LGraph graph;

    @Override
    public LayoutProcessorConfiguration<LayeredPhases, LGraph> getLayoutProcessorConfiguration(final LGraph graph) {
        return INTERMEDIATE_PROCESSING_CONFIGURATION;
    }

    @Override
    public void process(final LGraph graph, final IElkProgressMonitor monitor) {
        monitor.begin("Depth-first cycle removal", 1);
        
        this.graph = graph;
        List<LNode> nodes = graph.getLayerlessNodes();

        // initialize values for the algorithm
        int nodeCount = nodes.size();

        sources = new ArrayList<>();
        visited = new boolean[nodeCount];
        active = new boolean[nodeCount];
        edgesToBeReversed = new ArrayList<>();

        int index = 0;
        for (LNode node : nodes) {
            // The node id is used as index into our arrays
            node.id = index;
            if (Iterables.isEmpty(node.getIncomingEdges())) {
                sources.add(node);
            }
            index++;
        }

        // From every source node start a DFS
        for (LNode source : sources) {
            dfs(source);
        }

        // Start more DFS runs for all nodes that have not been visited yet. These must be part of a cycle since they
        // are not source nodes
        for (int i = 0; i < nodeCount; i++) {
            if (!visited[i]) {
                LNode n = nodes.get(i);
                assert n.id == i;
                dfs(n);
            }
        }

        // Reverse "back edges"
        for (LEdge edge : edgesToBeReversed) {
            edge.reverse(graph, true);
            graph.setProperty(InternalProperties.CYCLIC, true);
        }

        // Cleanup
        this.sources = null;
        this.visited = null;
        this.active = null;
        this.edgesToBeReversed = null;


        monitor.done();
    }

    private void dfs(final LNode n) {
        if (visited[n.id]) {
            return;
        }

        // We're now visiting the node, and it's active
        this.visited[n.id] = true;
        this.active[n.id] = true;

        HashMap<Integer, HashSet<LEdge>> modelOrderMap = new HashMap<Integer, HashSet<LEdge>>();
        boolean groupModelOrder = this.graph.getProperty(LayeredOptions.CONSIDER_MODEL_ORDER_GROUP_MODEL_ORDER_CB_GROUP_ORDER_STRATEGY) == GroupOrderStrategy.ENFORCED;

        
        // Construct tree set to efficiently find next node to look at by finding all connected nodes and sorting them.
        n.getOutgoingEdges().forEach(e -> {
            if (e.getTarget().getNode().getProperty(InternalProperties.MODEL_ORDER) == null) {
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
                if (modelOrderMap.containsKey(targetModelOrder)) {
                    modelOrderMap.get(targetModelOrder).add(e);
                } else {
                    modelOrderMap.put(targetModelOrder, new HashSet<LEdge>(Arrays.asList(e)));
                }
            }
        });
        // The model order serves as the key by which all elements are sorted.
        // Construct this by <number of nodes> * groupId + MO if ENFORCED, otherwise using model order is ok.
        TreeSet<Integer> modelOrderSet = new TreeSet<>(modelOrderMap.keySet());
        for (int key : modelOrderSet) {
            LEdge out = modelOrderMap.get(key).iterator().next();
            // Ignore self loops
            if (out.isSelfLoop()) {
                continue;
            }

            LNode target = out.getTarget().getNode();

            // If the edge connects to an active node, we have found a path from said active node back to itself since
            // active nodes are on our current path. That's a backward edge and needs to be reversed
            if (this.active[target.id]) {
                edgesToBeReversed.addAll(modelOrderMap.get(key));
            } else {
                dfs(target);
            }
        }

        // We're leaving this node
        this.active[n.id] = false;
    }

}