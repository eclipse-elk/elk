/*******************************************************************************
 * Copyright (c) 2013 - 2022 Kiel University and others.
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
package org.eclipse.elk.alg.mrtree.p1treeify;

import java.util.LinkedList;
import java.util.List;

import org.eclipse.elk.alg.mrtree.TreeLayoutPhases;
import org.eclipse.elk.alg.mrtree.graph.TEdge;
import org.eclipse.elk.alg.mrtree.graph.TGraph;
import org.eclipse.elk.alg.mrtree.graph.TNode;
import org.eclipse.elk.alg.mrtree.intermediate.IntermediateProcessorStrategy;
import org.eclipse.elk.alg.mrtree.options.InternalProperties;
import org.eclipse.elk.alg.mrtree.options.MrTreeOptions;
import org.eclipse.elk.alg.mrtree.options.TreeifyingOrder;
import org.eclipse.elk.core.alg.ILayoutPhase;
import org.eclipse.elk.core.alg.LayoutProcessorConfiguration;
import org.eclipse.elk.core.util.IElkProgressMonitor;

/**
 * This phase should create a tree out of a graph that is hardly a tree.
 * It should for example look for cycles in a tree and remove that cycles. This cycle breaking is
 * done through eliminating the edge that destroys the tree property and putting that edge into a
 * list. So this phase builds up a list with removed edges, so that the Mr. Tree algorithm can
 * operate on the newly constructed graph which is now a tree. The previously removed edges could be
 * reinserted during a post-processing. The search order can be DFS or BFS.
 * 
 * @author sor
 * @author sgu
 * @author msp
 * @author sdo
 */
public class DFSTreeifyer implements ILayoutPhase<TreeLayoutPhases, TGraph> {

    /** intermediate processing configuration. */
    private static final LayoutProcessorConfiguration<TreeLayoutPhases, TGraph> INTERMEDIATE_PROCESSING_CONFIG =
            LayoutProcessorConfiguration.<TreeLayoutPhases, TGraph>create()
                    .addAfter(TreeLayoutPhases.P3_NODE_PLACEMENT, IntermediateProcessorStrategy.DETREEIFYING_PROC);

    /**
     * {@inheritDoc}
     */
    public void process(final TGraph tGraph, final IElkProgressMonitor progressMonitor) {
        progressMonitor.begin("DFS Treeifying phase", 1);

        init(tGraph);
        collectEdges(tGraph);

        eliminated = null;
        visited = null;
        progressMonitor.done();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public LayoutProcessorConfiguration<TreeLayoutPhases, TGraph> getLayoutProcessorConfiguration(final TGraph graph) {
        return INTERMEDIATE_PROCESSING_CONFIG;
    }

    /** marks the nodes that have been visited during DFS. */
    private int[] visited;

    /** a list with the edges that can be removed. */
    private List<TEdge> eliminated;

    /**
     * Initializes all previously defined structures of a TGraph with useful values.
     * 
     * @param tGraph
     */
    private void init(final TGraph tGraph) {

        int size = tGraph.getNodes().size();
        eliminated = new LinkedList<TEdge>();
        visited = new int[size];

        /** for using id property initialize ids first */
        int id = 0;
        for (TNode node : tGraph.getNodes()) {
            node.id = id++;
        }
    }

    /**
     * Create a list with the edges to remove. These edges should be all edges that destroy the
     * tree-property. 
     * 
     * @param TGraph
     *            where to collect the edges
     */
    private void collectEdges(final TGraph tGraph) {
        TreeifyingOrder treeifyingOrder = tGraph.getProperty(MrTreeOptions.SEARCH_ORDER);

        // start DFS on every node in graph
        for (TNode tNode : tGraph.getNodes()) {
            // if the node has not been visited yet
            if (visited[tNode.id] == 0) {
                // call DFS or BFS on that node
                switch (treeifyingOrder) {
                case DFS:
                    dfs(tNode);
                    break;
                case BFS:
                    bfs(tNode);
                    break;
                }
                // if we come back to that node again, set the node as root
                visited[tNode.id] = 2;
            }
        }
        // remove the found edges out of graph structure
        for (TEdge tEdge : eliminated) {
            tEdge.getSource().getOutgoingEdges().remove(tEdge);
            tEdge.getTarget().getIncomingEdges().remove(tEdge);
        }
        // set the list of collected edges as a graph property
        tGraph.setProperty(InternalProperties.REMOVABLE_EDGES, eliminated);
    }

    /**
     * Perform a DFS on a given graph till all nodes of the graph have been visited.
     * 
     * @param node
     *            to start DFS
     */
    private void dfs(final TNode tNode) {
        // dfs starts on a node and marks that node as visited
        visited[tNode.id] = 1;

        // go to all child nodes of start node
        for (TEdge tEdge : tNode.getOutgoingEdges()) {
            TNode target = tEdge.getTarget();
            // if a child has been visited
            if (visited[target.id] == 1) {
                // put that edge to the list that contains the edges to remove
                eliminated.add(tEdge);
            } else if (visited[target.id] == 2) {
                // if a previous root can be visited from another node unmark the root property
                visited[target.id] = 1;
            } else {
                // recurse
                dfs(target);
            }
        }
    }
    
    /**
     * Perform a BFS on a given graph till all nodes of the graph have been visited.
     * 
     * @param node
     *            to start BFS
     */
    private void bfs(final TNode startNode) {
        LinkedList<TNode> nodeQueue = new LinkedList<TNode>();
        nodeQueue.add(startNode);
        
        do {
            TNode node = nodeQueue.removeFirst();
            visited[node.id] = 1;

            for (TEdge tEdge : node.getOutgoingEdges()) {
                TNode target = tEdge.getTarget();
                // if a child has been visited
                if (visited[target.id] == 1) {
                    // put that edge to the list that contains the edges to remove
                    eliminated.add(tEdge);
                } else if (visited[target.id] == 2) {
                    // if a previous root can be visited from another node unmark the root property
                    visited[target.id] = 1;
                } else {
                    // process the node in breadth-first order
                    nodeQueue.addLast(target);
                }
            }
        } while (!nodeQueue.isEmpty());
    }
    
}
