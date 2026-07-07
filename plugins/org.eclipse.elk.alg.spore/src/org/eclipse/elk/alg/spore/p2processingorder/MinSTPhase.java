/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.spore.p2processingorder;

import java.util.Map;

import org.eclipse.elk.alg.common.NaiveMinST;
import org.eclipse.elk.alg.common.Tree;
import org.eclipse.elk.alg.common.spore.InternalProperties;
import org.eclipse.elk.alg.common.spore.Node;
import org.eclipse.elk.alg.spore.SPOrEPhases;
import org.eclipse.elk.alg.spore.graph.Graph;
import org.eclipse.elk.core.alg.ILayoutPhase;
import org.eclipse.elk.core.alg.LayoutProcessorConfiguration;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.util.ElkUtil;
import org.eclipse.elk.core.util.IElkProgressMonitor;

import com.google.common.collect.Maps;

/**
 * This phase finds a minimum spanning tree for the edges of the {@link Graph} with the 
 * specified root {@link Graph#preferredRoot}.
 */
public class MinSTPhase implements ILayoutPhase<SPOrEPhases, Graph> {
    
    private Map<KVector, Node> nodeMap = Maps.newHashMap();

    @Override
    public LayoutProcessorConfiguration<SPOrEPhases, Graph> getLayoutProcessorConfiguration(final Graph graph) {
        return LayoutProcessorConfiguration.<SPOrEPhases, Graph>create();
    }
    
    @Override
    public void process(final Graph graph, final IElkProgressMonitor progressMonitor) {
        progressMonitor.begin("Minimum spanning tree construction", 1);
        
        KVector root;
        if (graph.preferredRoot != null) {
            root = graph.preferredRoot.originalVertex;
        } else {
            root = graph.vertices.get(0).originalVertex;
        }
        
        Tree<KVector> tTree;
        if (graph.getProperty(InternalProperties.DEBUG_SVG)) {
            tTree = NaiveMinST.createSpanningTree(graph.tEdges, root, graph.costFunction, 
                    ElkUtil.debugFolderPath("spore") + "20minst");
        } else {
            tTree = NaiveMinST.createSpanningTree(graph.tEdges, root, graph.costFunction);
        }
        
        // convert result to a Tree that can be used in the execution phase
        convert(tTree, graph);
        
        progressMonitor.done();
    }
    
    /**
     * Converts a {@link Tree} of {@link KVector}s to a {@link Tree} of {@link Node}s.
     * <p>Can be used in extending class.</p>
     * @param tTree a tree of vectors
     * @param graph the {@link Graph} that contains the {@link Tree} of {@link Node}s
     */
    protected void convert(final Tree<KVector> tTree, final Graph graph) {
        nodeMap.clear();
        graph.vertices.forEach(n -> nodeMap.put(n.originalVertex, n));
        Tree<Node> root = new Tree<Node>(nodeMap.get(tTree.node));
        addNode(root, tTree);
        graph.tree = root;
    }
    
    /**
     * Recursive function to find the node related to a vertex and add it to an existing tree.
     * @param s
     * @param t
     */
    private void addNode(final Tree<Node> s, final Tree<KVector> t) {
        for (Tree<KVector> tTree : t.children) {
            Tree<Node> child = new Tree<Node>(nodeMap.get(tTree.node));
            s.children.add(child);
            addNode(child, tTree);
        }
    }
}
