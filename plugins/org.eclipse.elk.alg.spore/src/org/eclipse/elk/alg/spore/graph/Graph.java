/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.spore.graph;

import java.util.List;
import java.util.Set;

import org.eclipse.elk.alg.common.ICostFunction;
import org.eclipse.elk.alg.common.TEdge;
import org.eclipse.elk.alg.common.Tree;
import org.eclipse.elk.alg.common.spore.Node;
import org.eclipse.elk.alg.spore.options.CompactionStrategy;
import org.eclipse.elk.alg.spore.options.TreeConstructionStrategy;
import org.eclipse.elk.graph.properties.MapPropertyHolder;

import com.google.common.collect.Lists;

/**
 * A graph structure for SPOrE.
 */
public class Graph extends MapPropertyHolder {
    /** The serial version UID. */
    private static final long serialVersionUID = -8005862731704572621L;
    
    // CHECKSTYLEOFF VisibilityModifier
    /** All vertices of this graph. */
    public List<Node> vertices = Lists.newArrayList();
    /** All edges of this graph. */
    public Set<TEdge> tEdges;
    /** Determines the kind of spanning tree to be generated. */
    public final TreeConstructionStrategy treeConstructionStrategy;
    /** Holds the tree structure for the processing order. */
    public Tree<Node> tree;
    /** Determines the compaction method applied to the spanning tree. */
    public final CompactionStrategy compactionStrategy;
    /** Callback function for the cost/weight of an edge. */
    public final ICostFunction costFunction;
    /** One of the vertices can be selected to be used as the root for the tree construction. */
    public Node preferredRoot;
    /** Restricts the translation of nodes to orthogonal directions in a compaction. */
    public boolean orthogonalCompaction;
    // CHECKSTYLEON VisibilityModifier
    
    /**
     * Constructor. The callback functions and the strategy for tree construction that are needed in the processing
     * of this graph are set at this point.
     * @param costFun a function returning a cost value for a {@link TEdge}
     * @param treeStrategy strategy for the spanning tree construction (e.g. MINIMUM_SPANNING_TREE)
     * @param compStrategy a strategy for the compaction
     */
    public Graph(final ICostFunction costFun, final TreeConstructionStrategy treeStrategy, 
            final CompactionStrategy compStrategy) {
        costFunction = costFun;
        treeConstructionStrategy = treeStrategy;
        compactionStrategy = compStrategy;
    }
}
