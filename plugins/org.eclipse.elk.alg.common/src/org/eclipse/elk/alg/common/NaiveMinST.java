/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.common;

import java.util.List;
import java.util.Map;
import java.util.Set;

import org.eclipse.elk.alg.common.utils.SVGImage;
import org.eclipse.elk.core.math.KVector;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;

/**
 * Simple minimum spanning tree calculation.
 * This class expects a list of edges that are pairs of vertices, a cost function for these edges, and a
 * root vertex to start the spanning tree. 
 * <p>precondition: The edge list E represents a connected graph (V,E) where V={w|(u,w) \in E or (w,v) \in E}
 * and root \in V.</p>
 * <p>postcondition: The returned tree structure connects all vertices in V once.</p>
 */
public final class NaiveMinST {
    
    /** Hidden constructor. */
    private NaiveMinST() { };
    
    /**
     * Creates a minimum spanning tree for a graph using the cost function specified in the constructor.
     * @param tEdges the edges of the graph
     * @param root the root node to start the spanning tree
     * @param costFunction a function returning a cost value for a {@link TEdge}
     * @param debugOutputFile file name for debug SVG. Debug output will be deactivated if this is null.
     * @return the spanning tree
     */
    public static Tree<KVector> createSpanningTree(final Set<TEdge> tEdges, final KVector root, 
            final ICostFunction costFunction, final String debugOutputFile) {
        // determine edge weights
        Map<TEdge, Double> weight = Maps.newHashMap();
        for (TEdge edge : tEdges) {
            weight.put(edge, costFunction.cost(edge));
        }
        
        // sort edges by weight
        List<TEdge> edgeList = Lists.newArrayList(tEdges);
        edgeList.sort((TEdge e1, TEdge e2) -> weight.get(e1).compareTo(weight.get(e2)));
        Set<TEdge> edges = Sets.newLinkedHashSet(edgeList); // preserves order
        
        
        // iteratively add cheapest edge where one node is contained in current tree and one is new
        Tree<KVector> minST = new Tree<KVector>(root);
        Map<KVector, Tree<KVector>> treeNodes = Maps.newHashMap();
        treeNodes.put(root, minST);
        
        // debug output ----------------------------------------------------------------------------------------------
        SVGImage svg = new SVGImage(debugOutputFile);
        // elkjs-exclude-start
        svg.addGroups("e", "t");
        for (TEdge e : edges) {
            svg.g("e").addLine(e.u.x, e.u.y, e.v.x, e.v.y, "stroke=\"black\" stroke-width=\"1\"");
            svg.g("t").addElementStr("<text x=\"" + (e.u.x + e.v.x) / 2 
                + "\" y=\"" + (e.u.y + e.v.y) / 2 + "\" fill=\"blue\""
                + " font-size=\"20px\">" + String.format("%.2f", weight.get(e)) + "</text>");
            }
        svg.isave();
        // elkjs-exclude-end
        // -----------------------------------------------------------------------------------------------------------
        
        while (!edges.isEmpty()) {
            TEdge nextEdge = null;
            KVector nextNode = null;
            KVector nodeInTree = null;
            double minWeight = Double.POSITIVE_INFINITY;
            for (TEdge edge : edges) {
                if (weight.get(edge) <= minWeight) {
                    if (treeNodes.containsKey(edge.u) && !treeNodes.containsKey(edge.v)) {
                        nextNode = edge.v;
                        nodeInTree = edge.u;
                        nextEdge = edge;
                        break; // because edges is sorted we don't have to look any further
                    }
                    if (treeNodes.containsKey(edge.v)) {
                        if (!treeNodes.containsKey(edge.u)) {
                            nextNode = edge.u;
                            nodeInTree = edge.v;
                            nextEdge = edge;
                            break; // because edges is sorted we don't have to look any further
                        }
                    }
                }
            }
            
            // Since the underlying graph is connected, the spanning tree is complete if no edge can be found that
            // connects a new node to the tree.
            if (nextEdge == null) {
                break;
            }
            
            // add the new node to the spanning tree
            Tree<KVector> subTree = new Tree<KVector>(nextNode);
            treeNodes.get(nodeInTree).children.add(subTree);
            treeNodes.put(nextNode, subTree);
            edges.remove(nextEdge);
            
            // debug output -------------------------------------------------------------------------------------------
            svg.g("e").addLine(nextEdge.u.x, nextEdge.u.y, nextEdge.v.x, nextEdge.v.y, 
                    "stroke=\"red\" stroke-width=\"3\"");
            svg.isave();
            // --------------------------------------------------------------------------------------------------------
        }
        
        return minST;
    }
    
    /**
     * Creates a minimum spanning tree for a graph using the cost function specified in the constructor.
     * @param tEdges the edges of the graph
     * @param root the root node to start the spanning tree
     * @param costFunction a function returning a cost value for a {@link TEdge}
     * @return the spanning tree
     */
    public static Tree<KVector> createSpanningTree(final Set<TEdge> tEdges, final KVector root, 
            final ICostFunction costFunction) {
        return createSpanningTree(tEdges, root, costFunction, null);
    }
}
