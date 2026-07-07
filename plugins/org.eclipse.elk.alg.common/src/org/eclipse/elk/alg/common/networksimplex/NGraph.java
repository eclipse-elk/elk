/*******************************************************************************
 * Copyright (c) 2016 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.common.networksimplex;

import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.options.Direction;
import org.eclipse.elk.graph.ElkEdge;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.util.ElkGraphUtil;
import org.eclipse.emf.common.util.URI;
import org.eclipse.emf.ecore.resource.Resource;
import org.eclipse.emf.ecore.resource.ResourceSet;
import org.eclipse.emf.ecore.resource.impl.ResourceSetImpl;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

/**
 * A graph structure used by the {@link NetworkSimplex} algorithm.
 * The class holds a list of nodes and provides some convenient methods. 
 */
public class NGraph {

    // SUPPRESS CHECKSTYLE NEXT 2 VisibilityModifier
    /** The nodes of the network simplex graph. */
    public List<NNode> nodes = Lists.newArrayList();
    
    /**
     * Converts this {@link NGraph} to a KGraph and writes it to the specified filed.
     * 
     * @param filePath a path to a file on the filesystem
     */
    public void writeDebugGraph(final String filePath) {
        
        ElkNode elkGraph = ElkGraphUtil.createGraph();
        elkGraph.setProperty(CoreOptions.DIRECTION, Direction.DOWN);
        Map<NNode, ElkNode> nodeMap = Maps.newHashMap();
        
        for (NNode nNode : nodes) {
            ElkNode elkNode = ElkGraphUtil.createNode(elkGraph);
            nodeMap.put(nNode, elkNode);
            
            ElkGraphUtil.createLabel(nNode.type + " " + nNode.layer, elkNode);
        }

        for (NNode nNode : nodes) {
            for (NEdge nEdge : nNode.getOutgoingEdges()) {
                ElkEdge elkEdge = ElkGraphUtil.createSimpleEdge(nodeMap.get(nEdge.source), nodeMap.get(nEdge.target));
                ElkGraphUtil.createLabel(nEdge.weight + " " + nEdge.delta, elkEdge);
            }
        }
        
        ResourceSet rs = new ResourceSetImpl();
        Resource r = rs.createResource(URI.createFileURI(filePath));
        r.getContents().add(elkGraph);
        try {
            r.save(Collections.emptyMap());
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
    
    /**
     * Checks if this {@link NGraph} is connected. If not, it identifies one representative per connected component
     * and connects each of them to a new artificial root node. 
     */
    public final NNode makeConnected() {
        int id = 0;
        for (NNode n : nodes) {
            n.internalId = id++;
        }
        List<NNode> ccRep = findConCompRepresentatives();
        NNode root = null;
        if (ccRep.size() > 1) {
            root = createArtificialRootAndConnect(ccRep);
        }
        return root;
    }
    
    /**
     * Creates and returns a new {@link NNode} and connects the passed {@code nodesToConnect}
     * to the new node with zero-weight, zero-delta edges.
     */
    private NNode createArtificialRootAndConnect(final List<NNode> nodesToConnect) {
        NNode root = NNode.of().create(this);
        for (NNode src : nodesToConnect) {
            NEdge.of()
               .delta(0)
               .weight(0)
               .source(root)
               .target(src)
               .create();
        }
        return root;
    }
    
    /**
     * Identifies and returns one representative {@link NNode} per connected component. 
     */
    private List<NNode> findConCompRepresentatives() {
        List<NNode> ccRep = Lists.newArrayList();
        boolean[] mark = new boolean[nodes.size()];
        Arrays.fill(mark, false);
        for (NNode node : nodes) {
            if (!mark[node.internalId]) {
                ccRep.add(node);
                dfs(node, mark);
            }
        }
        return ccRep;
    }
    
    private void dfs(final NNode node, final boolean[] mark) {
        if (mark[node.internalId]) {
            return;
        }
        mark[node.internalId] = true;
        for (NEdge edge : node.getConnectedEdges()) {
            NNode other = edge.getOther(node);
            dfs(other, mark);
        }
    }
    
    
    /**
     * Creates a topological ordering and checks for back edges.
     * 
     * @return true if the graph is acyclic, false if it is cyclic.
     */
    public boolean isAcyclic() {

        int id = 0;
        for (NNode n : nodes) {
            n.internalId = id++;
        }
        
        // initialize the number of incident edges for each node
        int[] incident = new int[nodes.size()];
        int[] layer = new int[nodes.size()];
        for (NNode node : nodes) {
            incident[node.internalId] += node.getIncomingEdges().size();
        }

        LinkedList<NNode> roots = Lists.newLinkedList();
        for (NNode node : nodes) {
            if (node.getIncomingEdges().isEmpty()) {
                roots.add(node);
            }
        }
        if (roots.isEmpty() && !nodes.isEmpty()) {
            return false;
        }
        while (!roots.isEmpty()) {
            NNode node = roots.poll();
            
            for (NEdge edge : node.getOutgoingEdges()) {
                NNode target = edge.getTarget();
                layer[target.internalId] = Math.max(layer[target.internalId], layer[node.internalId] + 1);
                incident[target.internalId]--;
                if (incident[target.internalId] == 0) {
                    roots.add(target);
                }
            }
        }
        
        // check for backward edges
        for (NNode node : nodes) {
            for (NEdge edge : node.getOutgoingEdges()) {
                if (layer[edge.target.internalId] <= layer[edge.source.internalId]) {
                    return false;
                }
            }
        }
        
        return true;
    }
}
