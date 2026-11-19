/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
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
package org.eclipse.elk.alg.disco.transform;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Function;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

import org.eclipse.elk.graph.ElkEdge;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.ElkPort;
import org.eclipse.elk.graph.util.ElkGraphUtil;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;

/**
 * Class for splitting an {@link ElkNode} into a list of its components using a simple depth-first search.
 */
public final class ElkGraphComponentsProcessor {

    ///////////////////////////////////////////////////////////////////////////////
    // Variables

    /** List of already visited nodes during dfs. */
    private static List<ElkNode> visited = Lists.newArrayList();
    /** Holds a set of incident nodes to each node of the graph regardless of edge directions. */
    private static Map<ElkNode, Set<ElkNode>> incidenceMap = Maps.newHashMap();
    /** Holds components after splitting. */
    private static List<List<ElkNode>> components = Lists.newArrayList();
    // private List<ElkNode> children;

    ///////////////////////////////////////////////////////////////////////////////
    // Constructor (unused)

    /**
     * This utility class isn't supposed to be instantiated.
     */
    private ElkGraphComponentsProcessor() {
        super();
    }

    ///////////////////////////////////////////////////////////////////////////////
    // public method

    /**
     * Splits the graph into a list of {@link ElkNode ElkNodes}, each representing a connected component.
     * 
     * @param graph
     *            Graph to split into connected components
     * @return Connected components as a list
     */
    public static List<List<ElkNode>> split(final ElkNode graph) {
        visited = Lists.newArrayList();
        incidenceMap = Maps.newHashMap();
        components = Lists.newArrayList();
        List<ElkNode> children = graph.getChildren();
        computeIncidences(children);
        for (ElkNode node : children) {
            if (!visited.contains(node)) {
                List<ElkNode> component = Lists.newArrayList();
                components.add(component);
                dfs(node, component);
            }
        }

        return components;
    }

    ///////////////////////////////////////////////////////////////////////////////
    // private helper methods

    /**
     * Computes the incidences between each node of the graph.
     * 
     * @param nodes
     *            graph to split as a collection of {@link ElkNode ElkNodes} without the surrounding parent node
     */
    private static <C extends Collection<ElkNode>> void computeIncidences(final C nodes) {
        Map<Boolean, List<ElkEdge>> edges;
        Set<ElkNode> adjacentNodes;
        List<ElkEdge> connectedToParentPort;
        Map<ElkPort, Set<ElkNode>> adjacentAndInsideParent = Maps.newHashMap();
        for (ElkNode node : nodes) {
            adjacentNodes = Sets.newHashSet();
            incidenceMap.put(node, adjacentNodes);

            Predicate<? super ElkEdge> sameHierarchyLevel = edge -> ElkGraphUtil.getSourceNode(edge).getParent()
                    .equals(ElkGraphUtil.getTargetNode(edge).getParent());

            edges = StreamSupport.stream(ElkGraphUtil.allIncomingEdges(node).spliterator(), true)
                    .collect(Collectors.partitioningBy(sameHierarchyLevel));
            addNodesToIncidenceSet(adjacentNodes, edges.get(true), edge -> ElkGraphUtil.getSourceNode(edge));

            // Filter because edges might also go one level down instead of up
            connectedToParentPort = edges.get(false).parallelStream().filter(
                    edge -> ElkGraphUtil.getSourceNode(edge).equals(ElkGraphUtil.getTargetNode(edge).getParent()))
                    .collect(Collectors.toList());

            ElkPort port;
            Set<ElkNode> nodesAtPort;

            for (ElkEdge incomingEdge : connectedToParentPort) {
                port = ElkGraphUtil.getSourcePort(incomingEdge);

                if (port != null) {
                    nodesAtPort = adjacentAndInsideParent.get(port);

                    if (nodesAtPort == null) {
                        nodesAtPort = getInnerNeighborsOfPort(port);
                        adjacentAndInsideParent.put(port, nodesAtPort);
                    }

                    adjacentNodes.addAll(nodesAtPort);
                }
            }

            edges = StreamSupport.stream(ElkGraphUtil.allOutgoingEdges(node).spliterator(), true)
                    .collect(Collectors.partitioningBy(sameHierarchyLevel));
            addNodesToIncidenceSet(adjacentNodes, edges.get(true), edge -> ElkGraphUtil.getTargetNode(edge));

            // Filter because edges might also go one level down instead of up
            connectedToParentPort = edges.get(false).parallelStream().filter(
                    edge -> ElkGraphUtil.getTargetNode(edge).equals(ElkGraphUtil.getSourceNode(edge).getParent()))
                    .collect(Collectors.toList());

            for (ElkEdge outgoingEdge : connectedToParentPort) {
                port = ElkGraphUtil.getTargetPort(outgoingEdge);

                if (port != null) {
                    nodesAtPort = adjacentAndInsideParent.get(port);

                    if (nodesAtPort == null) {
                        nodesAtPort = getInnerNeighborsOfPort(port);
                        adjacentAndInsideParent.put(port, nodesAtPort);
                    }

                    adjacentNodes.addAll(nodesAtPort);
                }
            }
        }
    }

    /**
     * Gets all edges connected to a port that are coming from or leading to the inside of its parent node.
     * 
     * @param port
     *            port to be inspected
     * @return inner edges
     */
    private static Set<ElkNode> getInnerNeighborsOfPort(final ElkPort port) {
        ElkNode portParent = port.getParent();
        Predicate<ElkEdge> inwardsPredicate = edge -> portParent.equals(ElkGraphUtil.getSourceNode(edge).getParent())
                || portParent.equals(ElkGraphUtil.getTargetNode(edge).getParent());
        Function<ElkEdge, ElkNode> nodeMapper =
                edge -> portParent.equals(ElkGraphUtil.getSourceNode(edge)) ? ElkGraphUtil.getTargetNode(edge)
                        : ElkGraphUtil.getSourceNode(edge);
        ArrayList<ElkEdge> allEdges = Lists.newArrayList();
        allEdges.addAll(port.getIncomingEdges());
        allEdges.addAll(port.getOutgoingEdges());
        return allEdges.parallelStream().filter(inwardsPredicate).map(nodeMapper).collect(Collectors.toSet());
    }

    /**
     * Extracts nodes from edges based on a given function {@code chooseNode} and adds them to a given set of
     * {@link ElkNode ElkNodes}.
     * 
     * @param incidentNodes
     *            Set to add incident nodes to
     * @param edges
     *            edges to test
     * @param chooseNode
     *            function determining which of the two nodes of an edge to add to the set
     */
    private static void addNodesToIncidenceSet(final Set<ElkNode> incidentNodes, final Collection<ElkEdge> edges,
            final Function<ElkEdge, ElkNode> chooseNode) {
        for (ElkEdge edge : edges) {
            incidentNodes.add(chooseNode.apply(edge));
        }
    }

    /**
     * Naive recursive implementation of depth-first search (dfs).
     * 
     * @param start Node to start with
     * @param component List to add all connected nodes to.
     */
    private static void dfs(final ElkNode start, final List<ElkNode> component) {
        visited.add(start);
        component.add(start);
        Set<ElkNode> adjacentNodes = incidenceMap.get(start);
        if (adjacentNodes != null) {
            for (ElkNode node : adjacentNodes) {
                if (!visited.contains(node)) {
                    dfs(node, component);
                }
            }
        }

    }
}
