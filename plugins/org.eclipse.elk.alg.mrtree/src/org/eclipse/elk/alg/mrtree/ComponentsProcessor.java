/*******************************************************************************
 * Copyright (c) 2013 - 2022 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.alg.mrtree;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.stream.Collectors;

import org.eclipse.elk.alg.mrtree.graph.TEdge;
import org.eclipse.elk.alg.mrtree.graph.TGraph;
import org.eclipse.elk.alg.mrtree.graph.TNode;
import org.eclipse.elk.alg.mrtree.options.InternalProperties;
import org.eclipse.elk.alg.mrtree.options.MrTreeOptions;
import org.eclipse.elk.core.math.ElkPadding;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.graph.properties.IProperty;

import com.google.common.collect.Lists;

/**
 * A processor that is able to split an input graph into connected components and to pack those
 * components after layout.
 * 
 * @author sor
 * @author sgu
 * @author msp
 * @author cds
 * @author sdo
 */
public class ComponentsProcessor {

    /** the incidence lists of the graph. */
    private List<TEdge>[] incidence;
    /** the visited array of the graph. */
    private boolean[] visited;

    /**
     * Initialize adjacency lists.
     * 
     * @param graph
     *            a force graph
     */
    @SuppressWarnings("unchecked")
    private void initialize(final TGraph graph) {
        int n = graph.getNodes().size();
        incidence = new List[n];
        visited = new boolean[n];

        // create incidence lists
        for (TNode node : graph.getNodes()) {
            incidence[node.id] = new LinkedList<TEdge>();
        }
        // add edges to incidence lists
        for (TEdge edge : graph.getEdges()) {
            incidence[edge.getSource().id].add(edge);
            incidence[edge.getTarget().id].add(edge);
        }
    }

    /**
     * Split the given graph into its connected components.
     * 
     * @param graph
     *            an input graph
     * @return a list of components that can be processed one by one
     */
    public List<TGraph> split(final TGraph graph) {
        Boolean separate = graph.getProperty(MrTreeOptions.SEPARATE_CONNECTED_COMPONENTS);
        if (separate == null || separate.booleanValue()) {
            initialize(graph);

            // perform DFS starting on each node, collecting connected components
            List<TGraph> components = new ArrayList<TGraph>();
            for (TNode node : graph.getNodes()) {
                TGraph comp = dfs(node, null);
                if (comp != null) {
                    comp.copyProperties(graph);
                    components.add(comp);
                }
            }
            incidence = null;
            visited = null;

            // redistribute identifier numbers to each component
            if (components.size() > 1) {
                for (TGraph comp : components) {
                    int id = 0;
                    for (TNode node : comp.getNodes()) {
                        node.id = id++;
                    }
                }
            }
            return components;
        }
        return Lists.newArrayList(graph);
    }

    /**
     * Perform a DFS starting on the given node and collect all nodes that are found in the
     * corresponding connected component.
     * 
     * @param node
     *            the starting node
     * @param graph
     *            a graph representing a connected component, or {@code null}
     * @return the connected component, or {@code null} if the node was already visited
     */
    private TGraph dfs(final TNode node, final TGraph graph) {
        // check if node was visited earlier
        if (!visited[node.id]) {
            // if not it is visited now
            visited[node.id] = true;
            // create new component if necessary
            TGraph component = graph;
            if (component == null) {
                component = new TGraph();
            }
            // add the node to the component
            component.getNodes().add(node);
            // follow all edges from this node
            for (TEdge edge : incidence[node.id]) {
                if (edge.getSource() != node) {
                    dfs(edge.getSource(), component);
                }
                if (edge.getTarget() != node) {
                    dfs(edge.getTarget(), component);
                }
                // add the edges an labels to the component
                component.getEdges().add(edge);
            }
            return component;
        }
        // if the node was visited before skip it
        return null;
    }

    /**
     * Pack the given components into a single graph.
     * 
     * @param components
     *            a list of components
     * @return a single graph that contains all components
     */
    @SuppressWarnings("unchecked")
    public TGraph pack(final List<TGraph> components) {
        if (components.size() == 1) {
            applyPaddingAndNormalizePositions(components.get(0));
            return components.get(0);
        } else if (components.size() <= 0) {
            return new TGraph();
        }

        // assign priorities and sizes
        for (TGraph graph : components) {
            int priority = 0;
            double minx = Integer.MAX_VALUE, miny = Integer.MAX_VALUE;
            double maxx = Integer.MIN_VALUE, maxy = Integer.MIN_VALUE;
            for (TNode node : graph.getNodes()) {
                priority += node.getProperty(MrTreeOptions.PRIORITY);
                minx = Math.min(minx, node.getPosition().x);
                miny = Math.min(miny, node.getPosition().y);
                maxx = Math.max(maxx, node.getPosition().x + node.getSize().x);
                maxy = Math.max(maxy, node.getPosition().y + node.getSize().y);
            }
            graph.setProperty(MrTreeOptions.PRIORITY, priority);
            graph.setProperty(InternalProperties.BB_UPLEFT, new KVector(minx, miny));
            graph.setProperty(InternalProperties.BB_LOWRIGHT, new KVector(maxx, maxy));
        }

        // sort the components by their priority and size
        Collections.sort(components, new Comparator<TGraph>() {
            public int compare(final TGraph graph1, final TGraph graph2) {
                int prio = graph2.getProperty(MrTreeOptions.PRIORITY)
                        - graph1.getProperty(MrTreeOptions.PRIORITY);
                if (prio == 0) {
                    KVector size1 = graph1.getProperty(InternalProperties.BB_LOWRIGHT).clone().sub(
                            graph1.getProperty(InternalProperties.BB_UPLEFT));
                    KVector size2 = graph2.getProperty(InternalProperties.BB_LOWRIGHT).clone().sub(
                            graph2.getProperty(InternalProperties.BB_UPLEFT));
                    return Double.compare(size1.x * size1.y, size2.x * size2.y);
                }
                return prio;
            }
        });

        TGraph result = new TGraph();
        result.copyProperties(components.get(0));

        // determine the maximal row width by the maximal box width and the total area
        double maxRowWidth = 0.0f;
        double totalArea = 0.0f;
        for (TGraph graph : components) {
            KVector size = graph.getProperty(InternalProperties.BB_LOWRIGHT).clone().sub(
                    graph.getProperty(InternalProperties.BB_UPLEFT));
            maxRowWidth = Math.max(maxRowWidth, size.x);
            totalArea += size.x * size.y;
        }
        maxRowWidth = Math.max(maxRowWidth,
                (float) Math.sqrt(totalArea) * result.getProperty(MrTreeOptions.ASPECT_RATIO));
        double spacing = result.getProperty(MrTreeOptions.SPACING_NODE_NODE).doubleValue();

        // place nodes iteratively into rows
        double xpos = 0, ypos = 0, highestBox = 0, broadestRow = spacing;
        for (TGraph graph : components) {
            KVector size = graph.getProperty(InternalProperties.BB_LOWRIGHT).clone().sub(
                    graph.getProperty(InternalProperties.BB_UPLEFT));
            if (xpos + size.x > maxRowWidth) {
                // place the graph into the next row
                xpos = 0;
                ypos += highestBox + spacing;
                highestBox = 0;
            }
            moveGraph(result, graph, xpos, ypos);
            broadestRow = Math.max(broadestRow, xpos + size.x);
            highestBox = Math.max(highestBox, size.y);
            xpos += size.x + spacing;
        }

        Map<IProperty<?>, Object> propMerge = new HashMap<IProperty<?>, Object>();

        Map<IProperty<?>, Object> debug = new HashMap<IProperty<?>, Object>();

        for (TGraph tGraph : components) {
            boolean debugMode = tGraph.getProperty(CoreOptions.DEBUG_MODE);
            Map<IProperty<?>, Object> propComp = tGraph.getAllProperties();
            for (Entry<IProperty<?>, Object> entry : propComp.entrySet()) {
                if (propMerge.containsKey(entry.getKey())) {
                    if (entry.getKey().getDefault() != entry.getValue()) {
                        if (debugMode && debug.containsKey(entry.getKey())) {
                            System.err.println("Found different values for property "
                                    + entry.getKey().getId() + " in components.");
                        } else {
                            propMerge.put(entry.getKey(), entry.getValue());
                            result.setProperty((IProperty<Object>) entry.getKey(), entry.getValue());
                            if (debugMode) {
                                debug.put(entry.getKey(), entry.getValue());
                            }
                        }
                    }
                } else {
                    propMerge.put(entry.getKey(), entry.getValue());
                    result.setProperty((IProperty<Object>) entry.getKey(), entry.getValue());
                }
            }
        }
        
        // Move the resulting graph to 0,0 and apply padding
        applyPaddingAndNormalizePositions(result);


        return result;
    }
    
    
    /**
     * Move the nodes so that minX and minY node coords are zero and apply padding.
     * 
     * @param g the graph to work on
     */
    private void applyPaddingAndNormalizePositions(final TGraph g) {
        ElkPadding padding = g.getProperty(MrTreeOptions.PADDING);
        g.setProperty(InternalProperties.BB_UPLEFT, new KVector(0, 0));
        moveGraph(new TGraph(), g, padding.getLeft() - g.getProperty(InternalProperties.GRAPH_XMIN), 
                padding.getTop() - g.getProperty(InternalProperties.GRAPH_YMIN));
    }


    /**
     * Move the source graph into the destination graph using a specified offset.
     * 
     * @param destGraph
     *            the destination graph
     * @param sourceGraph
     *            the source graph
     * @param offsetx
     *            x coordinate offset
     * @param offsety
     *            y coordinate offset
     */
    private void moveGraph(final TGraph destGraph, final TGraph sourceGraph, final double offsetx,
            final double offsety) {
        KVector graphOffset = new KVector(offsetx, offsety);
        graphOffset.sub(sourceGraph.getProperty(InternalProperties.BB_UPLEFT));

        for (TNode node : sourceGraph.getNodes()) {
            node.getPosition().add(graphOffset);
            destGraph.getNodes().add(node);
        }

        for (TEdge edge : sourceGraph.getEdges().stream().distinct().collect(Collectors.toList())) {
            for (KVector bendpoint : edge.getBendPoints()) {
                bendpoint.add(graphOffset);
            }
            destGraph.getEdges().add(edge);
        }
    }

}
