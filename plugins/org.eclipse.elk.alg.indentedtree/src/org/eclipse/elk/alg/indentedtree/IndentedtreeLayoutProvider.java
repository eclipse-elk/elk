/*******************************************************************************
 * Copyright (c) 2024 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.alg.indentedtree;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

import org.eclipse.elk.core.AbstractLayoutProvider;
import org.eclipse.elk.core.UnsupportedGraphException;
import org.eclipse.elk.core.math.ElkPadding;
import org.eclipse.elk.core.util.IElkProgressMonitor;
import org.eclipse.elk.graph.ElkEdge;
import org.eclipse.elk.graph.ElkEdgeSection;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.util.ElkGraphUtil;
import org.eclipse.elk.alg.indentedtree.options.IndentedtreeOptions;
import org.eclipse.elk.alg.indentedtree.options.InternalProperties;

/**
 * A layout algorithm for creating a filesystem-like layout. This assumes that the graph is a tree
 * and that no ports exist
 * 
 */
public class IndentedtreeLayoutProvider extends AbstractLayoutProvider {

    @Override
    public void layout(ElkNode layoutGraph, IElkProgressMonitor progressMonitor) {
        // Start progress monitor
        progressMonitor.begin("org.eclipse.elk.alg.indentedtree", 2);
        progressMonitor.log("Algorithm began");

        // Retrieve padding
        ElkPadding padding = layoutGraph.getProperty(IndentedtreeOptions.PADDING);

        // Get the list of nodes to lay out
        List<ElkNode> nodes = new ArrayList<>(layoutGraph.getChildren());

        // Create a sub monitor for node placement
        IElkProgressMonitor nodePlacingMonitor = progressMonitor.subTask(1);
        nodePlacingMonitor.begin("Node Spacing", nodes.size());

        // Positions
        double currY = padding.top;
        double maxX = 0.0;

        // Make an output to the debug log
        nodePlacingMonitor.log("left padding: " + padding.getLeft());
        nodePlacingMonitor.log("top padding: " + padding.getTop());

        nodePlacingMonitor.logGraph(layoutGraph, "No node placed yet");

        // find sources
        List<ElkNode> sources = new LinkedList<>();
        for (ElkNode node : nodes) {
            if (node.getIncomingEdges().isEmpty()) {
                sources.add(node);
            }
        }

        // declare outside of loop for later use when setting graph bounds
        double currentVerticalSpacing = 0.0;

        // iterate through every tree in the forest
        for (ElkNode source : sources) {
            List<ElkNode> orderedNodes = new LinkedList<>();
            // find ordered list of nodes
            dfs(source, orderedNodes, 0.0);
            // iterate through a single tree in order, as determined by depth-first-search
            for (ElkNode n : orderedNodes) {
                double currentDepth = padding.getLeft() + n.getProperty(InternalProperties.DEPTH);
                n.setX(currentDepth);
                n.setY(currY); // currY increases monotonically

                // update maxX in order to know the real bounds of the graph
                maxX = Math.max(maxX, currentDepth + n.getWidth());

                currentVerticalSpacing = n.getProperty(IndentedtreeOptions.VERTICAL_NODE_SPACING);
                // increase y counter by the node's height plus its vertical node spacing
                currY += n.getHeight() + currentVerticalSpacing;

                nodePlacingMonitor.log("currX: " + n.getProperty(InternalProperties.DEPTH));
                nodePlacingMonitor.logGraph(layoutGraph, n.getIdentifier() + " placed");
            }

        }

        // Close the sub monitor
        nodePlacingMonitor.done();
        progressMonitor.log("Node Placing done!");

        // Create sub monitor for edge routing
        IElkProgressMonitor edgeRoutingMonitor = progressMonitor.subTask(1);
        edgeRoutingMonitor.begin("Edge Routing", layoutGraph.getContainedEdges().size());
        edgeRoutingMonitor.logGraph(layoutGraph, "No edge routed yet");

        // Route the edges
        if (!layoutGraph.getContainedEdges().isEmpty()) {

            // iterate through all edges in the graph
            for (ElkEdge edge : layoutGraph.getContainedEdges()) {
                // obtain source and target
                ElkNode source = ElkGraphUtil.connectableShapeToNode(edge.getSources().get(0));
                ElkNode target = ElkGraphUtil.connectableShapeToNode(edge.getTargets().get(0));

                // create one section for the new edge
                ElkEdgeSection section = ElkGraphUtil.firstEdgeSection(edge, true, true);

                

                // if horizontal edge indentation is set, use that, otherwise, use half of the
                // parent's width
                double horizontalEdgeIndentation;
                if (source.hasProperty(IndentedtreeOptions.HORIZONTAL_EDGE_INDENTATION)) {
                    horizontalEdgeIndentation =
                            source.getProperty(IndentedtreeOptions.HORIZONTAL_EDGE_INDENTATION);
                } else {
                    horizontalEdgeIndentation = source.getWidth() / 2;
                }
                
                // throw error if edge would have to be placed far away from the parent node
                if (source.getWidth() < horizontalEdgeIndentation) {
                    throw new UnsupportedGraphException(
                            "Horizontal Edge Indentation is larger than its parent node's width");
                }

                // set starting location for new edge
                section.setStartLocation(source.getX() + horizontalEdgeIndentation,
                        source.getY() + source.getHeight());

                // end location: middle of left border of target
                section.setEndLocation(target.getX(), target.getY() + target.getHeight() / 2);

                // create a bend point at "longitude" of start point and "latitude" of end point
                ElkGraphUtil.createBendPoint(section, section.getStartX(), section.getEndY());

                edgeRoutingMonitor.log("bendPoint created at y=" + currY);
                edgeRoutingMonitor.logGraph(layoutGraph,
                        source.getIdentifier() + " -> " + target.getIdentifier());
            }

        }

        // Close the sub monitor
        edgeRoutingMonitor.done();

        progressMonitor.log("Edge Routing done!");

        // Set the size of the final diagram
        layoutGraph.setWidth(maxX + padding.right); // TODO should left/top padding be respected
                                                    // here?
        layoutGraph.setHeight(currY + padding.bottom - currentVerticalSpacing);

        // End the progress monitor
        progressMonitor.log("Algorithm executed");
        progressMonitor.logGraph(layoutGraph, "Final graph");
        progressMonitor.done();
    }

    /**
     * Depth-first-search through a tree
     * 
     * @param origin
     *                   the root node of the current subtree
     * @param list
     *                   the list of all nodes already visited
     * @param depth
     *                   the accumulated horizontal indentation of the origin
     * @return a list of all nodes in the (sub-)tree beginning with origin
     */
    private List<ElkNode> dfs(ElkNode origin, List<ElkNode> list, Double depth) {
        // add received depth to relativeIndentation and save to property
        depth = depth + origin.getProperty(IndentedtreeOptions.RELATIVE_INDENTATION);
        origin.setProperty(InternalProperties.DEPTH, depth);

        list.add(origin);
        // iterate through origin's "children"
        for (ElkEdge out : origin.getOutgoingEdges()) {
            // recursion
            dfs(ElkGraphUtil.connectableShapeToNode(out.getTargets().get(0)), list,
                    depth + (origin.hasProperty(IndentedtreeOptions.HORIZONTAL_EDGE_INDENTATION) ? origin.getProperty(IndentedtreeOptions.HORIZONTAL_EDGE_INDENTATION) : origin.getWidth()/2));
        }
        return list;
    }
}
