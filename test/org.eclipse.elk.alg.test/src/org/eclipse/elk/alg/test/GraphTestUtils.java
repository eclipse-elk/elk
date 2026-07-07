/*******************************************************************************
 * Copyright (c) 2016, 2020 Kiel University and others.
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
package org.eclipse.elk.alg.test;

import java.awt.geom.Rectangle2D;
import java.util.ArrayList;
import java.util.List;
import java.util.Objects;

import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.util.ElkUtil;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.ElkPort;
import org.eclipse.elk.graph.ElkShape;
import org.eclipse.elk.graph.util.ElkGraphUtil;

/**
 * Utility class for creating simple test graphs.
 */
public final class GraphTestUtils {

    private GraphTestUtils() { }
    
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Graph Traversal
    
    /**
     * Returns all nodes in the given graph, with or without the graph itself.
     */
    public static List<ElkNode> allNodes(final ElkNode graph, final boolean includeGraph) {
        List<ElkNode> nodes = new ArrayList<>();
        
        if (includeGraph) {
            nodes.add(graph);
        } else {
            nodes.addAll(graph.getChildren());
        }
        
        // We walk through the list, adding all children we find in the process. It's important here that the list's
        // size can grow during each iteration and thus cannot be cached
        for (int i = 0; i < nodes.size(); i++) {
            nodes.addAll(nodes.get(i).getChildren());
        }
        
        return nodes;
    }
    
    /**
     * Returns whether any of the shapes in the list overlap.
     */
    public static boolean haveOverlaps(final List<? extends ElkShape> shapes) {
        // Obtain absolute rectangles for us to do the comparisons on
        List<Rectangle2D> shapeRects = new ArrayList<>(shapes.size());
        for (ElkShape shape : shapes) {
            KVector pos = ElkUtil.absolutePosition(shape);
            shapeRects.add(new Rectangle2D.Double(pos.x, pos.y, shape.getWidth(), shape.getHeight()));
        }

        for (int first = 0; first < shapeRects.size(); first++) {
            Rectangle2D firstRect = shapeRects.get(first);
            
            for (int second = first + 1; second < shapeRects.size(); second++) {
                Rectangle2D secondRect = shapeRects.get(second);
                
                if (firstRect.intersects(secondRect)) {
                    return true;
                }
            }
        }
        
        return false;
    }
    
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Graph Creation

    /**
     * Creates a KGraph, represented by a parent node that contains the actual nodes and edges of the graph.
     * 
     * @return a parent node that contains graph elements
     */
    public static ElkNode createSimpleGraph() {
        // create parent node
        ElkNode parentNode = ElkGraphUtil.createGraph();

        // create child nodes with labels
        ElkNode childNode1 = ElkGraphUtil.createNode(parentNode);
        ElkGraphUtil.createLabel("node1", childNode1);
        ElkNode childNode2 = ElkGraphUtil.createNode(parentNode);
        ElkGraphUtil.createLabel("node2", childNode2);

        // create ports (optional)
        ElkPort port1 = ElkGraphUtil.createPort(childNode1);
        ElkPort port2 = ElkGraphUtil.createPort(childNode2);

        // create a simple edge
        ElkGraphUtil.createSimpleEdge(port1, port2);

        return parentNode;
    }

    /**
     * Returns a basic graph with three connected nodes.
     */
    public static ElkNode createThreeNodeGraph() {
        ElkNode layoutGraph = ElkGraphUtil.createGraph();
        ElkNode node1 = ElkGraphUtil.createNode(layoutGraph);
        ElkNode node2 = ElkGraphUtil.createNode(layoutGraph);
        ElkNode node3 = ElkGraphUtil.createNode(layoutGraph);

        ElkGraphUtil.createSimpleEdge(node1, node3);
        ElkGraphUtil.createSimpleEdge(node2, node3);

        // layout options
        for (ElkNode node : layoutGraph.getChildren()) {
            // Suppress the Checkstyle warnings, because the meaning is clear
            // SUPPRESS CHECKSTYLE NEXT 2 MagicNumber
            node.setWidth(50.0);
            node.setHeight(50.0);
        }

        return layoutGraph;
    }

    /**
     * @return a hierarchical graph composed of two {@link #createSimpleGraph()}s.
     */
    public static ElkNode createHierarchicalGraph() {
        // create parent node
        ElkNode parentNode = ElkGraphUtil.createGraph();

        ElkNode child1 = createSimpleGraph();
        child1.setParent(parentNode);
        ElkGraphUtil.createLabel("child1", child1);

        ElkNode child2 = createSimpleGraph();
        child2.setParent(parentNode);
        ElkGraphUtil.createLabel("child2", child2);

        return parentNode;
    }

    /**
     * Returns the {@link ElkNode} corresponding to the {@code childPath}.
     *
     * @param parent
     *            a parent in the hierarchy of the wanted child.
     * @param childPath
     *            a list of identifier ({@link ElkShape#getIdentifier()}) to go from {@code parent} to the wanted child.
     * @return the {@link ElkNode} corresponding to the {@code childPath}
     * @throws NullPointerException
     *             if {@code parent} is {@code null}.
     * @throws IllegalArgumentException
     *             if the {@code childPath} does not match to an existing path
     */
    public static ElkNode getChild(final ElkNode parent, final String... childPath) {
        Objects.requireNonNull(parent, "parent cannot be null");
        ElkNode result = null;
        ElkNode currentParent = parent;
        for (String identifier : childPath) {
            result = null;
            for (ElkNode child : currentParent.getChildren()) {
                if (identifier.equals(child.getIdentifier())) {
                    result = child;
                    break;
                }
            }
            if (result == null) {
                throw new IllegalArgumentException("The node \"" + currentParent.getIdentifier()
                        + "\" has no child with identifier \"" + identifier + "\".");
            } else {
                currentParent = result;
            }
        }
        return result;
    }

    /**
     * Returns the {@link ElkPort} of this {@code node} with the identifier {@code portId}.
     *
     * @param node
     *            the node containing the wanted port
     * @param portId
     *            the identifier of the wanted port
     * @return the {@link ElkPort} with the {@code portId}
     * @throws IllegalArgumentException
     *             if the {@code node} does not have a port with this identifier
     */
    public static ElkPort getPort(final ElkNode node, final String portId) {
        return node.getPorts().stream().filter(p -> portId.equals(p.getIdentifier())).findFirst()
                .orElseThrow(() -> new IllegalArgumentException(
                        "The node \"" + node.getIdentifier() + "\" has no port with identifier \"" + portId + "\"."));
    }

}
