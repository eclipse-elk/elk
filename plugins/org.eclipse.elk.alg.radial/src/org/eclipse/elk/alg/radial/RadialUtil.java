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
package org.eclipse.elk.alg.radial;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.graph.ElkEdge;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.ElkPort;
import org.eclipse.elk.graph.util.ElkGraphUtil;

import com.google.common.collect.Iterables;
import com.google.common.math.DoubleMath;

/** A class for smaller, independent calculation units. */
public final class RadialUtil {

    /** Constant for 2*PI. */
    private static final double TWO_PI = 2 * Math.PI;
    /** Constant for fuzzy compare. */
    private static final double EPSILON = 1e-10;

    /** Constructor should not be public visible. */
    private RadialUtil() {
        // Do nothing
    }

    /**
     * The nodes which succedes the given node in the tree. These are defined by the targets of all outgoing edges of
     * the node.
     * 
     * @param node
     *            A node.
     * @return List of neighbors with outgoing edges.
     */
    public static List<ElkNode> getSuccessors(final ElkNode node) {
        List<ElkNode> successors = new ArrayList<>();
        HashSet<ElkNode> children = new HashSet<ElkNode>(node.getChildren());
        for (ElkEdge outgoingEdge : ElkGraphUtil.allOutgoingEdges(node)) {
            if (!(outgoingEdge.getSources().get(0) instanceof ElkPort)) {
                ElkNode target = ElkGraphUtil.connectableShapeToNode(outgoingEdge.getTargets().get(0));
                if (!children.contains(target)) {
                    successors.add(target);
                }
            }
        }
        return successors;
    }

    /**
     * The nodes which succedes the given node in the tree. These are defined by the targets of all outgoing edges of
     * the node.
     * 
     * @param node
     *            A node.
     * @return Set of neighbors with outgoing edges.
     */
    public static Set<ElkNode> getSuccessorSet(final ElkNode node) {
        Set<ElkNode> successors = new HashSet<>();
        HashSet<ElkNode> children = new HashSet<ElkNode>(node.getChildren());
        for (ElkEdge outgoingEdge : ElkGraphUtil.allOutgoingEdges(node)) {
            if (!(outgoingEdge.getSources().get(0) instanceof ElkPort)) {
                ElkNode target = ElkGraphUtil.connectableShapeToNode(outgoingEdge.getTargets().get(0));
                if (!children.contains(target)) {
                    successors.add(target);
                }
            }
        }
        return successors;
    }

    /**
     * Computes the root node of a graph. That means in this case the root of a tree and not the surrounding root. A
     * hierarchical graph is always a tree.
     * 
     * @param graph
     * @return Root node of graph.
     */
    public static ElkNode findRoot(final ElkNode graph) {
        for (ElkNode child : graph.getChildren()) {
            Iterable<ElkEdge> incomingEdges = ElkGraphUtil.allIncomingEdges(child);
            if (!incomingEdges.iterator().hasNext()) {
                return child;
            }
        }
        return null;
    }

    /**
     * Given a node of the tree, calculate the root of the tree.
     * 
     * @param elkNode
     * @return
     */
    public static ElkNode findRootOfNode(final ElkNode elkNode) {
        ElkNode parent = RadialUtil.getTreeParent(elkNode);
        if (parent == null) {
            return elkNode;
        } else {
            return findRootOfNode(parent);
        }
    }

    /**
     * Computes the number of leaves that a node has.
     * 
     * @param node
     * @return number of leaves.
     */
    public static int getNumberOfLeaves(final ElkNode node) {
        int leafs = 0;
        List<ElkNode> successors = getSuccessors(node);
        if (successors.isEmpty()) {
            return 1;
        } else {
            for (final ElkNode child : successors) {
                leafs += getNumberOfLeaves(child);
            }
        }
        return leafs;
    }

    /**
     * Comparator for sorting a node by polar coordinates. Offset is set between 0 and 2*Pi. 0 starts the sorting in the
     * right middle of a node and sorts clockwise.
     * 
     * @param node
     * @param radialOffset
     * @return Comparator for polar sorting.
     */
    public static Comparator<ElkNode> createPolarComparator(final double radialOffset, final double nodeOffsetY) {
        Comparator<ElkNode> comparator = (node1, node2) -> {

            KVector position1 = node1.getProperty(CoreOptions.POSITION);

            double xPos1 = position1.x;
            double yPos1 = position1.y + nodeOffsetY;
            double arc1 = Math.atan2(yPos1, xPos1);
            if (arc1 < 0) {
                arc1 += TWO_PI;
            }
            arc1 += radialOffset;
            if (arc1 > TWO_PI) {
                arc1 -= TWO_PI;
            }

            KVector position2 = node2.getProperty(CoreOptions.POSITION);

            double xPos2 = position2.x;
            double yPos2 = position2.y + nodeOffsetY;
            double arc2 = Math.atan2(yPos2, xPos2);
            if (arc2 < 0) {
                arc2 += TWO_PI;
            }
            arc2 += radialOffset;
            if (arc2 > TWO_PI) {
                arc2 -= TWO_PI;
            }

            return DoubleMath.fuzzyCompare(arc1, arc2, EPSILON);
        };
        return comparator;
    }

    /**
     * Search for the largest diameter of all nodes.
     * 
     * @param graph
     * @return
     */
    public static double findLargestNodeInGraph(final ElkNode graph) {
        double largestChildSize = 0;

        for (ElkNode child : graph.getChildren()) {
            double width = child.getWidth();
            double height = child.getHeight();
            double diameter = Math.sqrt(width * width + height * height);
            largestChildSize = Math.max(diameter, largestChildSize);

            double largestChild = findLargestNodeInGraph(child);
            largestChildSize = Math.max(largestChild, largestChildSize);
        }
        return largestChildSize;
    }

    /**
     * Returns the next level of successor from a list of nodes.
     * 
     * @param nodes
     * @return
     */
    public static List<ElkNode> getNextLevelNodes(final List<ElkNode> nodes) {
        List<ElkNode> successors = new ArrayList<ElkNode>();
        for (ElkNode node : nodes) {
            List<ElkNode> nextLevelNodes = getSuccessors(node);
            successors.addAll(nextLevelNodes);
        }
        return successors;
    }
    
    /**
     * Get set of the nodes in the next level.
     * @param nodes The nodes in the current level
     * @return The next nodes
     */
    public static Set<ElkNode> getNextLevelNodeSet(final List<ElkNode> nodes) {
        Set<ElkNode> successors = new HashSet<ElkNode>();
        for (ElkNode node : nodes) {
            Set<ElkNode> nextLevelSet = getSuccessorSet(node);
            successors.addAll(nextLevelSet);
        }
        return successors;
    }

    /**
     * Returns the next level of successor from a list of nodes.
     * 
     * @param nodes
     * @return
     */
    public static List<ElkNode> getLastLevelNodes(final List<ElkNode> nodes) {
        List<ElkNode> parents = new ArrayList<ElkNode>();
        for (ElkNode node : nodes) {
            ElkNode parent = getTreeParent(node);
            parents.add(parent);
        }
        return parents;

    }

    /**
     * Shift a node such that it is centered on the radius.
     * 
     * @param node
     *            A node to center.
     * @param xPos
     *            The nodes new x Position on the radius.
     * @param yPos
     *            The nodes new y Position on the radius.
     */
    public static void centerNodesOnRadi(final ElkNode node, final double xPos, final double yPos) {
        double xPosition = xPos - node.getWidth() / 2;
        double yPosition = yPos - node.getHeight() / 2;

        node.setX(xPosition);
        node.setY(yPosition);
    }

    /**
     * Shift nodes such that the edge closest to the parent node, is closest now.
     * 
     * @param node
     *            A node to shift.
     * @param xPos
     *            The nodes new x Position.
     * @param yPos
     *            The nodes new y Position.
     */
    public static void shiftClosestEdgeToRadi(final ElkNode node, final double xPos, final double yPos) {
        // center root
        if (DoubleMath.fuzzyEquals(xPos, 0, EPSILON) && DoubleMath.fuzzyEquals(yPos, 0, EPSILON)) {
            node.setX(xPos - node.getWidth() / 2);
            node.setY(yPos - node.getHeight() / 2);
        } else {

            if (xPos < 0) {
                if (yPos < 0) {
                    // lower left rectangle
                    node.setX(xPos - node.getWidth());
                    node.setY(yPos);
                } else {
                    // upper left rectangle
                    node.setX(xPos - node.getWidth());
                    node.setY(yPos + node.getHeight());
                }
            } else {

                if (yPos < 0) {
                    // lower right rectangle
                    node.setX(xPos);
                    node.setY(yPos);
                } else {
                    // upper right rectangle
                    node.setX(xPos);
                    node.setY(yPos + node.getHeight());
                }
            }
        }
    }

    /**
     * Computes the nodes which is the predecessor in the tree.
     * 
     * @param node
     * @return
     */
    public static ElkNode getTreeParent(final ElkNode node) {
        Iterable<ElkEdge> iterator = ElkGraphUtil.allIncomingEdges(node);
        if (!Iterables.isEmpty(iterator)) {
            ElkEdge edgeFromParent = Iterables.get(iterator, 0);
            return ElkGraphUtil.connectableShapeToNode(edgeFromParent.getSources().get(0));
        } else {
            return null;
        }
    }

}
