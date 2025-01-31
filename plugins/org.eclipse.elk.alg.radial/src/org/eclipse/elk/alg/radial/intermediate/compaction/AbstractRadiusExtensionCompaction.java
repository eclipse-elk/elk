/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.alg.radial.intermediate.compaction;

import java.util.List;

import org.eclipse.elk.graph.ElkNode;

/**
 * The class provides basic logic for extending or compacting radii, like overlap calculation.
 */
public class AbstractRadiusExtensionCompaction {
    /** The step size with which the contraction takes place. Default is one. */
    private int compactionStep = 1;

    /** The spacing between nodes which has to be considered while contracting. */
    private double spacing;

    /** The root node of the graph. */
    private ElkNode root;

    /**
     * Contracts/Extends a list of nodes from the same radius by moving them along their incoming edge.
     * 
     * @param layerNodes
     *            List of nodes of one radius that shall be moved.
     * @param isContracting
     *            Determines if the layer shall be contracted or extended.
     */
    public void contractLayer(final List<ElkNode> layerNodes, final boolean isContracting) {
        for (ElkNode node : layerNodes) {

            // node center
            double xPos = node.getX() + node.getWidth() / 2;
            double yPos = node.getY() + node.getHeight() / 2;

            ElkNode treeParent = root;
            double parentX = treeParent.getX() + treeParent.getWidth() / 2;
            double parentY = treeParent.getY() + treeParent.getHeight() / 2;

            // vector of edge
            double x = xPos - parentX;
            double y = yPos - parentY;
            // vector length
            double length = Math.sqrt(x * x + y * y);

            // multiply with normalized vector
            x *= compactionStep / length;
            y *= compactionStep / length;

            if (isContracting) {
                xPos -= x;
                yPos -= y;
            } else {
                xPos += x;
                yPos += y;
            }

            node.setX(xPos - node.getWidth() / 2);
            node.setY(yPos - node.getHeight() / 2);
        }
    }
    
//    /**
//     * Moves children of a node by the amount of the given vector in the direction of their connection to the root.
//     * Since the parent/ancestor node might have moved since nodes in that layer overlapped all children have to be 
//     * moved by the same amount.
//     * @param node The node
//     * @param oldX The nodes old x position
//     * @param oldY The nodes old y position
//     */
//    public void moveChildrenAfterContractLayer(final ElkNode node, final double oldX, final double oldY) {
//        // Calculate direction of movement for
//        for (ElkEdge connection : node.getOutgoingEdges()) {
//            for (ElkConnectableShape target : connection.getTargets()) {
//                ElkNode child = ElkGraphUtil.connectableShapeToNode(target);
//                moveChild(child, movedDistance);
//            }
//        }     
//    }
    
    /**
     * Move the node by by the given distance in the direction of the root node to this node.
     * 
     * @param node The node
     * @param distance The distance to move
     */
    public void moveNode(final ElkNode node, final double distance) {
        // Move child distance into the direction of the connection to the root node
        // We have to measure from center to center since the different nodes might have different sizes.
        ElkNode root = this.root;
        double rootX = root.getX() + root.getWidth() / 2.0;
        double rootY = root.getY() + root.getHeight() / 2.0;
        double nodeX = node.getX() + node.getWidth() / 2.0;
        double nodeY = node.getY() + node.getHeight() / 2.0;
        double differenceX = nodeX - rootX;
        double differenceY = nodeY - rootY;
        // Calculate unit vector
        double length = Math.sqrt(differenceX * differenceX + differenceY * differenceY);
        double unitX = differenceX / length;
        double unitY = differenceY / length;
        // Move node by distance in direction of uni vector.
        node.setX(node.getX() + unitX * distance);
        node.setY(node.getY() + unitY * distance);
        
    }

    /**
     * Calculates if two nodes overlap with each other.
     * 
     * @param node1
     *            The first node.
     * @param node2
     *            A second node.
     * @return If a overlap exists.
     */
    public boolean overlap(final ElkNode node1, final ElkNode node2) {
        double x1 = node1.getX() - spacing / 2;
        double x2 = node2.getX() - spacing / 2;
        double y1 = node1.getY() - spacing / 2;
        double y2 = node2.getY() - spacing / 2;

        double width1 = node1.getWidth() + spacing;
        double width2 = node2.getWidth() + spacing;
        double height1 = node1.getHeight() + spacing;
        double height2 = node2.getHeight() + spacing;

        if ((x1 < x2 + width2 && x2 < x1) && (y1 < y2 + height2 && y2 < y1)) {
            // left upper and right lower corner overlap
            return true;
        } else if ((x2 < x1 + width1 && x1 < x2) && (y2 < y1 + height1 && y1 < y2)) {
            // right lower and left upper corner overlap
            return true;
        } else if ((x1 < x2 + width2 && x2 < x1) && (y1 < y2 && y2 < y1 + height1)) {
            // left lower and right upper corner overlap
            return true;
        } else if ((x2 < x1 + width1 && x1 < x2) && (y1 < y2 + height2 && y2 < y1)) {
            // right upper and left lower corner overlap
            return true;
        }
        return false;
    }

    /**
     * Calculate if the nodes of one radius are overlapping each other.
     * 
     * @param nodes
     *            List of nodes from one radius.
     * @return If the nodes of the radius are overlapping.
     */
    public boolean overlapLayer(final List<ElkNode> nodes) {
        boolean overlapping = false;
        if (nodes.size() < 2) {
            return false;
        }
        for (int i = 0; i < nodes.size(); i++) {
            if (i < nodes.size() - 1) {
                overlapping |= overlap(nodes.get(i), nodes.get(i + 1));
            } else {
                overlapping |= overlap(nodes.get(i), nodes.get(0));

            }
        }
        return overlapping;
    }

    /**
     * The step size which was chosen for the compaction/extension.
     * 
     * @return the compactionStep
     */
    public int getCompactionStep() {
        return compactionStep;
    }

    /**
     * Set the compaction step size which shall be used for the extension/compaction.
     * 
     * @param compactionStep
     *            the compactionStep to set
     */
    public void setCompactionStep(final int compactionStep) {
        this.compactionStep = compactionStep;
    }

    /**
     * The spacing between the nodes which shall be considered while
     * 
     * @return the spacing
     */
    public double getSpacing() {
        return spacing;
    }

    /**
     * @param spacing
     *            the spacing to set
     */
    public void setSpacing(final double spacing) {
        this.spacing = spacing;
    }

    /**
     * Commit the root node of the graph to this class.
     * 
     * @param root
     *            the root to set
     */
    public void setRoot(final ElkNode root) {
        this.root = root;
    }

    /**
     * Receive the saved root node.
     * 
     * @return root
     */
    public ElkNode getRoot() {
        return root;
    }
}