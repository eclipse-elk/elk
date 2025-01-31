/*******************************************************************************
 * Copyright (c) 2013 - 2022 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.alg.mrtree.p3place;

import java.util.LinkedList;
import org.eclipse.elk.alg.mrtree.TreeLayoutPhases;
import org.eclipse.elk.alg.mrtree.TreeUtil;
import org.eclipse.elk.alg.mrtree.graph.TGraph;
import org.eclipse.elk.alg.mrtree.graph.TNode;
import org.eclipse.elk.alg.mrtree.intermediate.IntermediateProcessorStrategy;
import org.eclipse.elk.alg.mrtree.options.InternalProperties;
import org.eclipse.elk.alg.mrtree.options.MrTreeOptions;
import org.eclipse.elk.core.alg.ILayoutPhase;
import org.eclipse.elk.core.alg.LayoutProcessorConfiguration;
import org.eclipse.elk.core.options.Direction;
import org.eclipse.elk.core.util.IElkProgressMonitor;
import com.google.common.collect.Iterables;

/**
 * The algorithm comes from
 * <ul>
 *   <li> John Q.Walker II, A Node-Positioning Algorithm for General Trees,
 *     <em>Software: Practice and Experience</em> 20(7), pp. 685-705, July 1990.</li>
 * </ul>
 * with some small fixes in the actual code.
 * 
 * <p>
 * This algorithm utilizes two concepts developed in previous positioning algorithms. First is the
 * concept of building subtrees as rigid units. When a node is moved, all of its descendants (if it
 * has any) are also moved--the entire subtree being thus treated as a rigid unit. A general tree is
 * positioned by building it up recursively from its leaves toward its root.
 * </p>
 * 
 * <p>
 * Second is the concept of using two fields for the positioning of each node. These two fields are:
 * <ul>
 *   <li> a preliminary x-coordinate
 *   <li> a modifier field.
 * </ul>
 * Two tree traversals are used to produce the final x-coordinate of a node. The first traversal
 * assigns the preliminary x-coordinate and modifier fields for each node; the second traversal
 * computes the final x-coordinate of each node by summing the node's preliminary x-coordinate with
 * the modifier fields of all of its ancestors. The final y-coordinate of the node is the height of
 * the node's ancestors levels and the height nodes's level and the adjust of the root location.
 * 
 * @author sor
 * @author sgu
 * @author sdo
 */
public class NodePlacer implements ILayoutPhase<TreeLayoutPhases, TGraph> {

    /** intermediate processing configuration. */
    private static final LayoutProcessorConfiguration<TreeLayoutPhases, TGraph> INTERMEDIATE_PROCESSING_CONFIG =
            LayoutProcessorConfiguration.<TreeLayoutPhases, TGraph>create()
                    .addBefore(TreeLayoutPhases.P2_NODE_ORDERING, IntermediateProcessorStrategy.ROOT_PROC)
                    .before(TreeLayoutPhases.P3_NODE_PLACEMENT)
                        .add(IntermediateProcessorStrategy.LEVEL_HEIGHT)
                        .add(IntermediateProcessorStrategy.NEIGHBORS_PROC)
                    .before(TreeLayoutPhases.P4_EDGE_ROUTING)
                        .add(IntermediateProcessorStrategy.DIRECTION_PROC)
                        .add(IntermediateProcessorStrategy.NODE_POSITION_PROC);

    private double spacing;

    /** Determine how to adjust all the nodes with respect to the location of the root. */
    private double xTopAdjustment = 0d; 
    private double yTopAdjustment = 0d;
    
    private Direction direction;

    /**
     * {@inheritDoc}
     */
    @Override
    public LayoutProcessorConfiguration<TreeLayoutPhases, TGraph> getLayoutProcessorConfiguration(final TGraph graph) {
        return INTERMEDIATE_PROCESSING_CONFIG;
    }

    /**
     * {@inheritDoc}
     */
    public void process(final TGraph tGraph, final IElkProgressMonitor progressMonitor) {
        progressMonitor.begin("Processor order nodes", 2);

        /** set the settings according to the user inputs */
        spacing = tGraph.getProperty(MrTreeOptions.SPACING_NODE_NODE).doubleValue();
        direction = tGraph.getProperty(MrTreeOptions.DIRECTION);
        
        // Set Direction to DOWN if its UNDEFINED
        if (direction == Direction.UNDEFINED) {
            direction = Direction.DOWN;
            tGraph.setProperty(MrTreeOptions.DIRECTION, direction);
        }

        /** find the root node of this component */
        LinkedList<TNode> roots = new LinkedList<TNode>();
        for (TNode tNode : tGraph.getNodes()) {
            if (tNode.getProperty(InternalProperties.ROOT)) {
                roots.add(tNode);
            }
        }
        TNode root = roots.getFirst();

        /** Do the preliminary positioning with a postorder walk. */
        firstWalk(root, 0);
        progressMonitor.worked(1);

        /** Do the final positioning with a preorder walk. */
        secondWalk(root, yTopAdjustment - (root.getProperty(InternalProperties.LEVELHEIGHT) / 2), xTopAdjustment);
        progressMonitor.worked(1);

        progressMonitor.done();
    }

    /**
     * In this first postorder walk, every node of the tree is assigned a preliminary x-coordinate
     * (held in property PRELIM). In addition, internal nodes are given modifiers, which will be
     * used to move their offspring to the right (held in property MODIFIER).
     * 
     * @param cN
     *            the root level of the tree
     * @param level
     *            the index of the passed level
     */
    private void firstWalk(final TNode cN, final int level) {
        cN.setProperty(InternalProperties.MODIFIER, 0d);
        TNode lS = cN.getProperty(InternalProperties.LEFTSIBLING);

        if (cN.isLeaf()) {
            if (lS != null) {
                /**
                 * Determine the preliminary x-coordinate based on: the preliminary x-coordinate of
                 * the left sibling, the separation between sibling nodes, and tHe mean size of left
                 * sibling and current node.
                 */
                double p = lS.getProperty(InternalProperties.PRELIM) + spacing + meanNodeWidth(lS, cN);
                cN.setProperty(InternalProperties.PRELIM, p);
            } else {
                /** No sibling on the left to worry about. */
                cN.setProperty(InternalProperties.PRELIM, 0d);
            }
        } else {
            /**
             * This Node is not a leaf, so call this procedure recursively for each of its
             * offspring.
             */
            for (TNode child : cN.getChildren()) {
                firstWalk(child, level + 1);
            }

            /**
             * Set the prelim and modifer for this node by determine the midpoint of its offsprings
             * and the middle node size of the node and its left sibling
             */
            TNode lM = Iterables.getFirst(cN.getChildren(), null);
            TNode rM = Iterables.getLast(cN.getChildren(), null);
            double midPoint = (rM.getProperty(InternalProperties.PRELIM) + lM
                    .getProperty(InternalProperties.PRELIM)) / 2f;

            if (lS != null) {
                /** This Node has a left sibling so its offsprings must be shifted to the right */
                double p = lS.getProperty(InternalProperties.PRELIM) + spacing + meanNodeWidth(lS, cN);
                cN.setProperty(InternalProperties.PRELIM, p);
                cN.setProperty(InternalProperties.MODIFIER, cN.getProperty(InternalProperties.PRELIM) - midPoint);
                /** shift the offsprings of this node to the right */
                apportion(cN, level);
            } else {
                /** No sibling on the left to worry about. */
                cN.setProperty(InternalProperties.PRELIM, midPoint);
            }
        }
    }

    /**
     * This method cleans up the positioning of small sibling subtrees, thus fixing the
     * "left-to-right gluing" problem evident in earlier algorithms. When moving a new subtree
     * farther and farther to the right, gaps may open up among smaller subtrees that were
     * previously sandwiched between larger subtrees. Thus, when moving the new, larger subtree to
     * the right, the distance it is moved is also apportioned to smaller, interior subtrees,
     * creating a pleasing aesthetic placement.
     * 
     * @param cN
     *            the root of the subtree
     * @param level
     *            the level of the root in the global tree
     */
    private void apportion(final TNode cN, final int level) {
        /** Initialize the leftmost and neighbor corresponding to the root of the subtree */
        TNode leftmost = Iterables.getFirst(cN.getChildren(), null);
        TNode neighbor = leftmost != null ? leftmost.getProperty(InternalProperties.LEFTNEIGHBOR) : null;
        int compareDepth = 1;
        /**
         * until this node and the neighbor to the left have nodes in the current level we have to
         * shift the current subtree
         */
        while (leftmost != null && neighbor != null) {
            /** Compute the location of Leftmost and where it should be with respect to Neighbor. */
            double leftModSum = 0;
            double rightModSum = 0;
            TNode ancestorLeftmost = leftmost;
            TNode ancestorNeighbor = neighbor;
            /** sum the modifiers of all ancestors according to the current level */
            for (int i = 0; i < compareDepth; i++) {
                ancestorLeftmost = ancestorLeftmost.getParent();
                ancestorNeighbor = ancestorNeighbor.getParent();
                rightModSum += ancestorLeftmost.getProperty(InternalProperties.MODIFIER);
                leftModSum += ancestorNeighbor.getProperty(InternalProperties.MODIFIER);
            }
            /**
             * Find the MoveDistance, and apply it to Node's subtree. Add appropriate portions to
             * smaller interior subtrees.
             */
            double prN = neighbor.getProperty(InternalProperties.PRELIM);
            double prL = leftmost.getProperty(InternalProperties.PRELIM);
            double mean = meanNodeWidth(leftmost, neighbor);
            double moveDistance = prN + leftModSum + spacing + mean - prL - rightModSum;

            if (0 < moveDistance) {
                /** Count interior sibling subtrees in LeftSiblings */
                TNode leftSibling = cN;
                int leftSiblings = 0;
                while (leftSibling != null && leftSibling != ancestorNeighbor) {
                    leftSiblings++;
                    leftSibling = leftSibling.getProperty(InternalProperties.LEFTSIBLING);
                }
                /** Apply portions to appropriate left sibling subtrees. */
                if (leftSibling != null) {
                    double portion = moveDistance / (double) leftSiblings;
                    leftSibling = cN;
                    while (leftSibling != ancestorNeighbor) {
                        double newPr = leftSibling.getProperty(InternalProperties.PRELIM) + moveDistance;
                        leftSibling.setProperty(InternalProperties.PRELIM, newPr);
                        double newMod = leftSibling.getProperty(InternalProperties.MODIFIER) + moveDistance;
                        leftSibling.setProperty(InternalProperties.MODIFIER, newMod);
                        moveDistance -= portion;
                        leftSibling = leftSibling.getProperty(InternalProperties.LEFTSIBLING);
                    }
                } else {
                    /**
                     * Don't need to move anything--it needs to be done by an ancestor because
                     * AncestorNeighbor and AncestorLeftmost are not siblings of each other.
                     */
                    return;
                }
            }
            /**
             * Determine the leftmost descendant of Node at the next lower level to compare its
             * positioning against that of its Neighbor.
             */
            compareDepth++;
            if (leftmost.isLeaf()) {
                leftmost = TreeUtil.getLeftMost(cN.getChildren(), compareDepth);
            } else {
                leftmost = Iterables.getFirst(leftmost.getChildren(), null);
            }
            neighbor = leftmost != null ? leftmost.getProperty(InternalProperties.LEFTNEIGHBOR) : null;
        }
    }

    /**
     * This function returns the mean width of the two passed nodes. It adds the width of the right
     * half of left hand node to the left half of right hand node. If all nodes are the same width,
     * this is a trivial calculation.
     * 
     * @param leftNode
     *            the left hand node
     * @param rightNode
     *            the right hand node
     * @return the sum of the width
     */
    private double meanNodeWidth(final TNode leftNode, final TNode rightNode) {
        double nodeWidth = 0d;
        if (leftNode != null) {
            if (direction.isVertical()) {
                nodeWidth += leftNode.getSize().x / 2d;
            } else {
                nodeWidth += leftNode.getSize().y / 2d;
            }
        }
        if (rightNode != null) {
            if (direction.isVertical()) {
                nodeWidth += rightNode.getSize().x / 2d;
            } else {
                nodeWidth += rightNode.getSize().y / 2d;
            }
        }
        return nodeWidth;
    }

    /**
     * During a second preorder walk, each node is given a final x-coordinate by summing its
     * preliminary x-coordinate and the modifiers of all the node's ancestors. The y-coordinate
     * depends on the height of the node's ancestors levels. If the actual position of an interior
     * node is right of its preliminary place, the subtree rooted at the node must be moved right to
     * center the sons around the father. Rather than immediately readjust all the nodes in the
     * subtree, each node remembers the distance to the provisional place in a modifier field
     * (property MODIFIER). In this second pass down the tree, modifiers are accumulated and applied
     * to every node.
     * 
     * @param tNode
     *            the root of the tree
     * @param yCoor
     *            the y coordinate of previous level
     * @param modsum
     *            the modifiers of all the node's ancestors
     */
    private void secondWalk(final TNode tNode, final double yCoor, final double modsum) {
        if (tNode != null) {
            // The x-position of the node is the sum of its prev x-coordinate and the modifiers of
            // all the node's ancestors and the adjust of the root location.
            double xTemp = tNode.getProperty(InternalProperties.PRELIM) + modsum;
            // The y-position of the node is the height of the node's ancestors levels and the
            // height nodes's level and the adjust of the root location.
            double yTemp = yCoor + (tNode.getProperty(InternalProperties.LEVELHEIGHT) / 2);
            // We do not check to see that xTemp and yTemp are of the proper size, because the
            // framework will take care of this.
            tNode.setProperty(InternalProperties.XCOOR, (int) Math.round(xTemp));
            tNode.setProperty(InternalProperties.YCOOR, (int) Math.round(yTemp));
            // Apply the modifier value for this node to all its offspring.
            if (!tNode.isLeaf()) {
                // This node got offsprings so we will step a level down and take care of them.
                secondWalk(Iterables.getFirst(tNode.getChildren(), null),
                        yCoor + tNode.getProperty(InternalProperties.LEVELHEIGHT) + spacing,
                        modsum + tNode.getProperty(InternalProperties.MODIFIER));
            }
            // Go ahead with the sibling to the right. This is a dfs so we just layout the current
            // subtree.
            if (tNode.getProperty(InternalProperties.RIGHTSIBLING) != null) {
                secondWalk(tNode.getProperty(InternalProperties.RIGHTSIBLING), yCoor, modsum);
            }
        }
    }
}
