/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.common.spore;

import org.eclipse.elk.alg.common.Tree;
import org.eclipse.elk.alg.common.utils.SVGImage;
import org.eclipse.elk.core.math.ElkRectangle;
import org.eclipse.elk.core.math.KVector;

import com.google.common.math.DoubleMath;

/**
 * This class implements a compaction method that expects a spanning tree of a graph.
 * While traversing the spanning tree in an depth-first manner the nodes of each subtree are shifted as one unit
 * towards the parent node of that subtree. 
 */
public final class DepthFirstCompaction {
    /** The root of the spanning tree. */
    private static Tree<Node> root;
    /** Restricts the translation of nodes to orthogonal directions. */
    private static boolean orthogonalCompaction = false;
    
    private static SVGImage svg;
    
    /** Hidden constructor. */
    private DepthFirstCompaction() { };

    /**
     * Executes the compaction of the graph represented by the spanning tree.
     * @param tree the spanning tree of a graph
     * @param orthogonal whether to restrict the movement to orthogonal translations
     * @param debugOutputFile file name for debug SVG. Debug output will be deactivated if this is null.
     */
    public static void compact(final Tree<Node> tree, final boolean orthogonal, 
            final String debugOutputFile) {
        svg = new SVGImage(debugOutputFile);
        orthogonalCompaction = orthogonal;
        root = tree;
        debugOut(tree, null);
        compactTree(root);
    }
    
    /**
     * Executes the compaction of the graph represented by the spanning tree.
     * @param tree the spanning tree of a graph
     * @param orthogonal whether to restrict the movement to orthogonal translations
     */
    public static void compact(final Tree<Node> tree, final boolean orthogonal) {
        compact(tree, orthogonal, null);
    }

    /**
     * Recursive function to compact a tree depth-first.
     * @param tree root of a tree
     */
    private static void compactTree(final Tree<Node> tree) {
        // first compact the children of the current node
        tree.children.forEach(DepthFirstCompaction::compactTree);
        
        // remove underlap between root and its children
        for (Tree<Node> child : tree.children) {
            
            // find the direction and distance the subtree has to be moved to close the gap to the root node
            KVector compactionVector = tree.node.vertex.clone().sub(child.node.vertex);
            
            if (orthogonalCompaction) {
                ElkRectangle rt = tree.node.rect;
                ElkRectangle rc = child.node.rect;
                // restrict the compaction vector to its larger coordinate
                if (Math.abs(compactionVector.x) >= Math.abs(compactionVector.y)) {
                    compactionVector.y = 0;
                    // scale compaction vector to orthogonal distance between tree and child
                    if (rc.y + rc.height > rt.y && rc.y < rt.y + rt.height) { // overlap in vertical dimension
                        compactionVector.scaleToLength(Math.max(rt.x - (rc.x + rc.width), rc.x - (rt.x + rt.width)));
                    }
                } else {
                    compactionVector.x = 0;
                 // scale compaction vector to orthogonal distance between tree and child
                    if (rc.x + rc.width > rt.x && rc.x < rt.x + rt.width) { // overlap in horizontal dimension
                        compactionVector.scaleToLength(Math.max(rt.y - (rc.y + rc.height), rc.y - (rt.y + rt.height)));
                    }
                }
            } else {
                compactionVector.scaleToLength(tree.node.underlap(child.node));
            }
            
            // determine the maximum distance the subtree can be moved without collision, i.e. the minimum underlap
            double minUnderlap = compactionVector.length();
            
            // find minimum underlap between any node in the child and any other node in the rest of the 
            // tree in compaction direction
            minUnderlap = getMinUnderlap(root, child, minUnderlap, compactionVector);
            
            // use the minimum underlap to move whole subtree
            compactionVector.scaleToLength(minUnderlap);
            translateSubtree(child, compactionVector);
            
            debugOut(tree, child);
        }
    }
    
    /**
     * Computes the minimum distance between nodes in a tree and nodes in a subtree of the tree.
     * 
     * @param tree the entire tree
     * @param child a subtree
     * @param currentMinUnderlap an initial value
     * @param compactionVector a KVector that determines the direction to calculate distances in
     * @return the minimum distance or underlap
     */
    private static double getMinUnderlap(final Tree<Node> tree, final Tree<Node> child, 
            final double currentMinUnderlap, final KVector compactionVector) {
        
        double minUnderlap = Math.min(currentMinUnderlap, 
                minUnderlapWithSubtree(tree.node, child, currentMinUnderlap, compactionVector));
        
        for (Tree<Node> c : tree.children) {
            if (c != child) {
                minUnderlap = Math.min(minUnderlap, getMinUnderlap(c, child, minUnderlap, compactionVector));
            }
        }
        
        return minUnderlap;
    }
    
    private static double minUnderlapWithSubtree(final Node r, final Tree<Node> tree, final double currentMinUnderlap,
            final KVector compactionVector) {
        
        double minUnderlap = currentMinUnderlap;
        for (Tree<Node> child : tree.children) {
            Node c = child.node;
            if (r.touches(c)) {
                // if they already touch, find out if they would be moved into each other
                if ((DoubleMath.fuzzyCompare(c.rect.x, r.rect.x + r.rect.width, 
                        InternalProperties.FUZZINESS) == 0
                    && compactionVector.x < 0)
                    || (DoubleMath.fuzzyCompare(c.rect.x + c.rect.width, r.rect.x, 
                            InternalProperties.FUZZINESS) == 0
                        && compactionVector.x > 0)
                    || (DoubleMath.fuzzyCompare(c.rect.y, r.rect.y + r.rect.height, 
                            InternalProperties.FUZZINESS) == 0
                        && compactionVector.y < 0)
                    || (DoubleMath.fuzzyCompare(c.rect.y + c.rect.height, r.rect.y, 
                            InternalProperties.FUZZINESS) == 0
                        && compactionVector.y > 0)) {
                    minUnderlap = 0;
                    break;
                }
            } else {
                // update minUnderlap
                minUnderlap = Math.min(minUnderlap, r.distance(c, compactionVector));
            }
            
            minUnderlap = Math.min(minUnderlap, minUnderlapWithSubtree(r, child, minUnderlap, compactionVector));
            
        }
        
        return minUnderlap;
    }
    
    private static void translateSubtree(final Tree<Node> tree, final KVector compactionVector) {
        tree.node.translate(compactionVector);
        tree.children.forEach(c -> translateSubtree(c, compactionVector));
    }

    
    //--- debug utils --------------------------------------
    
    private static void drawTree(final Tree<Node> t, final SVGImage svgImage, final Tree<Node> mark) {
        svgImage.g("rects").addRect(t.node.rect, "fill=\"none\" stroke=\"black\"");
        // SUPPRESS CHECKSTYLE NEXT 1 MagicNumber
        svgImage.g("centers").addCircle(t.node.vertex.x, t.node.vertex.y, 6, "fill=\"black\"");
        t.children.forEach(c -> {
            if (c == mark) {
                svgImage.g("edges").addLine(t.node.vertex.x, t.node.vertex.y, c.node.vertex.x, c.node.vertex.y, 
                        "stroke=\"red\" stroke-width=\"5\" stroke-dasharray=\"8,8\"");
            } else {
                svgImage.g("edges").addLine(t.node.vertex.x, t.node.vertex.y, c.node.vertex.x, c.node.vertex.y, 
                        "stroke=\"blue\" stroke-width=\"2\"");
            }
            drawTree(c, svgImage, mark);
        });
    }

    private static void debugOut(final Tree<Node> tree, final Tree<Node> mark) {
        svg.clear();
        svg.addGroups("rects", "root", "edges", "centers", "marks");
        svg.g("root").addRect(root.node.rect, "fill=\"#C4E3F3\" stroke=\"black\"");
        drawTree(root, svg, mark);
        svg.isave();
    }
}
