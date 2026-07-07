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

import java.util.Comparator;
import java.util.List;
import java.util.TreeSet;

import org.eclipse.elk.alg.common.compaction.Scanline;
import org.eclipse.elk.alg.common.compaction.Scanline.EventHandler;
import org.eclipse.elk.alg.common.utils.SVGImage;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.google.common.math.DoubleMath;

/**
 * Scanline-procedure to find overlapping nodes.
 */
public class ScanlineOverlapCheck {
    
    /** The function to apply to overlapping nodes. */
    private IOverlapHandler overlapHandler;
    
    /** SVG output for debugging. */
    private SVGImage svg;
    
    /**
     * Constructor.
     * 
     * @param dbg an {@link SVGImage} to write debug output to. null to disable output
     */
    public ScanlineOverlapCheck(final IOverlapHandler oh, final SVGImage dbg) {
        svg = dbg;
        svg.addGroups("n", "o", "e");
        
        overlapHandler = oh;
    }
    
    /**
     * Executes a single sweep of the scanline.
     * 
     * @param nodes the nodes to check
     */
    public void sweep(final List<Node> nodes) {
        
        svg.clearGroup("n");
        svg.clearGroup("o");
        svg.clearGroup("e");
        
        // add all nodes twice (once for the lower, once for the upper border)
        List<Timestamp> points = Lists.newArrayList();
        for (Node n : nodes) {
            points.add(new Timestamp(n, true));
            points.add(new Timestamp(n, false));
            
            svg.g("n").addRect(n.rect, "stroke=\"black\" stroke-width=\"1\" fill=\"none\"");
        }
        
        svg.isave();
        
        OverlapsScanlineHandler overlapsScanlineHandler = new OverlapsScanlineHandler();
        overlapsScanlineHandler.init();
        Scanline.execute(points, overlapsScanlineComparator, overlapsScanlineHandler);
        
        svg.isave();
    }
    
    /**
     * A timestamp in the sense of the scanline algorithm captures where a node begins and ends to determine overlaps.
     * The scanline traverses the plane from top (negative infinity) to bottom
     * (positive infinity). Every timestamp (or point) it encounters either represents the upper (y)
     * or lower (y+height) border of a rectangle. Note that y is also denoted as "low" and y+height
     * as "high".
     */
    private class Timestamp {
        private boolean low;
        private Node node;

        Timestamp(final Node node, final boolean low) {
            this.node = node;
            this.low = low;
        }
    }
    
    private Comparator<Timestamp> overlapsScanlineComparator = (p1, p2) ->  {
        
        // chose the corresponding coordinate
        double y1 = p1.node.rect.y;
        if (!p1.low) {
            y1 += p1.node.rect.height;
        }
        double y2 = p2.node.rect.y;
        if (!p2.low) {
            y2 += p2.node.rect.height;
        }

        // compare them ...
        int cmp = Double.compare(y1, y2);
        
        // ... and if they are equal, sort low->high
        if (cmp == 0) {
            if (!p1.low && p2.low) {
                return -1;
             } else if (!p2.low && p1.low) {
                 return 1;
             }    
        } 
        return cmp;
    };
    
    /**
     * Implements the scanline procedure.
     */
    private class OverlapsScanlineHandler implements EventHandler<Timestamp> {
        
        /** This data structure tracks which nodes are currently on the scanline. */
        private TreeSet<Node> intervalsSortLeft = Sets.newTreeSet((n1, n2) -> {
                    // sort left to right regarding left sides
                    int cmp = Double.compare(n1.rect.x, n2.rect.x);
                    if (cmp != 0) {
                        return cmp;
                    }
                    // if they have the same x-coordinate, sort by center point, which is unique
                    cmp = Double.compare(n1.originalVertex.x, n2.originalVertex.x);
                    if (cmp != 0) {
                        return cmp;
                    }
                    return Double.compare(n1.originalVertex.y, n2.originalVertex.y);
                });
        
        /**
         * Initialize internal data structures.
         * 
         * @param overlapEdges
         */
        public void init() {
            intervalsSortLeft.clear();
        }

        @Override
        public void handle(final Timestamp p) {
            if (p.low) {
                insert(p);
            } else {
                delete(p);
            }
        }
        
        private void insert(final Timestamp p) {
            
            boolean success = intervalsSortLeft.add(p.node);
            if (!success) {
                throw new IllegalStateException("Invalid hitboxes for scanline overlap calculation.");
            }
            
            // find overlaps with rectangles on current scanline and apply the overlapHandler
            boolean overlapsFound = false;
            for (Node other : intervalsSortLeft) { // have to start with first(left side) node because the 
                if (overlap(p.node, other)) {       // last(right side) that does not overlap is not necessarily the
                    overlapHandler.handle(p.node, other); // first candidate to overlap (there could be a wide node).
                    svg.g("o").addRect(p.node.rect, "stroke=\"none\" fill=\"red\" opacity=\"0.18\"");
                    svg.g("o").addRect(other.rect, "stroke=\"none\" fill=\"red\" opacity=\"0.18\"");
                    svg.g("e").addLine(p.node.vertex.x, p.node.vertex.y, other.vertex.x, other.vertex.y, 
                            "stroke=\"blue\"");
                    overlapsFound = true;
                } else {
                    if (overlapsFound) {
                        break; // there cannot be any more overlaps because the set is sorted
                    }
                }
            }
        }

        private void delete(final Timestamp p) {
            intervalsSortLeft.remove(p.node);
        }
    }
    
    private boolean overlap(final Node n1, final Node n2) {
        if (n1 == null || n2 == null || n1 == n2) {
            return false; 
        }
        return DoubleMath.fuzzyCompare(n1.rect.x, n2.rect.x + n2.rect.width, InternalProperties.FUZZINESS) < 0
                && DoubleMath.fuzzyCompare(n2.rect.x, n1.rect.x + n1.rect.width, InternalProperties.FUZZINESS) < 0;
    }
}
