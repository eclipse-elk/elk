/*******************************************************************************
 * Copyright (c) 2016 Kiel University and others.
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
package org.eclipse.elk.alg.layered.compaction.oned.algs;

import java.util.Comparator;
import java.util.List;
import java.util.TreeSet;

import org.eclipse.elk.alg.layered.compaction.oned.CNode;
import org.eclipse.elk.alg.layered.compaction.oned.OneDimensionalCompactor;
import org.eclipse.elk.alg.layered.compaction.recthull.Scanline;
import org.eclipse.elk.alg.layered.compaction.recthull.Scanline.EventHandler;

import java.util.function.Function;
import java.util.function.Predicate;
import java.util.ArrayList;

/**
 * Scanline-procedure to calculate the constraint graph in O(nlogn) time, n being the number of
 * rectangles.
 * 
 * The implementation is based on the Section 10.2.2.1.1.1 of the following book: 
 * <ul>
 *   <li>Lengauer, T. (1990). 
 *       Combinatorial algorithms for integrated circuit layout. John Wiley & Sons.
 *   </li>
 * </ul>
 */
public class ScanlineConstraintCalculator implements IConstraintCalculationAlgorithm {

    private static final double EPSILON = .5;
    
    // SUPPRESS CHECKSTYLE NEXT 2 VisibilityModifier
    /** The surrounding compactor object. */
    protected OneDimensionalCompactor compactor;
    
    @Override
    public void calculateConstraints(final OneDimensionalCompactor theCompactor) {
        this.compactor = theCompactor;
        
        // consider all nodes, do not add spacing
        sweep(n -> true, n -> 0.0);
    }

    /**
     * Executes a single sweep of the scanline, using {@link CNode}s that fulfill
     * <code>filterFun</code>.
     * 
     * @param filterFun
     *            a predicate determining for a node n whether it should be included
     * @param spacingFun
     *            a function returning a spacing value s for a node n.
     */
    protected void sweep(
            final Predicate<CNode> filterFun, 
            final Function<CNode, Double> spacingFun) {
        
        blowUpHitboxes(filterFun, spacingFun);
        
        // add all nodes twice (once for the lower, once for the upper border)
        List<Timestamp> points = new ArrayList<>();
        for (CNode n : compactor.cGraph.cNodes) {
            if (filterFun.test(n)) {
                points.add(new Timestamp(n, true));
                points.add(new Timestamp(n, false));
            }
        }
        
        // reset internal state
        constraintsScanlineHandler.reset();
        
        // execute the scanline
        Scanline.execute(points, constraintsScanlineComparator, constraintsScanlineHandler);
        
        normalizeHitboxes(filterFun, spacingFun);
    }
    
    private void blowUpHitboxes(
            final Predicate<CNode> filter, 
            final Function<CNode, Double> spacingFun) {

        for (CNode n : compactor.cGraph.cNodes) {
        
            if (!filter.test(n)) {
                continue;
            }
            
            double spacing = spacingFun.apply(n);
            
            if (spacing > 0) {
                if (!(compactor.direction.isHorizontal() && n.spacingIgnore.up)
                        && !(compactor.direction.isVertical() && n.spacingIgnore.left)) {
                    n.hitbox.y -= Math.max(0, spacing / 2 - EPSILON);
                }
                if (!(compactor.direction.isHorizontal() && n.spacingIgnore.down) 
                        && !(compactor.direction.isVertical() && n.spacingIgnore.right)) {
                    n.hitbox.height += Math.max(0, spacing - 2 * EPSILON);
                }
            }
        }
    }
    
    private void normalizeHitboxes(final Predicate<CNode> filter, 
            final Function<CNode, Double> spacingFun) {
        for (CNode n : compactor.cGraph.cNodes) {
            
            if (!filter.test(n)) {
                continue;
            }
            
            double spacing = spacingFun.apply(n);
            
            if (spacing > 0) {
                if (!(compactor.direction.isHorizontal() && n.spacingIgnore.up)
                        && !(compactor.direction.isVertical() && n.spacingIgnore.left)) {
                    n.hitbox.y += Math.max(0, spacing / 2 - EPSILON);
                }
                if (!(compactor.direction.isHorizontal() && n.spacingIgnore.down) 
                        && !(compactor.direction.isVertical() && n.spacingIgnore.right)) {
                    n.hitbox.height -= spacing - 2 * EPSILON;
                }
            }
        }
    }

    /**
     * A timestamp in the sense of the scanline algorithm to determine the constraints of the
     * compaction graph. The scanline traverses the plane from top (negative infinity) to bottom
     * (positive infinity). Every timestamp (or point) it encounters either represents the upper (y)
     * or lower (y+height) border of a rectangle. Note that y is also denoted as "low" and y+height
     * as "high".
     */
    private class Timestamp {
        private boolean low;
        private CNode node;

        Timestamp(final CNode node, final boolean low) {
            this.node = node;
            this.low = low;
        }
    }
  
    /**
     * Can be used to sort the passed {@link Timestamp}s either based on their 
     * {@link CNode}'s y-coordinate (low) or y-coordinate plus height (high).
     * If two positions are equal, the "high" value is sorted before the "low" value.
     * This implicates that <strong>no</strong> constraint is generated 
     * for element that barely share a common y-value.
     */
    private final Comparator<Timestamp> constraintsScanlineComparator = (p1, p2) ->  {
        
        // chose the right coordinate
        double y1 = p1.node.hitbox.y;
        if (!p1.low) {
            y1 += p1.node.hitbox.height;
        }
        double y2 = p2.node.hitbox.y;
        if (!p2.low) {
            y2 += p2.node.hitbox.height;
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
     * Implements the scanline procedure as discussed by Lengauer.
     */
    private class ConstraintsScanlineHandler implements EventHandler<Timestamp> {

        /**
         * Sorted set of intervals that allows to query for the left and right neighbor of an
         * interval. Sorted by the x coordinate of a {@link CNode}'s hitbox.
         * 
         * Note: the methods to query for neighbor elements (e.g., {@link TreeSet#lower(Object)})
         * are not constant time.
         */
        private TreeSet<CNode> intervals = new TreeSet<>((c1, c2) -> Double.compare(
                c1.hitbox.x + (c1.hitbox.width / 2), c2.hitbox.x + (c2.hitbox.width / 2)));
        /** Candidate array with possible constraints. */
        private CNode[] cand;
        
        /** 
         * Resets the internal data structures.
         */
        public void reset() {
            intervals.clear();
            cand = new CNode[compactor.cGraph.cNodes.size()];
            int index = 0;
            for (CNode n : compactor.cGraph.cNodes) {
                n.id = index++;
            }
        }

        @Override
        public void handle(final Timestamp p) {
            // System.out.println((p.low ? "insert " : "delete ") + " " + p.node.hitbox);
            
            if (p.low) {
                insert(p);
            } else {
                delete(p);
            }
        }

        private void insert(final Timestamp p) {
            intervals.add(p.node);

            cand[p.node.id] = intervals.lower(p.node);

            CNode right = intervals.higher(p.node);
            if (right != null) {
                cand[right.id] = p.node; 
            }
        }

        private void delete(final Timestamp p) {

            CNode left = intervals.lower(p.node);
            if (left != null && left == cand[p.node.id]) {
                // different groups?
                if (left.cGroup != null && left.cGroup != p.node.cGroup) {
                    left.constraints.add(p.node);
                }
            }

            CNode right = intervals.higher(p.node);
            if (right != null && cand[right.id] == p.node) {
                // different groups?
                if (right.cGroup != null && right.cGroup != p.node.cGroup) {
                    p.node.constraints.add(right);
                }
            }

            // we are done with you!
            intervals.remove(p.node);
        }
    }
    
    /** This compactor's instance of the scanline handler. */
    private final ConstraintsScanlineHandler constraintsScanlineHandler =
            new ConstraintsScanlineHandler();


}
