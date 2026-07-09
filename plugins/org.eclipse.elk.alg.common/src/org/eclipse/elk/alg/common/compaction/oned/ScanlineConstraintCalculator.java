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
package org.eclipse.elk.alg.common.compaction.oned;

import java.util.Comparator;
import java.util.List;
import java.util.TreeSet;

import org.eclipse.elk.alg.common.compaction.Scanline;
import org.eclipse.elk.alg.common.compaction.Scanline.EventHandler;

import com.google.common.base.Predicate;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;

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

    // SUPPRESS CHECKSTYLE NEXT 2 VisibilityModifier
    /** The surrounding compactor object. */
    protected OneDimensionalCompactor compactor;
    
    @Override
    public void calculateConstraints(final OneDimensionalCompactor theCompactor) {
        this.compactor = theCompactor;
        
        // consider all nodes
        sweep(n -> true);
    }

    /**
     * Executes a single sweep of the scanline, using {@link CNode}s that fulfill
     * <code>filterFun</code>.
     * 
     * @param filterFun
     *            a predicate determining for a node n whether it should be included
     */
    protected void sweep(final Predicate<CNode> filterFun) {
        
        // add all nodes twice (once for the lower, once for the upper border)
        List<Timestamp> points = Lists.newArrayList();
        for (CNode n : compactor.cGraph.cNodes) {
            if (filterFun.apply(n)) {
                points.add(new Timestamp(n, true));
                points.add(new Timestamp(n, false));
            }
        }
        
        // reset internal state
        constraintsScanlineHandler.reset();
        
        // execute the scanline
        Scanline.execute(points, constraintsScanlineComparator, constraintsScanlineHandler);
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
         * 
         * Further note that it must not happen that two equal elements (equal in terms of the {@code compareTo} method)
         * are to be added to the set at the same time. For the given problem this would happen if two hitboxes
         * overlap vertically and share the same x coordinate. The implementation of {@link TreeSet} would 
         * reject the adding of the second element. See the class's javadoc for further information. 
         */
        private TreeSet<CNode> intervals = Sets.newTreeSet((c1, c2) -> Double.compare(
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
            // System.out.println((p.low ? "insert " : "delete ") + " " + p.node);
            
            if (p.low) {
                insert(p);
            } else {
                delete(p);
            }
        }

        private void insert(final Timestamp p) {
            
            boolean success = intervals.add(p.node);
            if (!success) {
                throw new IllegalStateException("Invalid hitboxes for scanline constraint calculation.");
            }

            if (overlap(p.node, intervals.floor(p.node)) || overlap(p.node, intervals.ceiling(p.node))) {
                System.err.println(p.node + " has overlap.");
            }
            
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
        
    private boolean overlap(final CNode n1, final CNode n2) {
        if (n1 == null || n2 == null || n1 == n2) {
            return false; 
        }
        return CompareFuzzy.le(n1.hitbox.x, n2.hitbox.x + n2.hitbox.width) 
               && CompareFuzzy.le(n2.hitbox.x, n1.hitbox.x + n1.hitbox.width); 
    }
    
    /** This compactor's instance of the scanline handler. */
    private final ConstraintsScanlineHandler constraintsScanlineHandler =
            new ConstraintsScanlineHandler();


}
