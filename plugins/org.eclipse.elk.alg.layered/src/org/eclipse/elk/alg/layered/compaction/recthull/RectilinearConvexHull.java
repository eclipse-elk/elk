/*******************************************************************************
 * Copyright (c) 2016 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.layered.compaction.recthull;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.ListIterator;

import org.eclipse.elk.alg.layered.compaction.recthull.Point.Quadrant;
import org.eclipse.elk.alg.layered.compaction.recthull.Scanline.EventHandler;
import org.eclipse.elk.core.math.ElkRectangle;

import com.google.common.collect.Lists;

/**
 * The rectilinear convex hull is represented by four staircases.
 * 
 *   q1 | q2
 *  ---------
 *   q4 | q3
 */
public final class RectilinearConvexHull {

    /** This hull's determined points. */
    private List<Point> hull = Lists.newArrayList();

    @SuppressWarnings("unused")
    private Point xMax1 = null, xMax2 = null;
    private Point xMin1 = null, xMin2 = null;
    @SuppressWarnings("unused")
    private Point yMax1 = null, yMax2 = null;
    @SuppressWarnings("unused")
    private Point yMin1 = null, yMin2 = null;
    
    private static final boolean DEBUG = false;

    private RectilinearConvexHull() { }

    /**
     * @param points
     *            a cloud of points for which the rectilinear convex hull should be determined.
     * @return The {@link RectilinearConvexHull} for the passed points.
     */
    public static RectilinearConvexHull of(final Iterable<Point> points) {
        
        RectilinearConvexHull rch = new RectilinearConvexHull();

        for (Point p : points) {
            // equality important since the two extremal points can be equal
            if (rch.xMax1 == null || p.x >= rch.xMax1.x) {
                rch.xMax2 = rch.xMax1;
                rch.xMax1 = p;
            }
            if (rch.xMin1 == null || p.x <= rch.xMin1.x) {
                rch.xMin2 = rch.xMin1;
                rch.xMin1 = p;
            }
            if (rch.yMax1 == null || p.y >= rch.yMax1.y) {
                rch.yMax2 = rch.yMax1;
                rch.yMax1 = p;
            }
            if (rch.yMin1 == null || p.y <= rch.yMin1.y) {
                rch.yMin2 = rch.yMin1;
                rch.yMin1 = p;
            }
        }

        // q1 | q2
        // ---------
        // q4 | q3

        // Q1
        MaximalElementsEventHandler q1 = new MaximalElementsEventHandler(Quadrant.Q1);
        Scanline.execute(points, RIGHT_LOW_FIRST, q1);

        // Q4
        MaximalElementsEventHandler q4 = new MaximalElementsEventHandler(Quadrant.Q4);
        Scanline.execute(points, RIGHT_HIGH_FIRST, q4);

        // Q2
        MaximalElementsEventHandler q2 = new MaximalElementsEventHandler(Quadrant.Q2);
        Scanline.execute(points, LEFT_LOW_FIRST, q2);

        // Q3
        MaximalElementsEventHandler q3 = new MaximalElementsEventHandler(Quadrant.Q3);
        Scanline.execute(points, LEFT_HIGH_FIRST, q3);
        
        if (DEBUG) {
            // print clockwise
            System.out.println("The four stairs (no concaves):");
            System.out.println(q1.points);
            System.out.println(Lists.reverse(q2.points));
            System.out.println(q3.points);
            System.out.println(Lists.reverse(q4.points));
        }
        
        // the scanline algorithm detected all convex corners, 
        // now we add the concave corners to the lists
        addConcaveCorners(q1.points, Quadrant.Q1);
        addConcaveCorners(q2.points, Quadrant.Q2);
        addConcaveCorners(q3.points, Quadrant.Q3);
        addConcaveCorners(q4.points, Quadrant.Q4);

        // print clockwise
        if (DEBUG) {
            System.out.println("The four stairs (with concaves):");
            System.out.println(q1.points);
            System.out.println(Lists.reverse(q2.points));
            System.out.println(q3.points);
            System.out.println(Lists.reverse(q4.points));
        }
        
        // ... and add everything to the hull list
        // (in clockwise order, q1..q4)
        rch.getHull().clear();
        rch.getHull().addAll(q1.points);
        rch.getHull().addAll(Lists.reverse(q2.points));
        rch.getHull().addAll(q3.points);
        rch.getHull().addAll(Lists.reverse(q4.points));

        return rch;
    }
    
    /**
     * @return the points uniquely specifying the rectilinear convex hull. The points are ordered
     *         clock-wise starting with the Q1 staircase (Q1, Q2, Q3, Q4).
     */
    public List<Point> getHull() {
        return hull;
    }

    /**
     * @return the rectilinear convex hull represented by a set of non-overlapping {@link ElkRectangle}s.
     *         The rectangles are determined using a scanline from left-to-right creating a new
     *         rectangle each time a "vertical" segment is detected in the rectilinear convex hull.
     */
    public List<ElkRectangle> splitIntoRectangles() {

        RectangleEventHandler handler = new RectangleEventHandler();
        Scanline.execute(hull, RIGHT_SPECIAL_ORDER, handler);

        // append the lastly queued rectangle
        if (handler.queued != null) {
            handler.rects.add(handler.queued);
        }

        return handler.rects;
    }

    
    /**
     * @param pts
     *            a list of points representing convex points of one staircase of the rectilinear
     *            convex hull.
     * @param q
     *            the {@link Quadrant} for which the points were determined.
     */
    private static void addConcaveCorners(final List<Point> pts, final Quadrant q) {
        ListIterator<Point> pIt = pts.listIterator();
        Point last = pIt.next();
        while (pIt.hasNext()) {
            Point next = pIt.next();
            Point p = new Point(next.x, last.y, q);
            pIt.previous();
            pIt.add(p);
            pIt.next();
            p.convex = false;
            last = next;
        }
    }
    
    private static final Comparator<Point> RIGHT_HIGH_FIRST = (p1, p2) -> {
        if (p1.x == p2.x) {
            return Double.compare(p2.y, p1.y);
        } else {
            return Double.compare(p1.x, p2.x);
        }
    };

    private static final Comparator<Point> RIGHT_LOW_FIRST = (p1, p2) -> {
        if (p1.x == p2.x) {
            return Double.compare(p1.y, p2.y);
        } else {
            return Double.compare(p1.x, p2.x);
        }
    };

    private static final Comparator<Point> LEFT_HIGH_FIRST = (p1, p2) -> {
        if (p1.x == p2.x) {
            return Double.compare(p2.y, p1.y);
        } else {
            return Double.compare(p2.x, p1.x);
        }
    };

    private static final Comparator<Point> LEFT_LOW_FIRST = (p1, p2) -> {
        if (p1.x == p2.x) {
            return Double.compare(p1.y, p2.y);
        } else {
            return Double.compare(p2.x, p1.x);
        }
    };

    /**
     * An {@link EventHandler} that can be passed to the {@link Scanline} algorithm to 
     * determine the rectilinear convex hull.
     */
    private static class MaximalElementsEventHandler implements EventHandler<Point> {

        // SUPPRESS CHECKSTYLE NEXT 10 VisibilityModifier
        private Quadrant quadrant;
        /** The points this handler determined after the execution of the scanline. */
        public List<Point> points = Lists.newArrayList();
        private double maximalY;
        private Comparator<Double> compare;

        private static final Comparator<Double> DBL_CMP = (d1, d2) -> Double.compare(d1, d2);

        /**
         * The handler finds one of the four staircases of the rectilinear convex hull, depending on
         * the specified quadrant. Note however that for this handler to work properly the
         * {@link Scanline} must run left-to-right for quadrants Q1 and Q4 and right-to-left for
         * quadrants Q2 and Q3.
         * See the for specially defined comparators 
         * (such as {@link RectilinearConvexHull#LEFT_HIGH_FIRST) for more details.
         * 
         * @param quadrant
         *            the {@link Quadrant} for which this handler finds points on the hull.
         */
        MaximalElementsEventHandler(final Quadrant quadrant) {
            this.quadrant = quadrant;
            
            switch (quadrant) {
            case Q1: 
            case Q2: 
                compare = Collections.reverseOrder(DBL_CMP);
                maximalY = Double.POSITIVE_INFINITY;
                break;
            case Q3:
            case Q4:
                compare = DBL_CMP;
                maximalY = Double.NEGATIVE_INFINITY;
                break;
            }
        }

        @Override
        public void handle(final Point p) {
            if (compare.compare(p.y, maximalY) > 0) {
                points.add(new Point(p.x, p.y, quadrant));
                maximalY = p.y;
            }
        }
    }

    /*-----------------------------------------------------------------------------------------
     *                                Splitting into rectangles
     *-----------------------------------------------------------------------------------------*/
    

    /**
     * If the x coordinates are the same, Q1,Q4 < Q2,Q3, for Q1,Q2 convex smaller than concave. 
     */
    private static final Comparator<Point> RIGHT_SPECIAL_ORDER = (p1, p2) -> {

        if (p1.x == p2.x) {

            // same quadrant or both either left or right
            if (p1.quadrant == p2.quadrant
                    || Quadrant.isBothLeftOrBothRight(p1.quadrant, p2.quadrant)) {

                // for left quadrants concave before convex,
                // for right quadrants convex before concave
                int val = p1.quadrant.isLeft() ? 1 : -1;
                if (p1.convex && !p2.convex) {
                    return val;
                } else if (!p1.convex && p2.convex) {
                    return -val;
                }
            }

            return Integer.compare(p1.quadrant.ordinal(), p2.quadrant.ordinal());

        } else {
            return Double.compare(p1.x, p2.x);
        }
    };

    /**
     * A handler for the {@link Scanline} that detects vertical segments in the rectilinear convex
     * hull and uses them to create rectangles.
     */
    private class RectangleEventHandler implements EventHandler<Point> {

        private List<ElkRectangle> rects = Lists.newArrayList();

        private Point minY = null;
        private Point maxY = null;

        private double lastX = Math.min(xMin1.x, xMin2.x);

        private ElkRectangle queued = null;
        private Point queuedPnt = null;

        @Override
        public void handle(final Point p) {

            if (queued != null
                    && (p.x != queuedPnt.x 
                        || Quadrant.isOneLeftOneRight(queuedPnt.quadrant, p.quadrant))) {

                rects.add(queued);
                lastX = queued.x + queued.width;
                queued = null;
                queuedPnt = null;
            }

            if (p.quadrant.isUpper()) {
                minY = p;
            } else {
                maxY = p;
            }

            if ((p.quadrant == Quadrant.Q1 && !p.convex) || (p.quadrant == Quadrant.Q2 && p.convex)
                    || (p.quadrant == Quadrant.Q3 && p.convex)
                    || (p.quadrant == Quadrant.Q4 && !p.convex)) {

                // queue a possible rect
                if (minY != null && maxY != null) { // FIXME else?
                    ElkRectangle r = new ElkRectangle(lastX, minY.y, p.x - lastX, maxY.y - minY.y);
                    queued = r;
                    queuedPnt = p;
                }
            }
        }

    }

}
