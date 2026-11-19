/*******************************************************************************
 * Copyright (c) 2009, 2015 Kiel University and others.
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
package org.eclipse.elk.core.math;

import java.util.Iterator;
import java.util.ListIterator;

import com.google.common.math.DoubleMath;

/**
 * Mathematics utility class for the Eclipse Layout Kernel.
 * 
 * @author msp
 */
public final class ElkMath {

    /**
     * Hidden constructor to avoid instantiation.
     */
    private ElkMath() {
    }

    /** table of precomputed factorial values. */
    private static final long[] FACT_TABLE = { 1L, 1L, 2L, 6L, 24L, 120L, 720L, 5040L, 40320L,
            362880L, 3628800L, 39916800L, 479001600L, 6227020800L, 87178291200L, 1307674368000L,
            20922789888000L, 355687428096000L, 6402373705728000L, 121645100408832000L,
            2432902008176640000L };

    /**
     * The factorial of an integer x as long value. If x is negative or greater then 20 then
     * IllegalArgumentException. This method always returns the exact value for x between 0 and 20.
     * 
     * @param x
     *            an integer
     * @return the factorial of x
     * @throws IllegalArgumentException
     *             if x<0 or x>20
     */
    public static long factl(final int x) {
        // IllegalArgumentException when x not between 0 and 20
        if (x < 0 || x >= FACT_TABLE.length) {
            throw new IllegalArgumentException("The input must be between 0 and "
                    + FACT_TABLE.length);
        }
        // return the appropriate value from FACT_TABLE with index x
        return FACT_TABLE[x];

    }

    /**
     * The factorial of an integer x as double value. If x is negative then
     * IllegalArgumentException. If x>26 then the result is Infinity. This method returns the exact
     * value for small input values, and uses Stirling's approximation for large input values.
     * 
     * @param x
     *            an integer
     * @return the factorial of x
     * @throws IllegalArgumentException
     *             if x<0
     */
    public static double factd(final int x) {
        // IllegalArgumentException when x < 0
        if (x < 0) {
            throw new IllegalArgumentException("The input must be positive");
        } else if (x < FACT_TABLE.length) {
            return FACT_TABLE[x];
        } else {
            return Math.sqrt(2.0 * Math.PI * x) * (powf(x, x) / powd(Math.E, x));
        }
    }

    /**
     * The binomial coefficient of integers n and k as long value. If n is not positive or k is not
     * between 0 and n the result is IllegalArgumentException. This method always returns the exact
     * value, but may take very long for large input values.
     * 
     * @param n
     *            the upper integer
     * @param k
     *            the lower integer
     * @return n choose k
     * @throws IllegalArgumentException
     *             if n < 0 or n < 0 or k > n
     */
    public static long binomiall(final int n, final int k) {
        if (n < 0 || k < 0) {
            throw new IllegalArgumentException("k and n must be positive");
        } else if (k > n) {
            throw new IllegalArgumentException("k must be smaller than n");
        } else if (k == 0 || k == n) {
            return 1;
        } else if (n == 0) {
            return 0;
        } else if (n < FACT_TABLE.length) {
            return factl(n) / (factl(k) * factl(n - k));
        } else {
            return binomiall(n - 1, k - 1) + binomiall(n - 1, k);
        }
    }

    /**
     * The binomial coefficient of integers n and k as double value. If n is not positive or k is
     * not between 0 and n the result is IllegalArgumentException. This method returns the exact
     * value for small input values, and uses an approximation for large input values.
     * 
     * @param n
     *            the upper integer
     * @param k
     *            the lower integer
     * @return n choose k
     * @throws IllegalArgumentException
     *             if n < 0 or n < 0 or k > n
     */
    public static double binomiald(final int n, final int k) {
        if (n < 0 || k < 0) {
            throw new IllegalArgumentException("k and n must be positive");
        } else if (k > n) {
            throw new IllegalArgumentException("k must be smaller than n");
        } else if (k == 0 || k == n) {
            return 1;
        } else if (n == 0) {
            return 0;
        } else {
            return factd(n) / (factd(k) * factd(n - k));
        }
    }

    /**
     * The first argument raised to the power of the second argument.
     * 
     * @param a
     *            the base
     * @param b
     *            the exponent
     * @return a to the power of b
     */
    public static double powd(final double a, final int b) {
        double result = 1.0;
        double base = a;
        int exp = (b >= 0 ? b : -b);
        while (exp > 0) {
            if (exp % 2 == 0) {
                base *= base;
                exp /= 2;
            } else {
                result *= base;
                exp -= 1;
            }
        }
        if (b < 0) {
            return 1.0 / result;
        } else {
            return result;
        }
    }

    /**
     * The first argument raised to the power of the second argument.
     * 
     * @param a
     *            the base
     * @param b
     *            the exponent
     * @return a to the power of b
     */
    public static float powf(final float a, final int b) {
        float result = 1.0f;
        float base = a;
        int exp = (b >= 0 ? b : -b);
        while (exp > 0) {
            if (exp % 2 == 0) {
                base *= base;
                exp /= 2;
            } else {
                result *= base;
                exp -= 1;
            }
        }
        if (b < 0) {
            return 1.0f / result;
        } else {
            return result;
        }
    }

    /**
     * Compute a number of approximation points on the Bezier curve defined by the given control
     * points. The degree of the curve is derived from the number of control points. The array of
     * resulting curve points includes the target point, but does not include the source point of
     * the curve.
     * 
     * @param controlPoints
     *            the control points
     * @param resultSize
     *            number of returned curve points
     * @return points on the curve defined by the given control points
     */
    public static KVector[] approximateBezierSegment(final int resultSize,
            final KVector... controlPoints) {
        if (resultSize <= 0) {
            return new KVector[0];
        }
        KVector[] result = new KVector[resultSize];
        double dt = (1.0 / resultSize);
        double t = 0;
        for (int i = 0; i < resultSize; i++) {
            t += dt;
            result[i] = getPointOnBezierSegment(t, controlPoints);
        }
        return result;
    }

    /**
     * Compute a number of approximation points on the Bezier curve defined by the given control
     * points. The degree of the curve is derived from the number of control points. The array of
     * resulting curve points includes the target point, but does not include the source point of
     * the curve. The number of approximation points is derived from the given control points.
     * 
     * @param controlPoints
     *            the control points of the curve
     * @return points on the curve defined by the given control points
     */
    public static KVector[] approximateBezierSegment(final KVector... controlPoints) {
        // The number of approximation points simply equals the number of control points.
        // Although there might be more accurate approximations, this approach is the fastest.
        int approximationCount = controlPoints.length + 1;
        return approximateBezierSegment(approximationCount, controlPoints);
    }
    
    /**
     * Compute a point on the Bezier curve defined by the given control points. The value {@code t}
     * determines the position of the returned point.
     * 
     * @param t
     *            a value between 0 and 1, where 0 is the start point of the curve, and 1 is the
     *            end point
     * @param controlPoints
     *            the control points of the curve
     * @return the point at position {@code t}
     */
    public static KVector getPointOnBezierSegment(final double t, final KVector... controlPoints) {
        int n = controlPoints.length - 1;
        double px = 0;
        double py = 0;
        for (int j = 0; j <= n; j++) {
            KVector p = controlPoints[j];
            double factor = binomiald(n, j) * powd(1 - t, n - j) * powd(t, j);
            px += p.x * factor;
            py += p.y * factor;
        }
        return new KVector(px, py);
    }

    /**
     * Compute an approximation for the spline that is defined by the given control points. The
     * control points are interpreted as a series of cubic Bezier curves.
     * 
     * <p><em>Note:</em> As a more powerful alternative, you might consider using
     * {@link BezierSpline} instead.</p>
     * 
     * @param controlPoints
     *            control points of a piecewise cubic spline
     * @return a vector chain that approximates the spline
     */
    public static KVectorChain approximateBezierSpline(final KVectorChain controlPoints) {
        int ctrlPtCount = controlPoints.size();
        KVectorChain spline = new KVectorChain();
        ListIterator<KVector> controlIter = controlPoints.listIterator();
        KVector currentPoint = controlIter.next();
        spline.add(currentPoint);
        while (controlIter.hasNext()) {
            int remainingPoints = ctrlPtCount - controlIter.nextIndex();
            if (remainingPoints == 1) {
                spline.add(controlIter.next());
            } else if (remainingPoints == 2) {
                // calculate a quadratic bezier curve
                spline.addAll(approximateBezierSegment(currentPoint, controlIter.next(), controlIter.next()));
            } else {
                // calculate a cubic bezier curve
                KVector control1 = controlIter.next();
                KVector control2 = controlIter.next();
                KVector nextPoint = controlIter.next();
                spline.addAll(approximateBezierSegment(currentPoint, control1, control2, nextPoint));
                currentPoint = nextPoint;
            }
        }
        return spline;
    }

    /** degree of splines equation to find roots. */
    private static final int W_DEGREE = 5;

    /**
     * Calculate the distance from a cubic spline curve to the point {@code needle}.
     * 
     * @param start
     *            starting point
     * @param c1
     *            control point 1
     * @param c2
     *            control point 2
     * @param end
     *            end point
     * @param needle
     *            point to look for
     * @return distance from needle to curve
     */
    public static double distanceFromBezierSegment(final KVector start, final KVector c1,
            final KVector c2, final KVector end, final KVector needle) {
        double[] tCandidate = new double[W_DEGREE]; // possible roots
        KVector[] v = { start, c1, c2, end };

        // convert problem to 5th-degree Bezier form
        KVector[] w = convertToBezierForm(v, needle);

        // Find all possible roots of 5th-degree equation
        int nSolutions = findRoots(w, W_DEGREE, tCandidate, 0);

        // Compare distances of P5 to all candidates, and to t=0, and t=1
        // Check distance to beginning of curve, where t = 0
        double minDistance = needle.distance(start);
        double t = 0.0;

        // Find distances for candidate points
        for (int i = 0; i < nSolutions; i++) {
            KVector p = bezier(v, DEGREE, tCandidate[i], null, null);
            double distance = needle.distance(p);
            if (distance < minDistance) {
                minDistance = distance;
                t = tCandidate[i];
            }
        }

        // Finally, look at distance to end point, where t = 1.0
        double distance = needle.distance(end);
        if (distance < minDistance) {
            t = 1.0;
        }

        // Return the point on the curve at parameter value t
        KVector pn = new KVector(bezier(v, DEGREE, t, null, null));
        return Math.sqrt(pn.distance(needle));
    }

    /** cubic Bezier curves. */
    private static final int DEGREE = 3;
    /** precomputed "z" for cubics. */
    private static final double[][] CUBIC_Z = { { 1.0, 0.6, 0.3, 0.1 }, { 0.4, 0.6, 0.6, 0.4 },
            { 0.1, 0.3, 0.6, 1.0 }, };

    /**
     * Given a point and a Bezier curve, generate a 5th-degree Bezier-format equation whose solution
     * finds the point on the curve nearest the user-defined point.
     */
    private static KVector[] convertToBezierForm(final KVector[] v, final KVector pa) {
        KVector[] c = new KVector[DEGREE + 1]; // v(i) - pa
        KVector[] d = new KVector[DEGREE]; // v(i+1) - v(i)
        double[][] cdTable = new double[DEGREE][DEGREE + 1]; // dot product of c, d
        KVector[] w = new KVector[W_DEGREE + 1]; // ctl pts of 5th-degree curve

        // Determine the c's -- these are vectors created by subtracting
        // point pa from each of the control points
        for (int i = 0; i <= DEGREE; i++) {
            c[i] = new KVector(v[i].x - pa.x, v[i].y - pa.y);
        }

        // Determine the d's -- these are vectors created by subtracting
        // each control point from the next
        double s = DEGREE;
        for (int i = 0; i <= DEGREE - 1; i++) {
            d[i] = new KVector(s * (v[i + 1].x - v[i].x), s * (v[i + 1].y - v[i].y));
        }

        // Create the c,d table -- this is a table of dot products of the
        // c's and d's */
        for (int row = 0; row <= DEGREE - 1; row++) {
            for (int column = 0; column <= DEGREE; column++) {
                cdTable[row][column] = (d[row].x * c[column].x) + (d[row].y * c[column].y);
            }
        }

        // Now, apply the z's to the dot products, on the skew diagonal
        // Also, set up the x-values, making these "points"
        for (int i = 0; i <= W_DEGREE; i++) {
            w[i] = new KVector((double) (i) / W_DEGREE, 0.0);
        }

        int n = DEGREE;
        int m = DEGREE - 1;
        for (int k = 0; k <= n + m; k++) {
            int lb = Math.max(0, k - m);
            int ub = Math.min(k, n);
            for (int i = lb; i <= ub; i++) {
                int j = k - i;
                w[i + j].y = w[i + j].y + cdTable[j][i] * CUBIC_Z[j][i];
            }
        }

        return w;
    }

    /** maximum depth for recursion. */
    private static final int MAXDEPTH = 64;

    /**
     * Given a 5th-degree equation in Bernstein-Bezier form, find all of the roots in the interval
     * [0, 1].
     * 
     * @return the number of roots found.
     */
    private static int findRoots(final KVector[] w, final int degree, final double[] t,
            final int depth) {
        switch (crossingCount(w, degree)) {
        case 0: // No solutions here
            return 0;
        case 1: // Unique solution
            // Stop recursion when the tree is deep enough
            // if deep enough, return 1 solution at midpoint
            if (depth >= MAXDEPTH) {
                t[0] = (w[0].x + w[W_DEGREE].x) / 2.0;
                return 1;
            }
            if (controlPolygonFlatEnough(w, degree)) {
                t[0] = computeXIntercept(w, degree);
                return 1;
            }
            break;
        default: // nothing
        }

        // Otherwise, solve recursively after subdividing control polygon
        KVector[] left = new KVector[W_DEGREE + 1]; // New left and right
        KVector[] right = new KVector[W_DEGREE + 1]; // control polygons
        double[] leftT = new double[W_DEGREE + 1]; // Solutions from kids
        double[] rightT = new double[W_DEGREE + 1];

        // start in the middle of the bezier curve, t=0.5
        bezier(w, degree, 1.0 / 2, left, right);
        int leftCount = findRoots(left, degree, leftT, depth + 1);
        int rightCount = findRoots(right, degree, rightT, depth + 1);

        // Gather solutions together
        for (int i = 0; i < leftCount; i++) {
            t[i] = leftT[i];
        }
        for (int i = 0; i < rightCount; i++) {
            t[i + leftCount] = rightT[i];
        }

        // Send back total number of solutions */
        return leftCount + rightCount;
    }

    /** Flatness. */
    private static final double EPSILON = 1.0 * Math.pow(2, -MAXDEPTH - 1);

    /**
     * Check if the control polygon of a Bezier curve is flat enough for recursive subdivision to
     * bottom out.
     */
    private static boolean controlPolygonFlatEnough(final KVector[] v, final int degree) {

        // Find the perpendicular distance from each interior control point to
        // line connecting v[0] and v[degree]

        // Derive the implicit equation for line connecting first
        // and last control points
        double a = v[0].y - v[degree].y;
        double b = v[degree].x - v[0].x;
        double c = v[0].x * v[degree].y - v[degree].x * v[0].y;

        double abSquared = (a * a) + (b * b);
        double[] distance = new double[degree + 1]; // Distances from pts to line

        for (int i = 1; i < degree; i++) {
            // Compute distance from each of the points to that line
            distance[i] = a * v[i].x + b * v[i].y + c;
            if (distance[i] > 0.0) {
                distance[i] = (distance[i] * distance[i]) / abSquared;
            }
            if (distance[i] < 0.0) {
                distance[i] = -((distance[i] * distance[i]) / abSquared);
            }
        }

        // Find the largest distance
        double maxDistanceAbove = 0.0;
        double maxDistanceBelow = 0.0;
        for (int i = 1; i < degree; i++) {
            if (distance[i] < 0.0) {
                maxDistanceBelow = Math.min(maxDistanceBelow, distance[i]);
            }
            if (distance[i] > 0.0) {
                maxDistanceAbove = Math.max(maxDistanceAbove, distance[i]);
            }
        }

        // Implicit equation for zero line
        double a1 = 0.0;
        double b1 = 1.0;
        double c1 = 0.0;

        // Implicit equation for "above" line
        double a2 = a;
        double b2 = b;
        double c2 = c + maxDistanceAbove;

        double det = a1 * b2 - a2 * b1;
        double dInv = 1.0 / det;

        double intercept1 = (b1 * c2 - b2 * c1) * dInv;

        // Implicit equation for "below" line
        a2 = a;
        b2 = b;
        c2 = c + maxDistanceBelow;

        det = a1 * b2 - a2 * b1;
        dInv = 1.0 / det;

        double intercept2 = (b1 * c2 - b2 * c1) * dInv;

        // Compute intercepts of bounding box
        double leftIntercept = Math.min(intercept1, intercept2);
        double rightIntercept = Math.max(intercept1, intercept2);

        double error = (rightIntercept - leftIntercept) / 2;

        return error < EPSILON;
    }

    /**
     * Compute intersection of chord from first control point to last with 0-axis.
     */
    private static double computeXIntercept(final KVector[] v, final int degree) {
        double xnm = v[degree].x - v[0].x;
        double ynm = v[degree].y - v[0].y;
        double xmk = v[0].x;
        double ymk = v[0].y;

        double detInv = -1.0 / ynm;

        return (xnm * ymk - ynm * xmk) * detInv;
    }

    /**
     * Count the number of times a Bezier control polygon crosses the 0-axis. This number is >= the
     * number of roots.
     */
    private static int crossingCount(final KVector[] v, final int degree) {
        int nCrossings = 0;
        int sign = v[0].y < 0 ? -1 : 1;
        int oldSign = sign;
        for (int i = 1; i <= degree; i++) {
            sign = v[i].y < 0 ? -1 : 1;
            if (sign != oldSign) {
                nCrossings++;
            }
            oldSign = sign;
        }
        return nCrossings;
    }

    /**
     * Compute bezier curve.
     * 
     * @param c
     *            control points
     * @param degree
     *            degree of curve
     * @param t
     *            parameter for bezier function
     */
    private static KVector bezier(final KVector[] c, final int degree, final double t,
            final KVector[] left, final KVector[] right) {
        KVector[][] p = new KVector[W_DEGREE + 1][W_DEGREE + 1];

        for (int j = 0; j <= degree; j++) {
            p[0][j] = new KVector(c[j]);
        }

        for (int i = 1; i <= degree; i++) {
            for (int j = 0; j <= degree - i; j++) {
                p[i][j] = new KVector((1.0 - t) * p[i - 1][j].x + t * p[i - 1][j + 1].x, (1.0 - t)
                        * p[i - 1][j].y + t * p[i - 1][j + 1].y);
            }
        }

        if (left != null) {
            for (int j = 0; j <= degree; j++) {
                left[j] = p[j][0];
            }
        }

        if (right != null) {
            for (int j = 0; j <= degree; j++) {
                right[j] = p[degree - j][j];
            }
        }
        return p[degree][0];
    }

    /**
     * Determines the maximum for an arbitrary number of integers.
     * 
     * @param values
     *            integer values
     * @return the maximum of the given values, or {@code MIN_VALUE} if no values are given
     */
    public static int maxi(final int... values) {
        int max = Integer.MIN_VALUE;
        for (int i = 0; i < values.length; i++) {
            if (values[i] > max) {
                max = values[i];
            }
        }
        return max;
    }

    /**
     * Determines the minimum for an arbitrary number of integers.
     * 
     * @param values
     *            integer values
     * @return the minimum of the given values, or {@code MAX_VALUE} if no values are given
     */
    public static int mini(final int... values) {
        int min = Integer.MAX_VALUE;
        for (int i = 0; i < values.length; i++) {
            if (values[i] < min) {
                min = values[i];
            }
        }
        return min;
    }

    /**
     * Determines the average for an arbitrary number of integers.
     * 
     * @param values
     *            integer values
     * @return the average of the given values
     */
    public static int averagei(final int... values) {
        int avg = 0;
        for (int i = 0; i < values.length; i++) {
            avg += values[i];
        }
        return avg / values.length;
    }

    /**
     * Determines the maximum for an arbitrary number of long integers.
     * 
     * @param values
     *            integer values
     * @return the maximum of the given values, or {@code MIN_VALUE} if no values are given
     */
    public static long maxl(final long... values) {
        long max = Integer.MIN_VALUE;
        for (int i = 0; i < values.length; i++) {
            if (values[i] > max) {
                max = values[i];
            }
        }
        return max;
    }

    /**
     * Determines the minimum for an arbitrary number of long integers.
     * 
     * @param values
     *            integer values
     * @return the minimum of the given values, or {@code MAX_VALUE} if no values are given
     */
    public static long minl(final long... values) {
        long min = Integer.MAX_VALUE;
        for (int i = 0; i < values.length; i++) {
            if (values[i] < min) {
                min = values[i];
            }
        }
        return min;
    }

    /**
     * Determines the average for an arbitrary number of long integers.
     * 
     * @param values
     *            integer values
     * @return the average of the given values
     */
    public static long averagel(final long... values) {
        long avg = 0;
        for (int i = 0; i < values.length; i++) {
            avg += values[i];
        }
        return avg / values.length;
    }

    /**
     * Determines the maximum for an arbitrary number of floats.
     * 
     * @param values
     *            float values
     * @return the maximum of the given values, or {@code -MAX_VALUE} if no values are given
     */
    public static float maxf(final float... values) {
        float max = -Float.MAX_VALUE;
        for (int i = 0; i < values.length; i++) {
            if (values[i] > max) {
                max = values[i];
            }
        }
        return max;
    }

    /**
     * Determines the minimum for an arbitrary number of floats.
     * 
     * @param values
     *            float values
     * @return the minimum of the given values, or {@code MAX_VALUE} if no values are given
     */
    public static float minf(final float... values) {
        float min = Float.MAX_VALUE;
        for (int i = 0; i < values.length; i++) {
            if (values[i] < min) {
                min = values[i];
            }
        }
        return min;
    }

    /**
     * Determines the average for an arbitrary number of floats.
     * 
     * @param values
     *            float values
     * @return the average of the given values
     */
    public static float averagef(final float... values) {
        float avg = 0;
        for (int i = 0; i < values.length; i++) {
            avg += values[i];
        }
        return avg / values.length;
    }

    /**
     * Determines the maximum for an arbitrary number of doubles.
     * 
     * @param values
     *            double values
     * @return the maximum of the given values, or {@code -MAX_VALUE} if no values are given
     */
    public static double maxd(final double... values) {
        double max = -Double.MAX_VALUE;
        for (int i = 0; i < values.length; i++) {
            if (values[i] > max) {
                max = values[i];
            }
        }
        return max;
    }

    /**
     * Determines the minimum for an arbitrary number of doubles.
     * 
     * @param values
     *            double values
     * @return the minimum of the given values, or {@code MAX_VALUE} if no values are given
     */
    public static double mind(final double... values) {
        double min = Double.MAX_VALUE;
        for (int i = 0; i < values.length; i++) {
            if (values[i] < min) {
                min = values[i];
            }
        }
        return min;
    }

    /**
     * Determines the average for an arbitrary number of doubles.
     * 
     * @param values
     *            double values
     * @return the average of the given values
     */
    public static double averaged(final double... values) {
        double avg = 0;
        for (int i = 0; i < values.length; i++) {
            avg += values[i];
        }
        return avg / values.length;
    }
    
    /**
     * Limit the given integer to a specific range.
     * 
     * @param x an integer value
     * @param lower the lower limit
     * @param upper the upper limit
     * @return if x is beyond the limits, return the limit, otherwise just return x
     */
    public static int boundi(final int x, final int lower, final int upper) {
        if (x <= lower) {
            return lower;
        } else if (x >= upper) {
            return upper;
        }
        return x;
    }
    
    /**
     * Limit the given long to a specific range.
     * 
     * @param x an long value
     * @param lower the lower limit
     * @param upper the upper limit
     * @return if x is beyond the limits, return the limit, otherwise just return x
     */
    public static long boundl(final long x, final long lower, final long upper) {
        if (x <= lower) {
            return lower;
        } else if (x >= upper) {
            return upper;
        }
        return x;
    }
    
    /**
     * Limit the given float to a specific range.
     * 
     * @param x a float value
     * @param lower the lower limit
     * @param upper the upper limit
     * @return if x is beyond the limits, return the limit, otherwise just return x
     */
    public static float boundf(final float x, final float lower, final float upper) {
        if (x <= lower) {
            return lower;
        } else if (x >= upper) {
            return upper;
        }
        return x;
    }
    
    /**
     * Limit the given double to a specific range.
     * 
     * @param x a double value
     * @param lower the lower limit
     * @param upper the upper limit
     * @return if x is beyond the limits, return the limit, otherwise just return x
     */
    public static double boundd(final double x, final double lower, final double upper) {
        if (x <= lower) {
            return lower;
        } else if (x >= upper) {
            return upper;
        }
        return x;
    }
    
    /**
     * Clip the given vector to a rectangular box of given size.
     * 
     * @param v vector relative to the center of the box
     * @param width width of the rectangular box
     * @param height height of the rectangular box
     * 
     * @return {@code v}.
     */
    public static KVector clipVector(final KVector v, final double width, final double height) {
        double wh = width / 2, hh = height / 2;
        double absx = Math.abs(v.x), absy = Math.abs(v.y);
        double xscale = 1, yscale = 1;
        if (absx > wh) {
            xscale = wh / absx;
        }
        if (absy > hh) {
            yscale = hh / absy;
        }
        v.scale(Math.min(xscale, yscale));
        return v;
    }
    
    /**
     * Returns the signum function of the specified <tt>double</tt> value. The return value
     * is -1 if the specified value is negative; 0 if the specified value is zero; and 1 if the
     * specified value is positive. This is basically {@link Math#signum(double)} with an integer
     * return value.
     *
     * @return the signum function of the specified <tt>double</tt> value.
     */
    public static int signum(final double x) {
        if (x < 0) {
            return -1;
        }
        if (x > 0) {
            return 1;
        }
        return 0;
    }

    /**
     * Checks whether any of the four borders of {@code rect} intersects with any of the straight line segments of the
     * closed {@code path = p_1,...,p_n}. Also checks the closing segment {@code (p_n, p_1)}. If the path contains less
     * than two points, {@code false} is returned. If {@code rect} fully contains the shape, no intersection is assumed.
     * 
     * @param rect
     * @param path
     *            a {@link KVectorChain} describing a closed path of straight line segments.
     * @return {@code true} if {@code rect} intersects {@code path}, {@code false} otherwise.
     */
    public static boolean intersects(final ElkRectangle rect, final KVectorChain path) {
        if (path.size() < 2) {
            return false;
        }

        final Iterator<KVector> pathIt = path.iterator();
        KVector first = pathIt.next();
        KVector p1 = first;
        // check every segment
        while (pathIt.hasNext()) {
            KVector p2 = pathIt.next();
            if (intersects(rect, p1, p2)) {
                return true;
            }
            p1 = p2;
        }
        // check the closing segment
        if (intersects(rect, p1, first)) {
            return true;
        }

        return false;
    }
    
    /**
     * Checks whether the straight line {@code l} defined by {@code (p1, p2)} intersects with at least one of the four
     * segments defining the passed rectangle. If {@code l} is fully contained within {@code rect}'s bounding box, no
     * intersection assumed. If the line is parallel to, and lies on, one of the four borders, no intersection is
     * assumed.
     * 
     * @param rect
     * @param p1
     *            start point of a line
     * @param p2
     *            end point of a line
     * @return {@code true} if {@code (p1, p2)} intersects with {@code rect}, {@code false} otherwise.
     */
    public static boolean intersects(final ElkRectangle rect, final KVector p1, final KVector p2) {
        // simple cases first: fully contained
        if (contains(rect, p1, p2)) {
            return false;
        }
        // leaves the cases
        //  - where one point is inside and the other outside (don't use contains here, as point on border is 'outside')
        //  - where both points are outside
        // for that, check if (p1, p2) intersects with one of rect's borders
        return intersects(rect.getTopLeft(), rect.getTopRight(), p1, p2) 
            || intersects(rect.getTopRight(), rect.getBottomRight(), p1, p2)
            || intersects(rect.getBottomRight(), rect.getBottomLeft(), p1, p2)
            || intersects(rect.getBottomLeft(), rect.getTopLeft(), p1, p2);
    }

    /**
     * Double computations are potentially imprecise, we need a robust way of checking if a value is zero in
     * {@link #intersects(KVector, KVector, KVector, KVector)}, which is done using an epsilon comparison.
     * Note that other algorithms may use this epsilon and depend on its value being no less than 0.00001.
     */
    private static final double DOUBLE_EQ_EPSILON = 0.00001d;
    
    /**
     * Detects intersection of the two passed lines {@code (l11, l12)} and {@code (l21, l22)}.
     * 
     * <p>
     * Implementation based on https://stackoverflow.com/questions/4977491/determining-if-two-line-segments-intersect.
     * See also https://en.wikipedia.org/wiki/Line%E2%80%93line_intersection#Mathematics.
     * </p>
     * 
     * @see ElkMath#intersects2(KVector, KVector, KVector, KVector)
     * @param l11
     *            start point of first line
     * @param l12
     *            end point of first line
     * @param l21
     *            start point of second line
     * @param l22
     *            end point of second line
     * @return {@code true} if the two lines intersect, {@code false} otherwise. In particular, if the start or end
     *         point of one of the lines only "touches" the other line, no intersection is assumed.
     */
    public static boolean intersects(final KVector l11, final KVector l12, final KVector l21, final KVector l22) {
        KVector u0 = l11;
        KVector v0 = l12.clone().sub(l11);
        KVector u1 = l21;
        KVector v1 = l22.clone().sub(l21);
        double x00 = u0.x, y00 = u0.y;
        double x10 = u1.x, y10 = u1.y;
        double x01 = v0.x, y01 = v0.y;
        double x11 = v1.x, y11 = v1.y;
        
        double d = x11 * y01 - x01 * y11;
        if (DoubleMath.fuzzyEquals(0, d, DOUBLE_EQ_EPSILON)) {
            return false;
        }
        double s = (1 / d) * ((x00 - x10) * y01 - (y00 - y10) * x01);
        double t = (1 / d) * -(-(x00 - x10) * y11 + (y00 - y10) * x11);
        // System.out.println("d: " + d + ", s: " + s + ", t: " + t);
        
        // use < instead of <= to not recognize "touching" as intersection
        final boolean intersects =
                // 0 < s
                DoubleMath.fuzzyCompare(0, s, DOUBLE_EQ_EPSILON) < 0
                // s < 1
                && DoubleMath.fuzzyCompare(s, 1, DOUBLE_EQ_EPSILON) < 0
                // 0 < t
                && DoubleMath.fuzzyCompare(0, t, DOUBLE_EQ_EPSILON) < 0
                // t < 1
                && DoubleMath.fuzzyCompare(t, 1, DOUBLE_EQ_EPSILON) < 0;
        return intersects;
    }
    
    /**
     * Returns the intersection point of two line segments.
     * 
     * <p>
     * Approach by Ronald Goldman "Intersection of two lines in three-space"
     * <p>
     * The segments are defined as {@code p+t*r} and {@code q+u*s} with {@code 0 <= t,u <= 1}.
     * If the line segments are collinear and overlapping there is not a single intersection point, so of
     * the possible points the closest to the center of the second segment is returned.
     * 
     * @see ElkMath#intersects(KVector, KVector, KVector, KVector)
     * @param p
     * @param r
     * @param q
     * @param s
     * @return the point of intersection or null if none exists
     */
    public static KVector intersects2(final KVector p, final KVector r, final KVector q, final KVector s) {
        KVector pq = q.clone().sub(p);
        double pqXr = KVector.crossProduct(pq, r);
        double rXs = KVector.crossProduct(r, s);
        double t = KVector.crossProduct(pq, s) / rXs;
        double u = pqXr / rXs;
        if (rXs == 0) {
            if (pqXr == 0) { // segments are collinear: return point closest to center of s
                // CHECKSTYLEOFF MagicNumber
                KVector center = q.clone().add(s.clone().scale(0.5));
                double d1 = p.distance(center);
                double d2 = p.clone().add(r).distance(center);
                double l = s.length() * 0.5;
                // CHECKSTYLEON MagicNumber
                if (d1 < d2 && d1 <= l) {
                    return p.clone();
                }
                if (d2 <= l) {
                    return p.clone().add(r);
                }
                return null;
            } else { // segments are parallel
                return null;
            }
        } else {
            if (t >= 0 && t <= 1 && u >= 0 && u <= 1) { // segments intersect
                return p.clone().add(r.clone().scale(t));
            } else { // segments don't intersect
                return null;
            }
        }
    }
    
    /**
     * Returns the distance between the line segments given by {@code (a1,a2)} and {@code (b1,b2)}.
     * The distance is measured in the direction specified by vector {@code v} starting at segment {@code b}.
     * If the segments are oriented or spaced such that {@code v} cannot reach one segment from the other, 
     * {@code Double.POSITIVE_INFINITY} is returned.
     * 
     * @param a1 start point of first line segment
     * @param a2 end point of first line segment
     * @param b1 start point of second line segment
     * @param b2 end point of second line segment
     * @param v direction
     * @return direction dependent distance between two line segments
     */
    public static double distance(final KVector a1, final KVector a2, final KVector b1, final KVector b2, 
            final KVector v) {
        
        return Math.min(traceRays(a1, a2, b1, b2, v), 
                        traceRays(b1, b2, a1, a2, v.clone().negate()));
    }
    
    /**
     * Traces rays from the endpoints of line segment {@code b=(b1, b2)} in the direction {@code v} and 
     * with the length of {@code v}.
     * The rays' intersections with {@code a} yield the distance between {@code a} and {@code b}.
     * 
     * @param a1 start point of segment a
     * @param a2 end point of segment a
     * @param b1 start point of segment b
     * @param b2 end point of segment b
     * @param v the direction
     * @return the distance between {@code a} and {@code b} in direction {@code v} or infinity
     */
    private static double traceRays(final KVector a1, final KVector a2, final KVector b1, final KVector b2, 
            final KVector v) {
        double result = Double.POSITIVE_INFINITY;
        KVector intersection;
        boolean endpointHit = false;
        
        // check whether (b + v) intersects a to catch an edge case, where one ray exactly hits an endpoint
        // of a but the other ray is too short to reach a.
        intersection = intersects2(a1, a2.clone().sub(a1), b1.clone().add(v), b2.clone().sub(b1));
        boolean edgeCase = intersection != null && !(intersection.equalsFuzzily(a1) || intersection.equalsFuzzily(a2));
        
        // trace ray from point b1 to segment a in direction v
        intersection = intersects2(a1, a2.clone().sub(a1), b1, v);
        if (intersection != null) {
            if (intersection.equalsFuzzily(a1) == intersection.equalsFuzzily(a2) || edgeCase) {
                // update the distance result
                result = Math.min(result, intersection.sub(b1).length());
            } else {
                // ignore intersection if the ray exactly hits an endpoint of b
                // but track it to allow it only for one ray per segment
                endpointHit = true;
            }
        }
        
        // trace ray from point b2 to segment a in direction v
        intersection = intersects2(a1, a2.clone().sub(a1), b2, v);
        if (intersection != null) {
            if (endpointHit || intersection.equalsFuzzily(a1) == intersection.equalsFuzzily(a2) || edgeCase) {
                // update the distance result
                result = Math.min(result, intersection.sub(b2).length());
            }
        }
        
        return result;
    }
    
    /**
     * Checks whether every straight line segment of the path {@code p_1,...,p_n} 
     * is fully contained within {@code rect}. Also checks the closing segment {@code (p_n, p_1)}.
     * If the path contains less than two points, {@code false} is returned.
     * 
     * @param rect
     * @param path a {@link KVectorChain} describing a closed path of straight line segments. 
     * @return {@code true} if {@code rect} fully contains {@code path}, {@code false} otherwise.
     */
    public static boolean contains(final ElkRectangle rect, final KVectorChain path) {
        if (path.size() < 2) {
            return false;
        }
        
        final Iterator<KVector> pathIt = path.iterator();
        KVector first = pathIt.next();
        KVector p1 = first; 
        // check every segment
        while (pathIt.hasNext()) {
            KVector p2 = pathIt.next();
            if (!contains(rect, p1, p2)) {
                return false;
            }
            p1 = p2;
        }
        // check the closing segment
        if (!contains(rect, p1, first)) {
            return false;
        }
        
        return true;
    }

    /**
     * Check whether {@code rect} fully contains the straight line {@code l} defined by {@code (p1, p2)}. That is, 
     * it is checked if both {@code p1} and {@code p2} lie within the bounding box of {@code rect}. If one of the 
     * two points lies on the border of {@code rect}, the line is considered to be <em>not</em> contained. 
     * 
     * @param rect
     * @param p1 start point of a line 
     * @param p2 end point of a line
     * @return {@code true} if {@code rect} fully contains {@code (p1, p2)}, {@code false} otherwise.
     */
    public static boolean contains(final ElkRectangle rect, final KVector p1, final KVector p2) {
        return contains(rect, p1) && contains(rect, p2);        
    }
   
    /**
     * Check whether {@code rect} contains the point {@code p}. That is, it is checked if {@code p} lies within the
     * bounding box of {@code rect}. If the point lies on the border of {@code rect}, the line is considered to be
     * <em>not</em> contained.
     * 
     * @param rect
     * @param p
     * @return {@code true} if {@code rect} contains {@code p1}, {@code false} otherwise.
     */
    public static boolean contains(final ElkRectangle rect, final KVector p) {
        double minX = rect.x;
        double maxX = rect.x + rect.width;
        double minY = rect.y;
        double maxY = rect.y + rect.height;
        
        return (p.x > minX && p.x < maxX) && (p.y > minY && p.y < maxY);     
    }
    
    /**
     * Calculates shortest distance between two axis aligned rectangles.
     * If they intersect the returned distance becomes negative.
     * There are three cases to consider:
     * 1.
     *              +----+
     *    +----+....| r2 |
     *    | r1 |    +----+
     *    +----+
     *    
     * 2.
     *              +----+
     *              | r2 |
     *             .+----+
     *           .
     *    +----+
     *    | r1 |   
     *    +----+
     *    
     * 3. 
     *
     *        +----+
     *    +---|+r2 |
     *    | r1+|---+   
     *    +----+
     *    
     * In the first case the result is the axis aligned distance between the closest sides
     * and in the other cases it's the euclidean distance between the closest corners.
     * @param r1
     * @param r2
     * @return shortest distance between r1 and r2
     */
    public static double shortestDistance(final ElkRectangle r1, final ElkRectangle r2) {
        double rightDist = r2.x - (r1.x + r1.width);
        double leftDist = r1.x - (r2.x + r2.width);
        double topDist = r1.y - (r2.y + r2.height);
        double bottomDist = r2.y - (r1.y + r1.height);
        double horzDist = Math.max(leftDist, rightDist);
        double vertDist = Math.max(topDist, bottomDist);
        if (DoubleMath.fuzzyCompare(horzDist, 0, DOUBLE_EQ_EPSILON) >= 0 
                ^ DoubleMath.fuzzyCompare(vertDist, 0, DOUBLE_EQ_EPSILON) >= 0) { // case 1
            return Math.max(vertDist, horzDist);
        }
        if (DoubleMath.fuzzyCompare(horzDist, 0, DOUBLE_EQ_EPSILON) > 0) { // case 2
            return Math.sqrt(vertDist * vertDist + horzDist * horzDist);
        }
        // case 3
        return -Math.sqrt(vertDist * vertDist + horzDist * horzDist);
    }

}
