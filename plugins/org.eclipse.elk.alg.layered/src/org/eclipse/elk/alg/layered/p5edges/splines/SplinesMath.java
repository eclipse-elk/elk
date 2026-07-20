/*******************************************************************************
 * Copyright (c) 2014, 2015 Kiel University and others.
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
package org.eclipse.elk.alg.layered.p5edges.splines;

import java.text.DecimalFormat;

import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.options.PortSide;


/**
 * Mathematics utility class used in the splines routing.
 * 
 * @author tit
 * 
 */
public final class SplinesMath {
    
    /** Differences below this value are treated as zero. */
    private static final double EPSILON = 0.00000001;
    /** The representation for a null value in the output. */
    private static final String NULL_STRING = "(null)";
    // some class statics to increase performance //
    /** 1/2 * Pi. */
    public static final double HALF_PI = Math.PI / 2;
    /** 1/4 * Pi. */
    public static final double QUATER_PI = HALF_PI / 2;
    /** 3/2 * Pi. */
    public static final double THREE_HALF_PI = HALF_PI + HALF_PI + HALF_PI;
    /** 2 * Pi. */
    public static final double TWO_PI = 2.0 * Math.PI;
    /** Three. */
    public static final double THREE = 3;

    /**
     * Private default constructor to prevent instantiating this class.
     */
    private SplinesMath() {
        throw new AssertionError("Instantiating utility class SplinesMath.");
    }

    /**
     * Calculates the intersection point of two straight lines. If lines are parallel, null is
     * returned.
     * 
     * @param pt1 Point on first straight line.
     * @param pt2 Point on second straight line.
     * @param dirPT2 Direction of edge going from pt0 to intersection point
     * @param dirPT1 Direction of edge going from pt1 to intersection point
     * @return A new created KVector representing the intersection or null.
     */
    public static KVector intersect(final KVector pt1, final KVector pt2, final double dirPT1,
            final double dirPT2) {
        return intersect(pt1, new KVector(dirPT1).add(pt1), pt2, new KVector(dirPT2).add(pt2));
    }

    /**
     * Calculates the intersection point of two straight lines. If lines are parallel, null is
     * returned.
     * 
     * @param pt1 First point on first straight line.
     * @param pt2 Second point on first straight line.
     * @param pt3 First point on second straight line.
     * @param pt4 Second point on second straight line.
     * @return A new created KVector representing the intersection or null if lines are parallel.
     */
    public static KVector intersect(final KVector pt1, final KVector pt2, final KVector pt3,
            final KVector pt4) {
        final double x1 = pt1.x;
        final double y1 = pt1.y;
        final double x2 = pt2.x;
        final double y2 = pt2.y;
        final double x3 = pt3.x;
        final double y3 = pt3.y;
        final double x4 = pt4.x;
        final double y4 = pt4.y;

        final double divisor = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
        if (Math.abs(divisor) < EPSILON) {
            // The two lines are parallel or even coincident 
            return null;
        }

        final double first = x1 * y2 - y1 * x2;
        final double second = x3 * y4 - y3 * x4;
        
        final double newX = (first * (x3 - x4) - second * (x1 - x2)) / divisor;
        final double newY = (first * (y3 - y4) - second * (y1 - y2)) / divisor;
        
        return new KVector(newX, newY);
    }

    /**
     * Calculates the inner angle of two KVectors.
     * 
     * @param vec1 First vector.
     * @param vec2 Second Vector.
     * @return The inner angle of vec0 and vec1 in radians. Value is always between 0 and PI (180°).
     */
    public static double innerAngle(final KVector vec1, final KVector vec2) {
        return Math.acos((vec1.x * vec2.x + vec1.y * vec2.y) / (vec1.length() * vec2.length()));
    }
    
    /**
     * Calculates the inner angle of two directions in radians.
     * 
     * @param dir1 First direction.
     * @param dir2 Second direction.
     * @return The inner angle of dir0 and dir1 in radians. Value is always between 0 and PI (180°).
     */
    public static double innerAngle(final double dir1, final double dir2) {
        final double retval = Math.abs(dir1 - dir2);
        return retval % Math.PI;
    }

    
    // elkjs-exclude-start
    // DecimalFormat not available in gwt
    /** Option for debug formating. */
    private static final DecimalFormat DEC_FORMAT = new DecimalFormat("#0.0");
    
    /**
     * Constructs a string showing the coordinates of the given KVector.
     * 
     * @param vector
     *            The KVector to display.
     * @return A string representation of the KVector.
     */
    public static String convertKVectorToString(final KVector vector) {
        if (vector == null) {
            return NULL_STRING;
        }
        return "(" + DEC_FORMAT.format(vector.x) + "," + DEC_FORMAT.format(vector.y) + ")";
    }

    /**
     * Converts all given KVectors to a readable string with rounded results.
     * @param list The list of KVectors to create the string for.
     * @return The readable string.
     */
    public static String convertKVectorToString(final KVector... list) {
        if (list == null || list.length == 0) {
            return NULL_STRING;
        }
        final StringBuilder retVal = new StringBuilder();
        
        for (final KVector vector : list) {
            retVal.append(convertKVectorToString(vector)).append(", ");
        }

        return retVal.substring(0, retVal.length() - 2);
    }
    
    /**
     * Converts all given KVectors to a readable string with rounded results.
     * @param list The list of KVectors to create the string for.
     * @return The readable string.
     */
    public static String convertKVectorToString(final Iterable<KVector> list) {
        if (list == null || !list.iterator().hasNext()) {
            return NULL_STRING;
        }
        final StringBuilder retVal = new StringBuilder();
        
        for (final KVector vector : list) {
            retVal.append(convertKVectorToString(vector)).append(", ");
        }

        return retVal.substring(0, retVal.length() - 2);
    }
    // elkjs-exclude-end

    /**
     * Calculates a vector of given direction up to the coordinate, where it's perpendicular through 
     * the given point crosses the vector. The vector's length is returned.
     * 
     * @param direction
     *            Direction of the vector.
     * @param point
     *            Point that the perpendicular shall pass through.
     * @return The vector's length.
     */
    public static double lengthToOrthogonal(final double direction, final KVector point) {
        double angle = innerAngle(direction, point.toRadians());
        int factor = 1;
        
        if (angle > HALF_PI) {
            // point lays in the opposite direction
            factor = -1;
            angle = angle - QUATER_PI;
        }
        return factor * Math.cos(angle) * point.length();
    }

    /**
     * Converts a {@link PortSide} to the direction from a node's center to the given side in radian.
     * 
     * @param side The portSide to convert. 
     * @return The direction in radian.
     */
    public static double portSideToDirection(final PortSide side) {
        switch (side) {
        case NORTH:
            return SplinesMath.THREE_HALF_PI;
        case EAST:
            return 0.0;
        case SOUTH:
            return SplinesMath.HALF_PI;
        case WEST:
            return Math.PI;
        default:
            return 0.0;
        }
    }
    
    /**
     * Calculates the perpendicular to the given nodeSide through the AbsolutAnchor of the given port.
     * The distance of the AbsolutAnchor to the nodeSide on this perpendicular is returned.
     *
     * @param port The port to calculate the distance for.
     * @param side The side of the node to calculate the distance for.
     * @return The distance.
     */
    public static double distPortToNodeEdge(final LPort port, final PortSide side) {
        final KVector portPos = port.getPosition().clone().add(port.getAnchor());
        final KVector nodeSize = port.getNode().getSize();
        
        switch (side) {
        case NORTH:
            return -portPos.y;
        case EAST:
            return -portPos.x + nodeSize.x;
        case SOUTH:
            return -portPos.y + nodeSize.y;
        case WEST:
            return -portPos.x;
        default:
            return 0.0;
        }
    }
    
    /**
     * Checks if the given value lays between (or on) the two boundaries. No matter which of them is 
     * larger.
     * @param value The value to check.
     * @param boundary0 The first boundary. 
     * @param boundary1 The second boundary.
     * @return {@code true} if one boundary is {@code (<= value)} 
     * and the other one is {@code (>= value)}.
     */
    public static boolean isBetween(final int value, final int boundary0, final int boundary1) {
        return value < boundary0 
                ? boundary1 <= value 
                : value <= boundary1 || value == boundary0;
    }

    /**
     * Checks if the given value lays between (or on) the two boundaries. No matter which of them is 
     * larger.
     * @param value The value to check.
     * @param boundary0 The first boundary. 
     * @param boundary1 The second boundary.
     * @return {@code true} if one boundary is {@code (<= value)} 
     * and the other one is {@code (>= value)}.
     */
    public static boolean isBetween(final double value, final double boundary0, final double boundary1) {
        if (Math.abs(boundary0 - value) < EPSILON || Math.abs(boundary1 - value) < EPSILON) {
            return true;
        }
        return (boundary0 - value) > EPSILON 
                ? (value - boundary1) > EPSILON 
                : (boundary1 - value) > EPSILON;
    }
    
    /**
     * Returns the margin on the given side of a node.
     * @param node The node those margin to return.
     * @param side The port-side those margin to return.
     * @return The margin.
     */
    public static double getMarginOnPortSide(final LNode node, final PortSide side) {
        switch (side) {
        case NORTH:
            return node.getMargin().top;
        case EAST:
            return node.getMargin().right;
        case SOUTH:
            return node.getMargin().bottom;
        case WEST:
            return node.getMargin().left;
        default:
            return 0.0;
        }
    }
 }
