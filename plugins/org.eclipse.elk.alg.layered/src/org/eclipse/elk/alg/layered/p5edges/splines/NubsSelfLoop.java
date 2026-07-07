/*******************************************************************************
 * Copyright (c) 2014, 2015 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.layered.p5edges.splines;

import java.util.List;

import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.core.math.KVector;

/**
 * An extension of {@link NubSpline} for self-loops. Extends the parent class by methods for creation of self-loop
 * splines. Use the static create methods to obtain an instance to this class.
 */
public final class NubsSelfLoop extends NubSpline {
    /**
     * Fraction of the second/fourth control-point compared to the desired high. A bigger value will result in a
     * "higher" second/fourth control-point and a lower middle control-point.
     */
    private static final double FRACTION = 1.3;
    /** To avoid magic number warning. */
    private static final double HALF = 0.5;

    /** The position of the first self-loop label. */
    private KVector firstLabelPosition = new KVector();
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Constructors

    /**
     * A copy constructor.
     * 
     * @param nubs
     *            The NubsSelfLoop to copy.
     */
    private NubsSelfLoop(final NubsSelfLoop nubs) {
        super(nubs);
        firstLabelPosition = nubs.firstLabelPosition;
    }

    /**
     * Creates a new uniform and clamped NubsSelfLoop with the specified control points.
     * 
     * @param dimension
     *            The dimension if this NubsSelfLoop.
     * @param kVectors
     *            The control points of this NubsSelfLoop.
     */
    private NubsSelfLoop(final int dim, final KVector... kVectors) {
        super(true, dim, kVectors);
    }

    /**
     * Creates a new uniform and clamped NubsSelfLoop with the specified control points. Dimension must be > 0.
     * 
     * @param dimension
     *            The dimension if this NubsSelfLoop.
     * @param kVectors
     *            The control points of this NubsSelfLoop.
     */
    public NubsSelfLoop(final int dimension, final List<KVector> kVectors) {
        super(true, dimension, kVectors);
    }

    /**
     * Creates a new uniform NubsSelfLoop with the specified control points. Dimension must be > 0.
     * 
     * @param clamped
     *            True, if the NubsSelfLoop shall be clamped.
     * @param dimension
     *            The dimension if this NubsSelfLoop.
     * @param kVectors
     *            The control points of this NubsSelfLoop.
     */
    public NubsSelfLoop(final boolean clamped, final int dimension, final List<KVector> kVectors) {
        super(clamped, dimension, kVectors);
    }
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Accessors
    
    /**
     * Returns the first label's position. May be modified.
     */
    public KVector getFirstLabelPosition() {
        return firstLabelPosition;
    }

    /**
     * Sets the label position.
     * 
     * @param position
     *            The new position for the first label.
     */
    public void setFirstLabelPosition(final KVector position) {
        firstLabelPosition = position;
    }
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Self Loop Builders

    /**
     * Generates a new NUBS representing a self loop at one of the straight sides of a node.
     * 
     * @param source
     *            Source port this spline shall start at.
     * @param target
     *            Target port this spline shall end at.
     * @param length
     *            The estimated distance the spline shall go away from the side of the node.
     * @return A new NUBS representing the loop spline.
     */
    public static NubsSelfLoop createSideSelfLoop(final LPort source, final LPort target, final double length) {

        double direction = SplinesMath.portSideToDirection(source.getSide());
        KVector sourcePos = source.getPosition().clone().add(source.getAnchor());
        KVector targetPos = target.getPosition().clone().add(target.getAnchor());

        // calculate the three control points (next to source and target point)
        KVector firstCP = new KVector(sourcePos).add(new KVector(direction).scale(length));
        KVector thirdCP = new KVector(targetPos).add(new KVector(direction).scale(length));
        KVector mid = new KVector(firstCP).sub(thirdCP).scale(1.0 / 2.0);
        KVector secondCP = new KVector(thirdCP).add(mid).add(new KVector(direction).scale(mid.length()));

        // calculate the spline
        NubsSelfLoop nubs = new NubsSelfLoop(DIM, sourcePos, firstCP, secondCP, thirdCP, targetPos);

        return nubs;
    }

    /**
     * Generates a new NUBS representing a self loop at the corner of a node.
     * 
     * @param source
     *            The source port this spline shall start at.
     * @param target
     *            The target port this spline shall end at.
     * @param sourceHeight
     *            The estimated distance the spline shall go away from the source side.
     * @param targetHeight
     *            The estimated distance the spline shall go away from the target side.
     * @return A new NubSpline representing the loop spline.
     */
    public static NubsSelfLoop createCornerSelfLoop(final LPort source, final LPort target, final double sourceHeight,
            final double targetHeight) {
        
        // Radian direction of source and target.
        double sourceDir = SplinesMath.portSideToDirection(source.getSide());
        double targetDir = SplinesMath.portSideToDirection(target.getSide());
        
        // Position of source and target (port).
        KVector sourcePos = source.getPosition().clone().add(source.getAnchor());
        KVector targetPos = target.getPosition().clone().add(target.getAnchor());

        // first and last CP are simply located in a specified distance (sourceHeight/targetHeight)
        // and a specified direction (sourceDirection/targetDirection) from source/target.
        KVector firstCP = sourcePos.clone().add(new KVector(sourceDir).scale(FRACTION * sourceHeight));
        KVector thirdCP = targetPos.clone().add(new KVector(targetDir).scale(FRACTION * targetHeight));

        // second CP is located at a defined scale of the vector from corner to cpCorner
        double cornerX = 0.0;
        double cornerY = 0.0;

        switch (source.getSide()) {
        case WEST:
            cornerX = 2 * (sourcePos.x - sourceHeight) - HALF * (firstCP.x + thirdCP.x);
            break;
        case EAST:
            cornerX = 2 * (sourcePos.x + sourceHeight) - HALF * (firstCP.x + thirdCP.x);
            break;
        case NORTH:
            cornerY = 2 * (sourcePos.y - sourceHeight) - HALF * (firstCP.y + thirdCP.y);
            break;
        case SOUTH:
            cornerY = 2 * (sourcePos.y + sourceHeight) - HALF * (firstCP.y + thirdCP.y);
            break;
        default:
            break;
        }
        
        switch (target.getSide()) {
        case WEST:
            cornerX = 2 * (targetPos.x - targetHeight) - HALF * (thirdCP.x + firstCP.x);
            break;
        case EAST:
            cornerX = 2 * (targetPos.x + targetHeight) - HALF * (thirdCP.x + firstCP.x);
            break;
        case NORTH:
            cornerY = 2 * (targetPos.y - targetHeight) - HALF * (thirdCP.y + firstCP.y);
            break;
        case SOUTH:
            cornerY = 2 * (targetPos.y + targetHeight) - HALF * (thirdCP.y + firstCP.y);
            break;
        default:
            break;
        }

        KVector secondCP = new KVector(cornerX, cornerY);

        // Calculate the spline
        return new NubsSelfLoop(DIM, sourcePos, firstCP, secondCP, thirdCP, targetPos);
    }

    /**
     * Generates a new NUBS representing a self loop going from one node side to the opposite side. Basic bend points
     * on the way are already computed.
     * 
     * @param source
     *            the source port.
     * @param target
     *            the target port.
     * @param sourceBendPoint
     *            first bend point from source to target.
     * @param cornerBendPoints
     *            bend points between the first and the last.
     * @param targetBendPoint
     *            last bend point from source to target.
     */
    public static NubsSelfLoop createAcrossSelfLoop(final LPort source, final LPort target,
            final KVector sourceBendPoint, final List<KVector> cornerBendPoints, final KVector targetBendPoint) {

        KVector sourcePos = source.getPosition().clone().add(source.getAnchor());
        KVector targetPos = target.getPosition().clone().add(target.getAnchor());

        assert cornerBendPoints.size() == 2;
        KVector thirdCP = cornerBendPoints.get(0);
        KVector fourthCP = cornerBendPoints.get(1);

        // calculate middle (between second and third) control point
        KVector middleCP = fourthCP.clone().add(thirdCP).scale(HALF);

        // calculate the spline
        return new NubsSelfLoop(DIM, sourcePos, sourceBendPoint, thirdCP, middleCP, fourthCP, targetBendPoint,
                targetPos);
    }

    /**
     * Generates a new NUBS representing a self loop that spans three sides. Basic bend points on the way are already
     * computed.
     * 
     * @param source
     *            the source port.
     * @param target
     *            the target port.
     * @param sourceBendPoint
     *            first bend point from source to target.
     * @param cornerBendPoints
     *            bend points between the first and the last.
     * @param targetBendPoint
     *            last bend point from source to target.
     */
    public static NubsSelfLoop createThreeSideSelfLoop(final LPort source, final LPort target,
            final KVector sourceBendPoint, final List<KVector> cornerBendPoints, final KVector targetBendPoint) {
        
        KVector sourcePos = source.getPosition().clone().add(source.getAnchor());
        KVector targetPos = target.getPosition().clone().add(target.getAnchor());
        
        // SUPPRESS CHECKSTYLE NEXT 1 MagicNumber
        assert cornerBendPoints.size() == 3;
        KVector thirdCP = cornerBendPoints.get(0);
        KVector fourthCP = cornerBendPoints.get(1);
        KVector fifthCP = cornerBendPoints.get(2);

        // calculate the spline
        return new NubsSelfLoop(DIM, sourcePos, sourceBendPoint, thirdCP, fourthCP, fifthCP, targetBendPoint,
                targetPos);
    }

    /**
     * Generates a new NUBS representing a self loop that spans three sides. Basic bend points on the way are already
     * computed.
     * 
     * @param source
     *            the source port.
     * @param target
     *            the target port.
     * @param sourceBendPoint
     *            first bend point from source to target.
     * @param cornerBendPoints
     *            bend points between the first and the last.
     * @param targetBendPoint
     *            last bend point from source to target.
     */
    public static NubsSelfLoop createFourSideSelfLoop(final LPort source, final LPort target,
            final KVector sourceBendPoint, final List<KVector> cornerBendPoints, final KVector targetBendPoint) {
        
        KVector sourcePos = source.getPosition().clone().add(source.getAnchor());
        KVector targetPos = target.getPosition().clone().add(target.getAnchor());

        // SUPPRESS CHECKSTYLE NEXT 5 MagicNumber
        assert cornerBendPoints.size() == 4;
        KVector thirdCP = cornerBendPoints.get(0);
        KVector fourthCP = cornerBendPoints.get(1);
        KVector fifthCP = cornerBendPoints.get(2);
        KVector sixthCP = cornerBendPoints.get(3);

        // calculate the spline
        return new NubsSelfLoop(DIM, sourcePos, sourceBendPoint, thirdCP, fourthCP, fifthCP, sixthCP, targetBendPoint,
                targetPos);
    }

}
