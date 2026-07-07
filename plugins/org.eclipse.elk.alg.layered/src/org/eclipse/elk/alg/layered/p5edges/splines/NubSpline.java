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

import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.ListIterator;

import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.math.KVectorChain;

import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;

/**
 * Represents a Non Uniform B-Spline. This spline is not rational, thus there is no weight for the
 * control points.
 * The control points of the represented spline are stored as PolarCPs. This class also holds their
 * name in the polar form, that is described by Sederberg. 
 * For details of this notation, please refer Senderberg: "An introduction to B-Spline Curves" (2005).
 * 
 * A control point those polar coordinates are all the same is laying on the spline.
 * 
 * The knot vector does not store the irrelevant first and last vector traditionally added to it.
 * These two vectors are meaningless, thus don't need to be stored. Remember this if you want to 
 * directly work with the knot vector.
 * 
 * @author tit
 * 
 */
public class NubSpline {
    /** Default dimension of the spline. */
    protected static final int DIM = 3;
    /** Doubles with a difference less than this value will be assumed to be equal. */
    private static final double EPSILON = 0.000001;

    /** The knotVector of this NubSpline. (The traditionally added 0 at the beginning and 1 at the end
     *  of the knot vector is not added, as these elements are not relevant for the calculation) */
    private List<Double> knotVector = Lists.newArrayList();
    /** The control points of this NubSpline. */
    private List<PolarCP> controlPoints = Lists.newArrayList();

    /** The dimension of this NubSpline. All contained and constructible vectors have the same dim. */
    private int dimNUBS;
    /** Specifies of this NubSpline is uniform. To be exact this only depends on the knotVecor. 
     * A NubSpline is uniform if the difference of two adjacent knots in the knot vector is either 0 or 
     * equal to a constant valid for the whole spline.*/
    private boolean isUniform;
    /** Specifies of this NubSpline is clamped. To be exact this only depends on the knotVecor. 
     * A NubSpline is clamped, if the first and last vector of it's knot vector has the multiplicity
     * equal to it's dimension. (Remember the neglected first and last vector. Otherwise this would
     * be dimension +1.) */
    private boolean isClamped;
    /** The outer margins of this spline. */
    private Rectangle outerBox;
    /** Indicates if all inner knots of the knotVector have the multiplicity of dimNUBS. */
    private boolean isBezier;
    /** The value of parameter t (time), where the resulting spline would start. */
    private double minKnot;
    /** The value of parameter t (time), where the resulting spline would end. */
    private double maxKnot;

    // ########################################################################################
    // Constructors

    /**
     * A copy constructor.
     * 
     * @param nubSpline The NubSpline to copy.
     */
    public NubSpline(final NubSpline nubSpline) {
        dimNUBS = nubSpline.dimNUBS;
        isUniform = nubSpline.isUniform;
        isClamped = nubSpline.isClamped;
        outerBox = nubSpline.outerBox;
        isBezier = nubSpline.isBezier;
        knotVector = Lists.newLinkedList(nubSpline.knotVector);
        minKnot = nubSpline.minKnot;
        maxKnot = nubSpline.maxKnot;
        controlPoints = Lists.newLinkedList(nubSpline.controlPoints);
    }

    /**
     * Creates a new uniform NubSpline with the specified control points. The dimension must be > 0.
     * 
     * @param clamped {@code true}, if the NubSpline shall be clamped.
     * @param dimension The dimension of this NubSpline.
     * @param kVectors The control points of this NubSpline.
     */
    public NubSpline(final boolean clamped, final int dimension, final List<KVector> kVectors) {
        if (dimension < 1) {
            throw new IllegalArgumentException("The dimension must be at least 1!");
        }
        
        // fill the list of control-points to be at least equal to dimension + 1 
        for (int i = kVectors.size() - 1; i < dimension; i++) {
            kVectors.add(0, kVectors.get(0));
        }

        if (kVectors.size() < (dimension + 1)) {
            throw new IllegalArgumentException(
                    "At (least dimension + 1) control points are necessary!");
        } else {
            dimNUBS = dimension;
            isClamped = clamped;
            isUniform = true;
            isBezier = false;

            // create the knot vector
            createUniformKnotVector(clamped, kVectors.size() + dimNUBS - 1);

            final List<Double> polarCoordinate = Lists.newArrayList();
            final Iterator<Double> knotIter = knotVector.iterator();

            // the first (dimNUBS - 1) elements of the knotVector for the "sliding window" that
            // determines the polarCoordinates of the PolarCP.
            for (int i = 0; i < (dimNUBS - 1); i++) {
                polarCoordinate.add(knotIter.next());
            }

            // Create the PolarCPs
            // A sliding window over the knot vector determines the polar coordinates of the
            // PolarCPs.
            for (final KVector kVector : kVectors) {
                polarCoordinate.add(knotIter.next());
                controlPoints.add(new PolarCP(kVector, polarCoordinate));
                polarCoordinate.remove(0);
            }
        }
    }

    /**
     * Creates a new uniform and clamped NubSpline with the specified control points. Dimension must be > 0.
     * 
     * @param clamped {@code true}, if the NubSpline shall be clamped.
     * @param dimension The dimension if this NubSpline.
     * @param kVectors The control points of this NubSpline.
     */
    public NubSpline(final boolean clamped, final int dimension, final KVector... kVectors) {
        this(clamped, dimension, Lists.newArrayList(kVectors));
    }

    /**
     * Creates a new NubSpline with the given, precalculated values.
     * There is no error checking! So if the spline shall be clamped, uniform or bezier, you have to 
     * ensure the correctness of the knot vector on your own. 
     * 
     * @param clamped {@code true}, if the NubSpline is clamped.
     * @param uniform {@code true}, if the NubSpline is uniform.
     * @param bezier {@code true}, if the NubSpline is in bezier form.
     * @param dim The dimension of the NubSpline.
     * @param knotVec The knot vector if the NubSpline.
     * @param polarVectors The control points of the NubSpline already in polar form.
     */
    private NubSpline(final boolean clamped, final boolean uniform, final boolean bezier, 
            final int dim, final List<Double> knotVec, final List<PolarCP> polarVectors) {
        isClamped = clamped;
        isUniform = uniform;
        isBezier = bezier;
        dimNUBS = dim;
        knotVector = knotVec;
        controlPoints = polarVectors;
        minKnot = knotVec.iterator().next();
        maxKnot = Iterables.getLast(knotVec);
    }

    // ########################################################################################
    // Generators

    /**
     * Creates a new NubSpline, that is the derivation of the given NubSpline.
     * 
     * @param nubSpline The NubSpline to create the derivation from.
     * @return The derivation.
     */
    public static NubSpline generateDerivedNUBS(final NubSpline nubSpline) {
        final boolean newClamped = nubSpline.isClamped;
        final boolean newUniform = nubSpline.isUniform;
        final boolean newBezier = nubSpline.isBezier;
        final int oldDim = nubSpline.dimNUBS;
        final int newDim = oldDim - 1;
        final List<Double> oldKnotVector = nubSpline.knotVector;
        final List<Double> newKnotVector = Lists.newLinkedList(
                nubSpline.knotVector.subList(1, nubSpline.knotVector.size() - 1));
        final List<KVector> newControlPoints = Lists.newArrayList();

        // Calculate the new control points.
        for (int i = 0; i < nubSpline.controlPoints.size() - 1; i++) {
            final KVector newCP =
                    nubSpline.controlPoints.get(i + 1).getCp().clone()
                            .sub(nubSpline.controlPoints.get(i).getCp())
                            .scale(oldDim / (oldKnotVector.get(i + oldDim) - oldKnotVector.get(i)));
            newControlPoints.add(newCP);
        }

        // Create the PolarCPs
        final List<Double> polarCoordinate = Lists.newArrayList();
        final Iterator<Double> knotIter = newKnotVector.iterator();
        final List<PolarCP> newPolarVectors = Lists.newArrayList();

        // the first (dimNUBS - 1) elements of the knotVector for the "sliding window" that
        // determines the polarCoordinates of the PolarCP.
        for (int i = 0; i < (newDim - 1); i++) {
            polarCoordinate.add(knotIter.next());
        }

        // A sliding window over the knot vector determines the polar coordinates of the PolarCPs.
        for (final KVector kVector : newControlPoints) {
            polarCoordinate.add(knotIter.next());
            newPolarVectors.add(new PolarCP(kVector, polarCoordinate));
            polarCoordinate.remove(0);
        }

        // Create the new NubSpline.
        return new NubSpline(newClamped, newUniform, newBezier, newDim, newKnotVector, newPolarVectors);
    }
    
    /**
     * Creates a new NubSpline, that is a inverted copy of the given one.
     * The control-points will be inverted.
     * The distribution of the vectors on the knot-vector will be inverted, but the knot-vector will
     * be sorted ascending, again.
     *  
     * @param nubSpline The NubSpline to create a inverted copy from.
     * @return The inverted NubSpline.
     */
    public static NubSpline generateInvertedNUBS(final NubSpline nubSpline) {
        final List<Double> newKnotVector = Lists.newArrayList();
        final double maxVector =  nubSpline.knotVector.get(nubSpline.knotVector.size() - 1); 
        for (final Double vector : nubSpline.knotVector) {
            newKnotVector.add(0, maxVector - vector);
        }
        final List<KVector> newControlPoints = KVectorChain.reverse(nubSpline.getControlPoints());

        // Create the PolarCPs
        final List<Double> polarCoordinate = Lists.newArrayList();
        final Iterator<Double> knotIter = newKnotVector.iterator();
        final List<PolarCP> newPolarVectors = Lists.newArrayList();

        // the first (dimNUBS - 1) elements of the knotVector for the "sliding window" that
        // determines the polarCoordinates of the PolarCP.
        for (int i = 0; i < (nubSpline.dimNUBS - 1); i++) {
            polarCoordinate.add(knotIter.next());
        }

        // A sliding window over the knot vector determines the polar coordinates of the PolarCPs.
        for (final KVector kVector : newControlPoints) {
            polarCoordinate.add(knotIter.next());
            newPolarVectors.add(new PolarCP(kVector, polarCoordinate));
            polarCoordinate.remove(0);
        }

        // Create the new NubSpline.
        return new NubSpline(
                nubSpline.isClamped, 
                nubSpline.isUniform, 
                nubSpline.isBezier, 
                nubSpline.dimNUBS, 
                newKnotVector, 
                newPolarVectors);
    }
    
    // ########################################################################################
    // Getter & setter

    /**
     * Returns the dimension of this NubSpline.
     * 
     * @return The dimension.
     */
    public int getDim() {
        return dimNUBS;
    }

    /**
     * Returns and, if necessary, calculates the outerBox of this spline. The outerBox lays around 
     * all parts of the spline.
     * @return A Rectangle representing a square around this NubSpline. The Rectangle is only an
     *         approximation.
     */
    public Rectangle getOuterBox() {
        if (outerBox == null) {
            calculateOuterBox();
        }
        return outerBox;
    }
    
    /**
     * Sets the outer box of this NubSpline. Only use if you know what you are doing!
     * 
     * @param outerRectange The outerBox of this spline.
     */
    protected void setOuterBox(final Rectangle outerRectange) {
        outerBox = outerRectange;
    }

    /**
     * Returns the control points. The control points may be modified,
     * but the list itself is only a container.
     * 
     * @return The control points.
     */
    public KVectorChain getControlPoints() {
        final KVectorChain retVal = new KVectorChain();

        for (final PolarCP polarCP : controlPoints) {
            retVal.add(polarCP.getCp());
        }
        return retVal;
    }

    /**
     * Returns the i'th control-point. May be modified.
     * 
     * @param i Position of the control-point.
     * @return The control-point.
     */
    public KVector getControlPoint(final int i) {
        return controlPoints.get(i).getCp();
    }

    /**
     * Returns the number of controlPoints in this NubSpline.
     * 
     * @return The number of control points.
     */
    public int getControlPointsSize() {
        return controlPoints.size();
    }

    /**
     * Returns a copy of the knot vector of this NubSpline.
     * 
     * @return A copy of the knotVector.
     */
    public List<Double> getKnotVector() {
        return Lists.newLinkedList(knotVector);
    }


    /**
     * Searches the given knot in the knot-vector.
     * 
     * @param knot The knot to find.
     * @return The index of the first occurrence of the knot in the knot-vector or the size of the
     *         knot-vector, if the knot is not found.
     */
    private int getIndexInKnotVector(final double knot) {
        final ListIterator<Double> iter = knotVector.listIterator();
        double currentKnot = knot - 1;

        while (iter.hasNext()) {
            currentKnot = iter.next();
            final double diff = Math.abs(currentKnot - knot);
            if (diff < EPSILON) {
                return iter.nextIndex() - 1;
            }
        }

        return knotVector.size();
    }

    /**
     * Calculates the multiplicity of a given knot.
     * 
     * @param knotToCheck The knot those multiplicity to check.
     * @return The multiplicity of the knot.
     */
    private int getMultiplicity(final double knotToCheck) {
        final Iterator<Double> iter = knotVector.listIterator();
        double currentKnot;
        int count = 0;

        while (iter.hasNext()) {
            currentKnot = iter.next();
            final double diff = currentKnot - knotToCheck;
            if (diff > EPSILON) {
                return count;
            } else if (diff > -EPSILON) {
                count++;
            }
        }

        return count;
    }

    // ########################################################################################
    // Assisting methods

    /**
     * Creates a uniform knot-vector from 0.0 to 1.0.
     * 
     * @param clamped Specifies if the knot-vector shall be clamped.
     * @param size The number of knots to generate.
     */
    private void createUniformKnotVector(final boolean clamped, final int size) {
        if (size < (2 * dimNUBS)) {
            throw new IllegalArgumentException(
                    "The knot vector must have at least two time the dimension elements.");
        }
        double mySize;

        if (clamped) {
            minKnot = 0.0;
            maxKnot = 1.0;
            for (int i = 0; i < dimNUBS; i++) {
                knotVector.add(0.0);
            }
            mySize = size + 1 - 2 * dimNUBS;
        } else {
            mySize = size + 1;
            final double ddim = (double) dimNUBS;
            minKnot = ddim / (mySize + 1);
            maxKnot = (mySize - ddim) / mySize;
        }

        final double fraction = mySize;
        for (int i = 1; i < mySize; i++) {
            knotVector.add((double) i / fraction);
        }

        if (isClamped) {
            for (int i = 0; i < dimNUBS; i++) {
                knotVector.add(1.0);
            }
        }
    }

    /**
     * Calculates the outerBox as the extreme values of all control points.
     */
    private void calculateOuterBox() {
        outerBox = new Rectangle(getControlPoints());
    }

    /**
     * Converts a polar vector to a single t.
     * 
     * @param polar The polar coordinate to convert.
     * @return A t value, representing the polar coordinate.
     */
    private static double getTFromPolar(final Collection<Double> polar) {
        double sum = 0.0;
        for (final Double val : polar) {
            sum += val;
        }
        
        return sum / polar.size();
    }

    /**
     * Inserts a new knots to the knot-vector. Also calculates the new list of CPs. There is a
     * limitation on the value of the inserted knot:
     * <UL>
     * <LI>minKnot < knotToInsert < maxKnot</LI>
     * <LI>multiplicity of knot to insert may not exceed the dimension of the Nubs.</LI>
     * </UL>
     * 
     * @param knotToInsert The knot to insert.
     */
    private void insertKnot(final double knotToInsert, final int insertions) {
        final ListIterator<Double> iterKnot = knotVector.listIterator();

        // Skip elements not referring to a point on the drawn spline.
        if (isClamped) {
            for (int i = 0; i < dimNUBS; i++) {
                iterKnot.next();
            }
        } else {
            for (int i = 0; i < dimNUBS - 1; i++) {
                iterKnot.next();
            }
        }

        final ListIterator<PolarCP> iterCP = controlPoints.listIterator();
        double currentKnot = iterKnot.next();

        
        // we always insert the new knot AFTER all identical knots.
        while (currentKnot - knotToInsert < EPSILON) {
            currentKnot = iterKnot.next();
            iterCP.next();
        }
        // One step back on the knot-vector
        iterKnot.previous();

        // Insert the new knot.
        insertKnotAtCurrentPosition(insertions, knotToInsert, iterCP, iterKnot);
    }

    /**
     * Inserts a new knot to the knot-vector and calculates the new list of CPs.
     * The correct position in knot-vector and list of controlPoints must be pre-calculated! 
     * 
     * 
     * @param insertions How many times to insert the knot?
     * @param knotToInsert The knot to insert.
     * @param iterCP A iterator over the CPs. Iterator must point to the position before the first CP,
     *            that will be used for the calculation.
     * @param iterKnot A iterator over the knot-vector. Iterator must point to position where to insert
     *            the new knot.
     */
    private void insertKnotAtCurrentPosition(final int insertions, final double knotToInsert,
            final ListIterator<PolarCP> iterCP, final ListIterator<Double> iterKnot) {
        final int multiplicity = getMultiplicity(knotToInsert);
        for (int i = 0; i < insertions; i++) {
            // Insert the new knot to the knotVector.
            iterKnot.add(knotToInsert);

            // We will first construct the new CPs and than add them.
            final List<PolarCP> newCPs = Lists.newArrayList();
            // The first CP we need for the calculation.
            PolarCP secondCP = iterCP.next();

            for (int j = multiplicity + i; j < dimNUBS; j++) {
                // The second CP we need for the calculation.
                final PolarCP firstCP = secondCP;
                secondCP = iterCP.next();

                // Calculate the new CP.
                newCPs.add(new PolarCP(firstCP, secondCP, knotToInsert));
            }

            // move to the insertion position, and on the way delete all CPs we have used for two
            // calculations as they don't belong to the new list of CPs.
            for (int j = multiplicity + i; j < dimNUBS; j++) {
                iterCP.previous();
                if (j > multiplicity + i) {
                    iterCP.remove();
                }
            }

            // now we can add the new CPs
            for (final PolarCP cp : newCPs) {
                iterCP.add(cp);
            }

            // Move back to the position in front of the first new CP, if there will be more
            // insertions
            if (i < insertions - 1) {
                for (int j = multiplicity + i; j < dimNUBS; j++) {
                    iterCP.previous();
                }
            }
        }
    }

    // ########################################################################################
    // Calculating extrema

    /**
     * Returns the t-value of the first vertical point on a NubSpline.
     * 
     * @param nubSpline The spline to check. 
     * @param accuracy The desired accuracy.
     * @param maxRecursion The maximum number of recursions.
     * @return The t-value of the maximum of the spline.
     */
    protected static KVector getFirstVerticalPoint(
            final NubSpline nubSpline, final double accuracy, final int maxRecursion) {
        
        final NubSpline firstDerive = generateDerivedNUBS(nubSpline);
        double currentAccuracy = Double.POSITIVE_INFINITY;
        KVector currentVector = null;
        int loopCount = 0;
        double knot = 0.0;
        
        while (currentAccuracy > accuracy && loopCount < maxRecursion) {
            knot = getZeroXOfControlPoligon(firstDerive);
            currentVector = firstDerive.getPointOnCurve(knot, true);
            currentAccuracy = Math.abs(currentVector.x);
            loopCount++;
        }
        return nubSpline.getPointOnCurve(knot, false);
    }
    
    /**
     * Returns the t-value of the last vertical point on a NubSpline. 
     * This is the value  
     * 
     * @param nubSpline The spline to check. 
     * @param accuracy The desired accuracy.
     * @param maxRecursion The maximum number of recursions.
     * @return The t-value of the maximum of the spline.
     */
    protected static KVector getLastVerticalPoint(
            final NubSpline nubSpline, final double accuracy, final int maxRecursion) {
        
        final NubSpline firstDerive = generateInvertedNUBS(generateDerivedNUBS(nubSpline));
        double currentAccuracy = Double.POSITIVE_INFINITY;
        KVector currentVector = null;
        int loopCount = 0;
        double knot = 0.0;
        
        while (currentAccuracy > accuracy && loopCount < maxRecursion) {
            knot = getZeroXOfControlPoligon(firstDerive);
            currentVector = firstDerive.getPointOnCurve(knot, true);
            currentAccuracy = Math.abs(currentVector.x);
            loopCount++;
        }
        final Double maxVal = nubSpline.getKnotVector().get(nubSpline.getKnotVector().size() - 1);
        return nubSpline.getPointOnCurve(maxVal - knot, false);
    }

    /**
     * Returns the t-value of the first horizontal point on a NubSpline.
     * 
     * @param nubSpline The spline to check. 
     * @param accuracy The desired accuracy.
     * @param maxRecursion The maximum number of recursions.
     * @return The t-value of the maximum of the spline.
     */
    protected static KVector getFirstHorizontalPoint(
            final NubSpline nubSpline, final double accuracy, final int maxRecursion) {
        
        final NubSpline firstDerive = generateDerivedNUBS(nubSpline);
        double currentAccuracy = Double.POSITIVE_INFINITY;
        KVector currentVector = null;
        int loopCount = 0;
        double knot = 0.0;
        
        while (currentAccuracy > accuracy && loopCount < maxRecursion) {
            knot = getZeroYOfControlPoligon(firstDerive);
            currentVector = firstDerive.getPointOnCurve(knot, true);
            currentAccuracy = Math.abs(currentVector.y);
            loopCount++;
        }
        return nubSpline.getPointOnCurve(knot, false);
    }
    
    /**
     * Returns the t-value of the last horizontal point on a NubSpline.
     * 
     * @param nubSpline The spline to check. 
     * @param accuracy The desired accuracy.
     * @param maxRecursion The maximum number of recursions.
     * @return The t-value of the maximum of the spline.
     */
    protected static KVector getLastHorizontalPoint(
            final NubSpline nubSpline, final double accuracy, final int maxRecursion) {
        
        final NubSpline firstDerive = generateInvertedNUBS(generateDerivedNUBS(nubSpline));
        double currentAccuracy = Double.POSITIVE_INFINITY;
        KVector currentVector = null;
        int loopCount = 0;
        double knot = 0.0;
        
        while (currentAccuracy > accuracy && loopCount < maxRecursion) {
            knot = getZeroYOfControlPoligon(firstDerive);
            currentVector = firstDerive.getPointOnCurve(knot, true);
            currentAccuracy = Math.abs(currentVector.y);
            loopCount++;
        }
        final Double maxVal = nubSpline.getKnotVector().get(nubSpline.getKnotVector().size() - 1);
        return nubSpline.getPointOnCurve(maxVal - knot, false);
    }

    /**
     * Calculates a point on this NubSpline at the position of value. This version does not modify this
     * Nubs. 
     * There is a limitation on the value
     * to check:
     * <UL>
     * <LI>If the NubSpline is clamped:</LI>
     * <UL>
     * <LI>knotVector[0] <= t <= knotVector[knotVector.size() - 1]</LI>
     * </UL>
     * <LI>If the NubSpline is unclamped:</LI>
     * <UL>
     * <LI>knotVector[dim - 1] <= t <= knotVector[knotVector.size() - dim]</LI>
     * 
     * @param t The position on this NubSpline to calculate. Must be in the range of the knot vector.
     * @return The position on this NubSpline.
     */
    protected KVector getPointOnCurve(final double t) {
        return getPointOnCurve(t, false);
    }
    

    /**
     * Calculates a point on this NubSpline at the position of value. There is a limitation on the value
     * to check:
     * <UL>
     * <LI>If the NubSpline is clamped:</LI>
     * <UL>
     * <LI>knotVector[0] <= t <= knotVector[knotVector.size() - 1]</LI>
     * </UL>
     * <LI>If the NubSpline is unclamped:</LI>
     * <UL>
     * <LI>knotVector[dim - 1] <= t <= knotVector[knotVector.size() - dim]</LI>
     * 
     * @param t The position on this NubSpline to calculate. Must be in the range of the knot vector.
     * @param modify Determines, if current Nubs is allowed to be modified.
     * @return The position on this NubSpline.
     */
    protected KVector getPointOnCurve(final double t, final boolean modify) {
        final int multiplicity = getMultiplicity(t);
        // If the multiplicity is already equal to dim, there is already a CP at the position t.
        if (multiplicity == dimNUBS) {
            return getControlPoint(getIndexInKnotVector(t));
        }

        // If the multiplicity is less than dim, we have to insert new knots. 
        if (modify) {
            insertKnot(t, dimNUBS - multiplicity);
            return getControlPoint(getIndexInKnotVector(t));
        } else {
            // To not modify this NubSpline, we will do so on a copy.
            final NubSpline copy = new NubSpline(this);
            copy.insertKnot(t, dimNUBS - multiplicity);
            return copy.getControlPoint(copy.getIndexInKnotVector(t));
        }
    }
    

    /**
     * Returns a progression value (t) of the first zero (Y-value) in the control polygon. 
     * If the polygon starts at zero, this zero is ignored and the next zero is returned. 
     * If there is no zero found, 0.0 is returned.
     * 
     * @param nubSpline The NubSpline those zero of control polygon to calculate.
     * @return The progression value t. min[knotVector] < t <= max[knotVector]
     */
    private static double getZeroYOfControlPoligon(final NubSpline nubSpline) {
        final Iterator<PolarCP> iter = nubSpline.controlPoints.iterator();
        // the CPs we are working with.
        PolarCP firstCP, secondCP;
        // Determines if one of the CPs is bigger than zero.
        boolean firstIsPositive, secondIsPositive;
        // Determines if one of the CPs is smaller than zero.
        boolean firstIsNegative, secondIsNegative;
        // the y values of the CPs we are working with.
        double firstY = 0.0, secondY;
        
        // initialize
        secondCP = iter.next();
        secondY = secondCP.getCp().y;
        secondIsPositive = secondY > EPSILON;
        secondIsNegative = secondY < -EPSILON;

        while (iter.hasNext()) {
            // Proceed to next line between CPs:
            //update first CP.
            firstCP = secondCP;
            firstY = secondY;
            firstIsPositive = secondIsPositive;
            firstIsNegative = secondIsNegative;

            // update second CP.
            secondCP = iter.next();
            secondY = secondCP.getCp().y;
            secondIsPositive = secondY > EPSILON;
            secondIsNegative = secondY < -EPSILON;

            // Check if there is a zero at the second CP.
            if (!(secondIsPositive || secondIsNegative)) {
                return getTFromPolar(secondCP.getPolarCoordinate());
            }

            // check for zero in the line between the two CPs
            if (firstIsPositive && secondIsNegative || firstIsNegative && secondIsPositive) {
                double factor;
                factor = firstY / (firstY - secondY);
                final double t1 = getTFromPolar(firstCP.getPolarCoordinate());
                final double t2 = getTFromPolar(secondCP.getPolarCoordinate());
                return factor * t1 + (1 - factor) * t2;
            }
        }
        return 0.0;
  }    
    
    /**
     * Returns a progression value (t) of the first zero (X-value) in the control polygon. If the polygon
     * starts at zero, this zero is ignored and the next zero is returned. If there is no zero found, 
     * 0.0 is returned.
     * 
     * @param nubSpline The NubSpline those zero of control polygon to calculate.
     * @return The progression value t. min[knotVector] < t <= max[knotVector]
     */
    private static double getZeroXOfControlPoligon(final NubSpline nubSpline) {
        final Iterator<PolarCP> iter = nubSpline.controlPoints.iterator();
        // the CPs we are working with.
        PolarCP firstCP, secondCP;
        // Determines if one of the CPs is bigger than zero.
        boolean firstIsPositive, secondIsPositive;
        // Determines if one of the CPs is smaller than zero.
        boolean firstIsNegative, secondIsNegative;
        // the x values of the CPs we are working with.
        double firstX = 0.0, secondX;
        
        // initialize
        secondCP = iter.next();
        secondX = secondCP.getCp().x;
        secondIsPositive = secondX > EPSILON;
        secondIsNegative = secondX < -EPSILON;

        while (iter.hasNext()) {
            // Proceed to next line between CPs:
            // Update first CP.
            firstCP = secondCP;
            firstX = secondX;
            firstIsPositive = secondIsPositive;
            firstIsNegative = secondIsNegative;
            
            // Update second CP.
            secondCP = iter.next();
            secondX = secondCP.getCp().x;
            secondIsPositive = secondX > EPSILON;
            secondIsNegative = secondX < -EPSILON;

            // Check if there is a zero at the second CP.
            if (!(secondIsPositive || secondIsNegative)) {
                return getTFromPolar(secondCP.getPolarCoordinate());
            }

            // check for zero in the line between the two CPs
            if (firstIsPositive && secondIsNegative || firstIsNegative && secondIsPositive) {
                double factor;
                factor = firstX / (firstX - secondX);
                final double t1 = getTFromPolar(firstCP.getPolarCoordinate());
                final double t2 = getTFromPolar(secondCP.getPolarCoordinate());
                return factor * t1 + (1 - factor) * t2;
            }
        }
        return 0.0;
  }

    // ########################################################################################
    // Bezier

    /**
     * Converts this NubSpline to a bezier spline. All inner knots of the knotVector get the multiplicity
     * of dimNUBS.
     */
    public void toBezier() {
        final ListIterator<Double> iterKnot = knotVector.listIterator();
        final ListIterator<PolarCP> iterCP = controlPoints.listIterator();

        // Unclamped knotVectors have (dim-1) leading and trailing knots that are not repeated and
        // don't need to be repeated.
        // Clamped knotVectors have (dim) leading and trailing knots that are already repeated. We
        // skip them for performance.
        if (isClamped) {
            for (int i = 0; i < dimNUBS; i++) {
                iterKnot.next();
            }
        } else {
            for (int i = 0; i < dimNUBS - 1; i++) {
                iterKnot.next();
                iterKnot.remove();
            }
        }

        double currentKnot = iterKnot.next();
        // Iterate over all knots those multiplicity we possibly have to increase.
        while (maxKnot - currentKnot > EPSILON) {
            final double knotToCount = currentKnot;
            int occurrence = 0;

            // Count occurrences of knotToCount.
            while (Math.abs(currentKnot - knotToCount) < EPSILON) {
                occurrence++;
                currentKnot = iterKnot.next();
                iterCP.next();
            }

            // insert new knots, if multiplicity is not as expected (dimNUBS)
            if (occurrence < dimNUBS) {
                iterKnot.previous();
                insertKnotAtCurrentPosition(dimNUBS - occurrence, knotToCount, iterCP, iterKnot);
                iterKnot.next();
            }

            // Proceed to next elements.
            iterCP.previous();
        }

        if (!isClamped) {
            for (int i = 0; i < dimNUBS - 1; i++) {
                iterKnot.next();
                iterKnot.remove();
            }
        }
        isClamped = true;
        isBezier = true;
    }

    /**
     * Returns all bezier control points needed to create an equivalent bezier spline. If this NubSpline
     * still not holds the requirements for a bezier spline, it is converted to to so.
     * 
     * @param withSourceVector
     *            Specify if the source vector shall be included in the vector chain.
     * @param withTargetVector
     *            Specify if the target vector shall be included in the vector chain.
     * @return The bezier control-points.
     */
    public KVectorChain getBezierCP(final boolean withSourceVector, final boolean withTargetVector) {
        if (!isBezier) {
            toBezier();
        }
        final KVectorChain retVal = new KVectorChain();
        final Iterator<PolarCP> iter = controlPoints.iterator();

        if (!withSourceVector) {
            iter.next();
        }

        while (iter.hasNext()) {
            retVal.add(iter.next().getCp());
        }

        if (!withTargetVector) {
            retVal.removeLast();
        }

        return retVal;
    }

    /**
     * Returns all bezier control points needed to create an equivalent bezier spline. If this NubSpline
     * still not holds the requirements for a bezier spline, it is converted to to so. The source and
     * target vector are not included in the list of control points.
     * 
     * @return The bezier control points without source and target vector.
     */
    public KVectorChain getBezierCP() {
        return getBezierCP(false, false);
    }
    
    // ########################################################################################
    // Other stuff

    /**
     * Generates NubSpline CPs for a nice curve.
     * @return The nice CPs
     */
    public static NubSpline generateNiceCurve() {
        final KVectorChain niceChain = new KVectorChain();
        // CHECKSTYLEOFF Magic Numbers    
        final KVector vector1 = new KVector(5.5, 23.0);
        final KVector vector2 = new KVector(2.5, 12.0);
        final KVector vector3 = new KVector(5.0, 10.0);
        final KVector vector4 = new KVector(5.0, 9.0);
        final KVector vector5 = new KVector(3.8, 5.5);
        final KVector vector6 = new KVector(7.0, 4.0);
        final KVector vector7 = new KVector(7.0, 3.5);
        final KVector vector8 = new KVector(6.0, 2.2);
        final KVector vector9 = new KVector(8.0, 0.5);
        //CHECKSTYLE.ON: Magic Numbers     

        niceChain.add(vector1);
        niceChain.add(vector1);
        niceChain.add(vector1);
        niceChain.add(vector2);
        niceChain.add(vector3);
        niceChain.add(vector4);
        niceChain.add(vector5);
        niceChain.add(vector6);
        niceChain.add(vector7);
        niceChain.add(vector8);
        niceChain.add(vector9);
        niceChain.add(mirrorOnX(vector8, vector9.x));
        niceChain.add(mirrorOnX(vector7, vector9.x));
        niceChain.add(mirrorOnX(vector6, vector9.x));
        niceChain.add(mirrorOnX(vector5, vector9.x));
        niceChain.add(mirrorOnX(vector4, vector9.x));
        niceChain.add(mirrorOnX(vector3, vector9.x));
        niceChain.add(mirrorOnX(vector2, vector9.x));
        niceChain.add(mirrorOnX(vector1, vector9.x));
        niceChain.add(mirrorOnX(vector1, vector9.x));
        niceChain.add(mirrorOnX(vector1, vector9.x));
        
        niceChain.scale(5.0);
        niceChain.offset(-50.0, -80.0);

        // Create the new NubSpline.
        return new NubSpline(true, 3, niceChain);
    }
    
    /**
     * Creates a copy of the original KVector, mirrored on a line laying on the given x-value.
     * @param original The original KVector.
     * @param xVal The mirror, laying parallel to the y-axis, at the given x-value.
     * @return
     */
    private static KVector mirrorOnX(final KVector original, final double xVal) {
        return new KVector(xVal + xVal - original.x, original.y);
    }
    
    @Override
    public String toString() {
        return controlPoints.toString();
    }
        
    // ########################################################################################
    // Class for polar form of control-points.
    
    /**
     * Represents a control point of this NubSpline including it's polar coordinate.
     * 
     * @author tit
     * 
     */
    protected static final class PolarCP {
        /** The control point. */
        private KVector cp;
        /** The polar coordinate. (size is always = dim) */
        private List<Double> polarCoordinate;

        /**
         * Constructs a new PolarCP from the given parameters.
         * 
         * @param controlPoint The control-point of the PolarCP.
         * @param polarCoordinate The polar-coordinate of the PolarCP.
         */
        PolarCP(final KVector controlPoint, final List<Double> polarCoordinate) {
            setCp(controlPoint.clone());
            setPolarCoordinate(Lists.newLinkedList(polarCoordinate));
        }

        /**
         * Copy constructor. All fields are created newly.
         * 
         * @param polarCP The PolarCP to copy.
         */
        PolarCP(final PolarCP polarCP) {
            setCp(polarCP.getCp().clone());
            setPolarCoordinate(Lists.newLinkedList(polarCP.getPolarCoordinate()));
        }

        /**
         * This Constructor is calculating the new PolarCP from two given once. There are
         * requirements on the polarCoordinates that must be met, but are not checked (for
         * performance):
         * <ul>
         * <li>firstCP.polarCoordinates = (xW)</li>
         * <li>secondCP.polarCoordinates = (Wy)</li>
         * <li>resulting polarCoordinate = (Wz)</li>
         * </ul>
         * Also (Wz) is properly ordered.
         * 
         * @param firstCP First polarCP.
         * @param secondCP Second polarCP.
         * @param newKnot The new knotValue.
         */
        PolarCP(final PolarCP firstCP, final PolarCP secondCP, final double newKnot) {
            final double firstFactor = firstCP.polarCoordinate.iterator().next();
            final double secondFactor = Iterables.getLast(secondCP.polarCoordinate);

            // vectorC = ((b-c) * vectorA + (c-a) * vectorB) / (b-a)
            final KVector aScaled = firstCP.cp.clone().scale(secondFactor - newKnot);
            final KVector bScaled = secondCP.cp.clone().scale(newKnot - firstFactor);
            final KVector total = aScaled.add(bScaled);
            total.scale(1.0 / (secondFactor - firstFactor));

            cp = total;
            polarCoordinate = Lists.newArrayList();

            // Specifies, if the newKnot still needs to be added.
            boolean needsToBeAdded = true;

            // From the firstCP, we will take the W. (See method comment)
            final Iterator<Double> iter = firstCP.polarCoordinate.iterator();
            iter.next();
            while (iter.hasNext()) {
                final double nextKnot = iter.next();
                if (needsToBeAdded && nextKnot - newKnot > EPSILON) {
                    polarCoordinate.add(newKnot);
                    needsToBeAdded = false;
                }
                polarCoordinate.add(nextKnot);
            }
            if (needsToBeAdded) {
                polarCoordinate.add(newKnot);
            }
        }

        /**
         * Returns the current control-point. May be modified.
         * 
         * @return The control-point.
         */
        public KVector getCp() {
            return cp;
        }

        /**
         * Sets the control-point. Vector is not copied!
         * 
         * @param cp
         *            The new control-point.
         */
        public void setCp(final KVector cp) {
            this.cp = cp;
        }

        /**
         * Returns the polar-coordinate.
         * 
         * @return The polar-coordinate. May be modified.
         */
        public List<Double> getPolarCoordinate() {
            return polarCoordinate;
        }

        /**
         * Sets the polar-coordinate.
         */
        public void setPolarCoordinate(final List<Double> polarCoordinate) {
            this.polarCoordinate = polarCoordinate;
        }

        // elkjs-exclude-start
        @Override
        public String toString() {
            return polarCoordinate + " " + SplinesMath.convertKVectorToString(cp);
        }
        // elkjs-exclude-end
    }
}
