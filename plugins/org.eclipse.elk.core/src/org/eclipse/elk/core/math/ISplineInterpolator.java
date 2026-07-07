/*******************************************************************************
 * Copyright (c) 2010, 2015 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.core.math;

import java.util.LinkedList;

/**
 * Global interface for any Spline interpolator.
 * 
 * @author uru
 */
public interface ISplineInterpolator {

    /**
     * returns a piecewise bezierspline.
     * 
     * @param points
     *            as an array, see implementing class if to prefer Vector or List implementation.
     * @return piecewise bezierspline
     */
    BezierSpline interpolatePoints(KVector[] points);

    /**
     * returns a piecewise bezierspline.
     * 
     * @param points
     *            as an array, see implementing class if to prefer Vector or List implementation.
     * @param startVec
     *            tangent vector specifying to head out of the first node
     * @param endVec
     *            tangent vector specifying to head into the last node
     * @param tangentScale
     *            if true, the tangent is scaled depending on the distance to the next ctr point, if
     *            false the tangent is used as passed
     * @return piecewise bezierspline
     */
    BezierSpline interpolatePoints(KVector[] points, KVector startVec, KVector endVec,
            boolean tangentScale);

    /**
     * returns a piecewise bezierspline.
     * 
     * @param points
     *            as an array, see implementing class if to prefer Vector or List implementation.
     * @return piecewise bezierspline
     */
    BezierSpline interpolatePoints(LinkedList<KVector> points);

    /**
     * returns a piecewise bezierspline.
     * 
     * @param points
     *            as an array, see implementing class if to prefer Vector or List implementation.
     * @param startVec
     *            tangent vector specifying to head out of the first node
     * @param endVec
     *            tangent vector specifying to head into the last node
     * @param tangentScale
     *            if true, the tangent is scaled depending on the distance to the next ctr point, if
     *            false the tangent is used as passed
     * @return piecewise bezierspline
     */
    BezierSpline interpolatePoints(LinkedList<KVector> points, KVector startVec, KVector endVec,
            boolean tangentScale);
}
