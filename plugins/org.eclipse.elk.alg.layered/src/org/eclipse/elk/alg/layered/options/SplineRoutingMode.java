/*******************************************************************************
 * Copyright (c) 2018 Kiel University and others.
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
package org.eclipse.elk.alg.layered.options;

/**
 * Specifies the strategy to be used for creating spline edge routes.
 */
public enum SplineRoutingMode {

    /** Uses computed long edge dummies as hint for the spline paths. Ensures that the splines 
     * do not overlap with nodes. On the downside, the spline paths feel rather orthogonal-ish.
     */
    CONSERVATIVE,
    /** Basically the same as {@link #CONSERVATIVE}. Uses softer curves where edges attach to the nodes,
     *  at the risk of overlapping other graph elements in that area. */
    CONSERVATIVE_SOFT,
    /** Still uses computed long edge dummies as hints but uses far less control points to define the spline paths.
     *  While this makes the splines curvier, the chances of splines overlapping nodes increas significantly.
     */
    SLOPPY,
}
