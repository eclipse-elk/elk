/*******************************************************************************
 * Copyright (c) 2018, 2020 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.rectpacking.options;

/**
 * Specifies the strategy employed to pack rectangles sequentially during {@link AreaApproximation}.
 */
public enum OptimizationGoal {
    /** Aspect ratio-driven packing heuristic. */
    ASPECT_RATIO_DRIVEN,
    /** Max scale-driven packing heuristic. */
    MAX_SCALE_DRIVEN,
    /** Area driven packing heuristic. */
    AREA_DRIVEN;
}
