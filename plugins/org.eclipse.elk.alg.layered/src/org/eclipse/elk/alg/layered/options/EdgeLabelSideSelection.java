/*******************************************************************************
 * Copyright (c) 2012, 2017 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.layered.options;

/**
 * Possible ways to determine whether edge labels should be placed above or below their edge.
 * 
 * @see org.eclipse.elk.alg.layered.intermediate.LabelSideSelector
 */
public enum EdgeLabelSideSelection {
    
    /** Labels are always placed above their edge. */
    ALWAYS_UP,
    /** Labels are always placed below their edge. */
    ALWAYS_DOWN,
    /** Labels are always placed above their edge if the edge points rightwards, or below if it doesn't. */
    DIRECTION_UP,
    /** Labels are always placed below their edge if the edge points leftwards, or above if it doesn't. */
    DIRECTION_DOWN,
    /** Determine the side automatically, default to UP if in doubt. */
    SMART_UP,
    /** Determine the side automatically, default to DOWN if in doubt. */
    SMART_DOWN;
    
    /**
     * Transposes this label side selection strategy and returns the result. Downward selections are turned upwards
     * while upward selections are turned downwards.
     *
     * @return the transformed edge label side selection.
     */
    public EdgeLabelSideSelection transpose() {
        switch (this) {
        case ALWAYS_UP:
            return ALWAYS_DOWN;
        case ALWAYS_DOWN:
            return ALWAYS_UP;
        case DIRECTION_UP:
            return DIRECTION_DOWN;
        case DIRECTION_DOWN:
            return DIRECTION_UP;
        case SMART_UP:
            return SMART_DOWN;
        case SMART_DOWN:
            return SMART_UP;
        default:
            assert false;
            return null;
        }
    }

}
