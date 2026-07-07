/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.disco.graph;

/**
 * Represents the four cardinal directions.
 */
public enum DCDirection {
    /**
     * The four cardinal directions.
     */
    NORTH, EAST, SOUTH, WEST;

    private boolean horizontal;

    static {
        NORTH.horizontal = false;
        EAST.horizontal = true;
        SOUTH.horizontal = false;
        WEST.horizontal = true;
    }

    /**
     * Returns whether the given direction is horizontal.
     * 
     * @return true, if {@code EAST} or {@code WEST}; false, otherwise.
     */
    public boolean isHorizontal() {
        return horizontal;
    }
}
