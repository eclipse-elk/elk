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
package org.eclipse.elk.core.options;

/**
 * Definition of layout directions. To be accessed using {@link CoreOptions#DIRECTION}.
 * 
 * @author msp
 */
public enum Direction {

    /** undefined layout direction. */
    UNDEFINED,
    /** rightward horizontal layout. */
    RIGHT,
    /** leftward horizontal layout. */
    LEFT,
    /** downward vertical layout. */
    DOWN,
    /** upward vertical layout. */
    UP;
    
    
    /**
     * Checks if this layout direction is horizontal. (that is, left or right) An undefined layout
     * direction is not horizontal.
     * 
     * @return {@code trure} if the layout direction is horizontal.
     */
    public boolean isHorizontal() {
        return this == LEFT || this == RIGHT;
    }

    /**
     * Checks if this layout direction is vertical. (that is, up or down) An undefined layout
     * direction is not vertical.
     * 
     * @return {@code trure} if the layout direction is vertical.
     */
    public boolean isVertical() {
        return this == UP || this == DOWN;
    }
    
    /**
     * @return the opposite direction of {@code this}. For instance, if this is {@link #LEFT},
     *         return {@link #RIGHT}.
     */
    public Direction opposite() {
        switch (this) {
        case LEFT:
            return RIGHT;
        case RIGHT:
            return LEFT;
        case UP:
            return DOWN;
        case DOWN:
            return UP;
        default:
            return UNDEFINED;
        }
    }
    
}
