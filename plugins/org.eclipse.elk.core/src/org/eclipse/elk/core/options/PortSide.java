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

import java.util.Collections;
import java.util.EnumSet;
import java.util.Set;

import com.google.common.collect.Sets;

/**
 * Definition of port sides on a node. Besides defining the actual available port sides, this enumeration
 * also defines convenience methods as well as constants with all possible port side combinations. The
 * latter are named in a clockwise manner, starting at the northern side.
 * 
 * @author msp
 * @author cds
 */
public enum PortSide {
    
    /** the side is undefined. */
    UNDEFINED,
    /** top side. */
    NORTH,
    /** right side. */
    EAST,
    /** bottom side. */
    SOUTH,
    /** left side. */
    WEST;
    
    
    //////////////////////////////////////////////////////////////////////////////////////////////
    // Port Side Combinations

    /** Immutable set containing the given port sides. */
    public static final Set<PortSide> SIDES_NONE = Collections.unmodifiableSet(
            EnumSet.noneOf(PortSide.class));
    /** Immutable set containing the given port sides. */
    public static final Set<PortSide> SIDES_NORTH = Sets.immutableEnumSet(
            PortSide.NORTH);
    /** Immutable set containing the given port sides. */
    public static final Set<PortSide> SIDES_EAST = Sets.immutableEnumSet(
            PortSide.EAST);
    /** Immutable set containing the given port sides. */
    public static final Set<PortSide> SIDES_SOUTH = Sets.immutableEnumSet(
            PortSide.SOUTH);
    /** Immutable set containing the given port sides. */
    public static final Set<PortSide> SIDES_WEST = Sets.immutableEnumSet(
            PortSide.WEST);
    /** Immutable set containing the given port sides. */
    public static final Set<PortSide> SIDES_NORTH_SOUTH = Sets.immutableEnumSet(
            PortSide.NORTH, PortSide.SOUTH);
    /** Immutable set containing the given port sides. */
    public static final Set<PortSide> SIDES_EAST_WEST = Sets.immutableEnumSet(
            PortSide.EAST, PortSide.WEST);
    /** Immutable set containing the given port sides. */
    public static final Set<PortSide> SIDES_NORTH_WEST = Sets.immutableEnumSet(
            PortSide.NORTH, PortSide.WEST);
    /** Immutable set containing the given port sides. */
    public static final Set<PortSide> SIDES_NORTH_EAST = Sets.immutableEnumSet(
            PortSide.NORTH, PortSide.EAST);
    /** Immutable set containing the given port sides. */
    public static final Set<PortSide> SIDES_SOUTH_WEST = Sets.immutableEnumSet(
            PortSide.SOUTH, PortSide.WEST);
    /** Immutable set containing the given port sides. */
    public static final Set<PortSide> SIDES_EAST_SOUTH = Sets.immutableEnumSet(
            PortSide.EAST, PortSide.SOUTH);
    /** Immutable set containing the given port sides. */
    public static final Set<PortSide> SIDES_NORTH_EAST_WEST = Sets.immutableEnumSet(
            PortSide.NORTH, PortSide.EAST, PortSide.WEST);
    /** Immutable set containing the given port sides. */
    public static final Set<PortSide> SIDES_EAST_SOUTH_WEST = Sets.immutableEnumSet(
            PortSide.EAST, PortSide.SOUTH, PortSide.WEST);
    /** Immutable set containing the given port sides. */
    public static final Set<PortSide> SIDES_NORTH_SOUTH_WEST = Sets.immutableEnumSet(
            PortSide.NORTH, PortSide.SOUTH, PortSide.WEST);
    /** Immutable set containing the given port sides. */
    public static final Set<PortSide> SIDES_NORTH_EAST_SOUTH = Sets.immutableEnumSet(
            PortSide.NORTH, PortSide.EAST, PortSide.SOUTH);
    /** Immutable set containing the given port sides. */
    public static final Set<PortSide> SIDES_NORTH_EAST_SOUTH_WEST = Sets.immutableEnumSet(
            PortSide.NORTH, PortSide.EAST, PortSide.SOUTH, PortSide.WEST);
    
    
    //////////////////////////////////////////////////////////////////////////////////////////////
    // Convenience Methods
    
    /**
     * Returns the next side in clockwise order.
     * 
     * @return the next side in clockwise order
     */
    public PortSide right() {
        switch (this) {
        case NORTH:
            return EAST;
        case EAST:
            return SOUTH;
        case SOUTH:
            return WEST;
        case WEST:
            return NORTH;
        default:
            return UNDEFINED;
        }
    }
    
    /**
     * Returns the next side in counter-clockwise order.
     * 
     * @return the next side in counter-clockwise order
     */
    public PortSide left() {
        switch (this) {
        case NORTH:
            return WEST;
        case EAST:
            return NORTH;
        case SOUTH:
            return EAST;
        case WEST:
            return SOUTH;
        default:
            return UNDEFINED;
        }
    }
    
    /**
     * Returns the opposed side.
     * 
     * @return the opposed side
     */
    public PortSide opposed() {
        switch (this) {
        case NORTH:
            return SOUTH;
        case EAST:
            return WEST;
        case SOUTH:
            return NORTH;
        case WEST:
            return EAST;
        default:
            return UNDEFINED;
        }
    }
    
    /**
     * Returns whether this and the given port side are adjacent port sides in that they share a corner of their node.
     * 
     * @param other the port side to check with.
     * @return {@code true} if the two sides share a corner.
     */
    public boolean areAdjacent(final PortSide other) {
        if (this == UNDEFINED) {
            return false;
        } else {
            return this.left() == other || this.right() == other;
        }
    }
    
    /**
     * Get the port side that corresponds to the given direction.
     * 
     * @param direction a direction
     * @return the corresponding port side
     */
    public static PortSide fromDirection(final Direction direction) {
        switch (direction) {
        case UP:
            return NORTH;
        case RIGHT:
            return EAST;
        case DOWN:
            return SOUTH;
        case LEFT:
            return WEST;
        default:
            return UNDEFINED;
        }
    }
    
    public static boolean isVertical(final PortSide side) {
        return side == NORTH || side == SOUTH;
    }
    
    public static boolean isHorizontal(final PortSide side) {
        return side == WEST || side == EAST;
    }
    
}
