/*******************************************************************************
 * Copyright (c) 2019 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.layered.intermediate.loops;

import java.util.Set;

import org.eclipse.elk.core.options.PortSide;

/**
 * The different types of self loops.
 */
public enum SelfLoopType {
    /** Connects ports that are all on the same side. */
    ONE_SIDE,
    /** Connects ports on two adjacent sides. */
    TWO_SIDES_CORNER,
    /** Connects ports on two opposing sides. */
    TWO_SIDES_OPPOSING,
    /** Connects ports spread out over three sides. */
    THREE_SIDES,
    /** Connects ports spread out over all four sides. */
    FOUR_SIDES;
    
    /**
     * Determines the self loop type for a self loop involving ports on the given sides. The result only depends on the
     * port configuration. In theory, a one-side self loop could be routed around the whole node to reach its target
     * port, thus causing it to occupy all four sides.
     */
    public static SelfLoopType fromPortSides(final Set<PortSide> portSides) {
        if (portSides.contains(PortSide.UNDEFINED)) {
            throw new IllegalArgumentException("Port sides must not contain UNDEFINED");
        }
        
        // CHECKSTYLEOFF MagicNumber
        switch (portSides.size()) {
        case 1:
            return ONE_SIDE;
            
        case 2:
            // Check if we have opposing sides or not
            boolean eastWest = portSides.contains(PortSide.EAST) && portSides.contains(PortSide.WEST);
            boolean northSouth = portSides.contains(PortSide.NORTH) && portSides.contains(PortSide.SOUTH);
            
            if (eastWest || northSouth) {
                return TWO_SIDES_OPPOSING;
            } else {
                return TWO_SIDES_CORNER;
            }
            
        case 3:
            return THREE_SIDES;
            
        case 4:
            return FOUR_SIDES;
            
        default:
            return null;
        }
        // CHECKSTYLEON MagicNumber
    }
}