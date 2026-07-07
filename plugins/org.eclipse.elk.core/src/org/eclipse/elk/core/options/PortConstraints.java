/*******************************************************************************
 * Copyright (c) 2009, 2015 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.core.options;

/**
 * Definition of port constraints. To be accessed using {@link CoreOptions#PORT_CONSTRAINTS}.
 * 
 * @author msp
 */
public enum PortConstraints {

    /** Undefined constraints. */
    UNDEFINED,
    /** All ports are free. */
    FREE,
    /** The side is fixed for each port. */
    FIXED_SIDE,
    /** The side is fixed for each port, and the order of ports is fixed for each side. */
    FIXED_ORDER,
    /**
     * The side is fixed for each port, the order or ports is fixed for each side and
     * the relative position of each port must be preserved. That means if the node is
     * resized by factor x, the port's position must also be scaled by x.
     */
    FIXED_RATIO,
    /** The exact position is fixed for each port. */
    FIXED_POS;
    
    
    /**
     * Returns whether the position of the ports is fixed. Note that this is not true
     * if port ratios are fixed.
     * 
     * @return true if the position is fixed
     */
    public boolean isPosFixed() {
        return this == FIXED_POS;
    }
    
    /**
     * Returns whether the ratio of port positions is fixed. Note that this is not true
     * if the port positions are fixed.
     * 
     * @return true if the ratio is fixed
     */
    public boolean isRatioFixed() {
        return this == FIXED_RATIO;
    }
    
    /**
     * Returns whether the order of ports is fixed.
     * 
     * @return true if the order of ports is fixed
     */
    public boolean isOrderFixed() {
        return this == FIXED_ORDER || this == FIXED_RATIO || this == FIXED_POS;
    }
    
    /**
     * Returns whether the sides of ports are fixed.
     * 
     * @see PortSide
     * @return true if the port sides are fixed
     */
    public boolean isSideFixed() {
        return this != FREE && this != UNDEFINED;
    }
    
}
