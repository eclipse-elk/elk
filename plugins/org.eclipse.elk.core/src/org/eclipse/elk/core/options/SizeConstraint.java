/*******************************************************************************
 * Copyright (c) 2012, 2015 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.core.options;

import java.util.EnumSet;

/**
 * Things to take into account when determining the size of a node. Each item of this enumeration
 * corresponds to something a layout algorithm should pay attention to when calculating node sizes.
 * Usually, one will use a combination of these values in an {@code EnumSet} instance, with the
 * empty set meaning that node sizes are fixed. <b>This enumeration is not set directly on
 * {@link CoreOptions#SIZE_CONSTRAINT}; instead, an {@code EnumSet} over this enumeration is used
 * there.</b>
 * 
 * <p><i>Note:</i> Layout algorithms may only support a subset of these options.</p>
 *
 * @see SizeOptions
 * @author msp
 * @author cds
 */
public enum SizeConstraint {
    
    /**
     * The number of ports and their position should be taken into account when determining the
     * size of nodes. The sum of port widths and heights and the minimum spacing between ports
     * is a lower bound for the node size.
     */
    PORTS,
    
    /**
     * Ports labels are taken into account when determining the size of nodes. Depending on where
     * the labels are positioned the node will be made large enough to avoid overlaps and to try
     * to place labels in as unambiguous a way as possible. Setting this option doesn't make any
     * sense if the {@link #PORTS} option is not set as well.
     */
    PORT_LABELS,
    
    /**
     * A node's labels are taken into account.
     */
    NODE_LABELS,
    
    /**
     * If set, a node's size will be at least the minimum size set on it. If no minimum size is set,
     * the behavior depends on whether the {@link SizeOptions#DEFAULT_MINIMUM_SIZE} constraint is
     * set as well.
     */
    MINIMUM_SIZE;
    
    
    /**
     * Returns an empty enum set over this enumeration, which corresponds to fixed size constraints.
     * 
     * @return set over this enumeration representing fixed size constraints.
     */
    public static EnumSet<SizeConstraint> fixed() {
        return EnumSet.noneOf(SizeConstraint.class);
    }
    
    /**
     * Returns a set containing the {@link #MINIMUM_SIZE} constraint.
     * 
     * @return set with minimum size constraint.
     */
    public static EnumSet<SizeConstraint> minimumSize() {
        return EnumSet.of(MINIMUM_SIZE);
    }
    
    /**
     * Returns a set containing the common combination of {@link #MINIMUM_SIZE} and {@link #PORTS}.
     * 
     * @return set with minimum size constraint in combination with ports.
     */
    public static EnumSet<SizeConstraint> minimumSizeWithPorts() {
        return EnumSet.of(PORTS, MINIMUM_SIZE);
    }
    
    /**
     * Returns a set containing all options defined in this enumeration, effectively giving
     * the layout algorithm as much freedom as possible in determining the node size.
     * 
     * @return set with all available options.
     */
    public static EnumSet<SizeConstraint> free() {
        return EnumSet.allOf(SizeConstraint.class);
    }
    
}
