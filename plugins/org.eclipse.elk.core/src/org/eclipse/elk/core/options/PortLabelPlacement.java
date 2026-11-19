/*******************************************************************************
 * Copyright (c) 2012, 2020 Kiel University and others.
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

import java.util.EnumSet;
import java.util.Set;

import com.google.common.collect.Sets;

/**
 * Options for controlling how port labels are placed by layout algorithms. The corresponding layout
 * option will usually accept an {@link EnumSet} over this enumeration, theoretically allowing
 * arbitrary and even contradictory subsets of options to be set. Note that you are restricted to
 * the following combinations if you want your choice to make sense:
 * <ul>
 *   <li>Exactly one of the {@link #INSIDE} and {@link #OUTSIDE} options (if neither is set, 
 *       label positions are not computed but left untouched).</li>
 *   <li>Exactly one of the {@link #ALWAYS_SAME_SIDE}, {@link #SPACE_EFFICIENT} options.</li>
 * </ul>
 * 
 * Additionally {@link #NEXT_TO_PORT_IF_POSSIBLE} can be used to indicate that, if possible, 
 * a label should be placed center-aligned with the corresponding port. 
 * 
 * <p><b>This enumeration is not set directly on {@link CoreOptions#PORT_LABELS_PLACEMENT}; instead,
 * an {@code EnumSet} over this enumeration is used there.</b></p>
 * 
 * <p><i>Note:</i> Layout algorithms may only support a subset of these options.</p>
 */
public enum PortLabelPlacement {

    /** Port labels are placed outside of the node. */
    OUTSIDE,

    /** Port labels are placed inside of the node. */
    INSIDE,

    /**
     * Not in all cases port labels are placed <i>next to</i> the port (that is, center-aligned). This option can be set
     * to direct the layout algorithm to place the label next to the port if no edge would cross the label.
     * 
     * Cases in which this behavior is not the default:
     * <ul>
     *   <li>Outside labels of both hierarchical and non-hierarchical nodes,</li>
     *   <li>Inside labels of hierarchical nodes.</li>
     * </ul>
     * 
     * Cases in which the behavior is the default anyway:
     * <ul>
     *   <li>Inside labels of non-hierarchical nodes.</li>
     * </ul>
     */
    NEXT_TO_PORT_IF_POSSIBLE,

    /**
     * The port labels shall always be placed on the same side relative to their corresponding port. For
     * {@link PortSide#WEST} and {@link PortSide#EAST} this is <i>below</i> the port and for {@link PortSide#NORTH} and
     * {@link PortSide#SOUTH} it is <i>right</i> of the port.
     * 
     * <p>
     * Note: the option does not apply to inside port labels (unless a hierarchical edge connects).
     */
    ALWAYS_SAME_SIDE,

    /**
     * The port labels shall always be placed on the same side relative to their corresponding port. For
     * {@link PortSide#WEST} and {@link PortSide#EAST} this is <i>above</i> the port and for {@link PortSide#NORTH} and
     * {@link PortSide#SOUTH} it is <i>left</i> of the port.
     * 
     * <p>
     * Note: the option does not apply to inside port labels (unless a hierarchical edge connects).
     */
    ALWAYS_OTHER_SAME_SIDE,

    /**
     * Unless there are exactly two ports at a given port side, outside port labels are usually all placed to the same
     * side of their port. For example, if there are three northern ports, all of their labels will be placed to the
     * right of their ports. If this option is active, the leftmost label will be placed to the left of its port while
     * the others stay on the right side (and similar for the other port sides). This allows the node to be smaller
     * because the node size doesn't have to accommodate as many port labels, but it breaks symmetry.
     * 
     * <p>
     * Note: the option does not apply to inside port labels (unless a hierarchical edge connects).
     */
    SPACE_EFFICIENT;

    /**
     * @return a set of {@link PortLabelPlacement} enum values indicating fixed port label positions.
     */
    public static EnumSet<PortLabelPlacement> fixed() {
        return EnumSet.noneOf(PortLabelPlacement.class);
    }

    /**
     * @return a set of {@link PortLabelPlacement} enum values indicating default inside port label positions.
     */
    public static EnumSet<PortLabelPlacement> inside() {
        return EnumSet.of(PortLabelPlacement.INSIDE);
    }

    /**
     * @return a set of {@link PortLabelPlacement} enum values indicating default outside port label positions.
     */
    public static EnumSet<PortLabelPlacement> outside() {
        return EnumSet.of(PortLabelPlacement.OUTSIDE);
    }

    /**
     * @return whether the passed placement represents a fixed port label placement, i.e. whether neither
     *         {@link #INSIDE} nor {@link #OUTSIDE} are included.
     */
    public static boolean isFixed(final Set<PortLabelPlacement> placement) {
        return !placement.contains(PortLabelPlacement.INSIDE) && !placement.contains(PortLabelPlacement.OUTSIDE);
    }
    
    /**
     * @return whether the passed {@value placement} is considered valid as described in the class's javadoc.
     */
    public static boolean isValid(final Set<PortLabelPlacement> placement) {
        final Set<PortLabelPlacement> validInsideOutside =
                EnumSet.of(PortLabelPlacement.INSIDE, PortLabelPlacement.OUTSIDE);
        if (Sets.intersection(validInsideOutside, placement).size() > 1) {
            return false;
        }

        final Set<PortLabelPlacement> validPosition =
                EnumSet.of(PortLabelPlacement.ALWAYS_SAME_SIDE,
                        PortLabelPlacement.ALWAYS_OTHER_SAME_SIDE,
                        PortLabelPlacement.SPACE_EFFICIENT);
        if (Sets.intersection(validPosition, placement).size() > 1) {
            return false;
        }

        return true;
    }
}

