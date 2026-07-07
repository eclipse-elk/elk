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
import java.util.Set;

import com.google.common.collect.Sets;

/**
 * Options for controlling how node labels are placed by layout algorithms. The corresponding layout
 * option will usually accept an {@link EnumSet} over this enumeration, theoretically allowing
 * arbitrary and even contradictory subsets of options to be set. Note that you are restricted to
 * the following combinations if you want your choice to make sense:
 * <ul>
 *   <li>Exactly one of the {@link #INSIDE} and {@link #OUTSIDE} options.</li>
 *   <li>Exactly one of the {@link #H_LEFT}, {@link #H_CENTER}, and {@link #H_RIGHT} options.</li>
 *   <li>Exactly one of the {@link #V_TOP}, {@link #V_CENTER}, and {@link #V_BOTTOM} options.</li>
 * </ul>
 * <p><b>This enumeration is not set directly on {@link CoreOptions#NODE_LABEL_PLACEMENT}; instead,
 * an {@code EnumSet} over this enumeration is used there.</b></p>
 * 
 * <p><i>Note:</i> Layout algorithms may only support a subset of these options.</p>
 *
 * @author msp
 * @author cds
 */
public enum NodeLabelPlacement {
    
    /**
     * Horizontal left placement.
     */
    H_LEFT,
    
    /**
     * Horizontal center placement.
     */
    H_CENTER,
    
    /**
     * Horizontal right placement.
     */
    H_RIGHT,
    
    /**
     * Vertical top placement.
     */
    V_TOP,
    
    /**
     * Vertical center placement.
     */
    V_CENTER,
    
    /**
     * Vertical bottom placement.
     */
    V_BOTTOM,
    
    /**
     * Place node labels on the inside of the node. This should usually be combined with size
     * constraints to ensure the node is big enough to accomodate its labels.
     */
    INSIDE,
    
    /**
     * Place node labels on the outside of the node.
     */
    OUTSIDE,
    
    /**
     * If set, the default behaviour is changed to give the horizontal placement priority over the
     * vertical one.
     * 
     * <h2>Inside Placement</h2>
     * 
     * <p>The default behaviour is for top / bottom labels to contribute to the top / bottom insets
     * of a node, thus causing the node's children to be moved down / leave space to the node's bottom.</p>
     * 
     * <p>If this option is set, the default behaviour is overridden and the horizontal placement options
     * are evaluated first. This makes a left top label occupy space on the left side of the node, causing
     * the node's children to be moved right instead of down.</p>
     * 
     * <p><em>Note:</em> Since this changes the placement of a node's children, this option is only
     * interpreted if set on the node itself, not on the individual labels.
     * 
     * 
     * <h2>Outside Placement</h2>
     * 
     * <p>The default behaviour is to first choose the northern or southern node side for placement of
     * the label, according to the {@link #V_TOP} and {@link #V_BOTTOM} options. The horizontal placement
     * is then constrained to the width of the node, thus resulting in labels that are placed above or
     * below the node, but never to its left or right side.</p>
     * 
     * <p>If this option is set, the default behaviour is overridden and the horizontal placement options
     * are evaluated first, thus causing the eastern or western node side to be chosen for placement of
     * the label. The vertical placement options are then constrained to the height of the node, thus
     * resulting in labels that are placed to the left or right of the node.</p>
     */
    H_PRIORITY;
    
    
    /**
     * Returns the enumeration value related to the given ordinal.
     * 
     * @param i ordinal value
     * @return the related enumeration value
     */
    public static NodeLabelPlacement valueOf(final int i) {
        return values()[i];
    }
    

    /**
     * Returns an empty enum set over this enumeration, which prevents the layout algorithm from
     * changing the label's coordinates.
     * 
     * @return set over this enumeration representing fixed node label placement constraints.
     */
    public static EnumSet<NodeLabelPlacement> fixed() {
        return EnumSet.noneOf(NodeLabelPlacement.class);
    }
    
    /**
     * Returns a node label placement to place the node label inside the node, left-aligned on top.
     * 
     * @return node label placement for inside top left placement.
     */
    public static EnumSet<NodeLabelPlacement> insideTopLeft() {
        return EnumSet.of(INSIDE, V_TOP, H_LEFT);
    }
    
    /**
     * Returns a node label placement to place the node label inside the node, centered at its top.
     * 
     * @return node label placement for inside top center placement.
     */
    public static EnumSet<NodeLabelPlacement> insideTopCenter() {
        return EnumSet.of(INSIDE, V_TOP, H_CENTER);
    }
    
    /**
     * Returns a node label placement to place the node label inside the node, right-aligned on top.
     * 
     * @return node label placement for inside top right placement.
     */
    public static EnumSet<NodeLabelPlacement> insideTopRight() {
        return EnumSet.of(INSIDE, V_TOP, H_RIGHT);
    }
    
    /**
     * Returns a node label placement to place the node label centered inside the node.
     * 
     * @return node label placement for inside centered placement.
     */
    public static EnumSet<NodeLabelPlacement> insideCenter() {
        return EnumSet.of(INSIDE, V_CENTER, H_CENTER);
    }
    
    /**
     * Returns a node label placement to place the node label inside the node, left-aligned on bottom.
     * 
     * @return node label placement for inside top bottom placement.
     */
    public static EnumSet<NodeLabelPlacement> insideBottomLeft() {
        return EnumSet.of(INSIDE, V_BOTTOM, H_LEFT);
    }
    
    /**
     * Returns a node label placement to place the node label inside the node, centered at its bottom.
     * 
     * @return node label placement for inside bottom center placement.
     */
    public static EnumSet<NodeLabelPlacement> insideBottomCenter() {
        return EnumSet.of(INSIDE, V_BOTTOM, H_CENTER);
    }
    
    /**
     * Returns a node label placement to place the node label inside the node, right-aligned on bottom.
     * 
     * @return node label placement for inside top right placement.
     */
    public static EnumSet<NodeLabelPlacement> insideBottomRight() {
        return EnumSet.of(INSIDE, V_BOTTOM, H_RIGHT);
    }
    
    /**
     * Returns a node label placement to place the node label outside the node, left-aligned on top.
     * 
     * @return node label placement for outside top left placement.
     */
    public static EnumSet<NodeLabelPlacement> outsideTopLeft() {
        return EnumSet.of(OUTSIDE, V_TOP, H_LEFT);
    }
    
    /**
     * Returns a node label placement to place the node label outside the node, centered on top.
     * 
     * @return node label placement for outside top center placement.
     */
    public static EnumSet<NodeLabelPlacement> outsideTopCenter() {
        return EnumSet.of(OUTSIDE, V_TOP, H_CENTER);
    }
    
    /**
     * Returns a node label placement to place the node label outside the node, right-aligned on top.
     * 
     * @return node label placement for outside top left placement.
     */
    public static EnumSet<NodeLabelPlacement> outsideTopRight() {
        return EnumSet.of(OUTSIDE, V_TOP, H_RIGHT);
    }
    
    /**
     * Returns a node label placement to place the node label outside the node, left-aligned on bottom.
     * 
     * @return node label placement for outside top left placement.
     */
    public static EnumSet<NodeLabelPlacement> outsideBottomLeft() {
        return EnumSet.of(OUTSIDE, V_BOTTOM, H_LEFT);
    }

    /**
     * Returns a node label placement to place the node label outside the node, centered below it.
     * 
     * @return node label placement for outside bottom center placement.
     */
    public static EnumSet<NodeLabelPlacement> outsideBottomCenter() {
        return EnumSet.of(OUTSIDE, V_BOTTOM, H_CENTER);
    }
    
    /**
     * Returns a node label placement to place the node label outside the node, right-aligned on bottom.
     * 
     * @return node label placement for outside top left placement.
     */
    public static EnumSet<NodeLabelPlacement> outsideBottomRight() {
        return EnumSet.of(OUTSIDE, V_BOTTOM, H_RIGHT);
    }
    
    
    /**
     * @return whether the passed {@value placement} is considered valid as described in the class's javadoc.
     */
    public static boolean isValid(final Set<NodeLabelPlacement> placement) {
        final Set<NodeLabelPlacement> validInsideOutside =
                EnumSet.of(NodeLabelPlacement.INSIDE, NodeLabelPlacement.OUTSIDE);
        if (Sets.intersection(validInsideOutside, placement).size() > 1) {
            return false;
        }

        final Set<NodeLabelPlacement> validHorizontal =
                EnumSet.of(NodeLabelPlacement.H_LEFT, NodeLabelPlacement.H_CENTER, NodeLabelPlacement.H_RIGHT);
        if (Sets.intersection(validHorizontal, placement).size() > 1) {
            return false;
        }

        final Set<NodeLabelPlacement> validVertical =
                EnumSet.of(NodeLabelPlacement.V_TOP, NodeLabelPlacement.V_CENTER, NodeLabelPlacement.V_BOTTOM);
        if (Sets.intersection(validVertical, placement).size() > 1) {
            return false;
        }

        return true;
    }
}
