/*******************************************************************************
 * Copyright (c) 2015, 2017 Kiel University and others.
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
package org.eclipse.elk.alg.common.nodespacing.internal;

import java.util.EnumSet;
import java.util.List;
import java.util.Set;

import org.eclipse.elk.alg.common.nodespacing.cellsystem.ContainerArea;
import org.eclipse.elk.alg.common.nodespacing.cellsystem.HorizontalLabelAlignment;
import org.eclipse.elk.alg.common.nodespacing.cellsystem.VerticalLabelAlignment;
import org.eclipse.elk.core.options.NodeLabelPlacement;
import org.eclipse.elk.core.options.PortSide;

import com.google.common.collect.Lists;

/**
 * Enumeration over all possible label placements and associated things.
 *
 * @author csp
 * @see NodeLabelPlacement
 */
@SuppressWarnings("unchecked")
public enum NodeLabelLocation {
    /** Outside top left. */
    OUT_T_L(HorizontalLabelAlignment.LEFT,
            VerticalLabelAlignment.BOTTOM,
            ContainerArea.BEGIN,
            ContainerArea.BEGIN,
            EnumSet.of(
                    NodeLabelPlacement.OUTSIDE,
                    NodeLabelPlacement.V_TOP,
                    NodeLabelPlacement.H_LEFT)),
            
    /** Outside top center. */
    OUT_T_C(HorizontalLabelAlignment.CENTER,
            VerticalLabelAlignment.BOTTOM,
            ContainerArea.BEGIN,
            ContainerArea.CENTER,
            EnumSet.of(
                    NodeLabelPlacement.OUTSIDE,
                    NodeLabelPlacement.V_TOP,
                    NodeLabelPlacement.H_CENTER),
            EnumSet.of(
                    NodeLabelPlacement.OUTSIDE,
                    NodeLabelPlacement.V_TOP,
                    NodeLabelPlacement.H_CENTER,
                    NodeLabelPlacement.H_PRIORITY)),
    
    /** Outside top right. */
    OUT_T_R(HorizontalLabelAlignment.RIGHT,
            VerticalLabelAlignment.BOTTOM,
            ContainerArea.BEGIN,
            ContainerArea.END,
            EnumSet.of(
                    NodeLabelPlacement.OUTSIDE,
                    NodeLabelPlacement.V_TOP,
                    NodeLabelPlacement.H_RIGHT)),
    
    /** Outside bottom left. */
    OUT_B_L(HorizontalLabelAlignment.LEFT,
            VerticalLabelAlignment.TOP,
            ContainerArea.END,
            ContainerArea.BEGIN,
            EnumSet.of(
                    NodeLabelPlacement.OUTSIDE,
                    NodeLabelPlacement.V_BOTTOM,
                    NodeLabelPlacement.H_LEFT)),
    
    /** Outside bottom center. */
    OUT_B_C(HorizontalLabelAlignment.CENTER,
            VerticalLabelAlignment.TOP,
            ContainerArea.END,
            ContainerArea.CENTER,
            EnumSet.of(
                    NodeLabelPlacement.OUTSIDE,
                    NodeLabelPlacement.V_BOTTOM,
                    NodeLabelPlacement.H_CENTER),
            EnumSet.of(
                    NodeLabelPlacement.OUTSIDE,
                    NodeLabelPlacement.V_BOTTOM,
                    NodeLabelPlacement.H_CENTER,
                    NodeLabelPlacement.H_PRIORITY)),
    
    /** Outside bottom right. */
    OUT_B_R(HorizontalLabelAlignment.RIGHT,
            VerticalLabelAlignment.TOP,
            ContainerArea.END,
            ContainerArea.END,
            EnumSet.of(
                    NodeLabelPlacement.OUTSIDE,
                    NodeLabelPlacement.V_BOTTOM,
                    NodeLabelPlacement.H_RIGHT)),

    /** Outside left top. */
    OUT_L_T(HorizontalLabelAlignment.RIGHT,
            VerticalLabelAlignment.TOP,
            ContainerArea.BEGIN,
            ContainerArea.BEGIN,
            EnumSet.of(
                    NodeLabelPlacement.OUTSIDE,
                    NodeLabelPlacement.H_LEFT,
                    NodeLabelPlacement.V_TOP,
                    NodeLabelPlacement.H_PRIORITY)),
    
    /** Outside left center. */
    OUT_L_C(HorizontalLabelAlignment.RIGHT,
            VerticalLabelAlignment.CENTER,
            ContainerArea.CENTER,
            ContainerArea.BEGIN,
            EnumSet.of(
                    NodeLabelPlacement.OUTSIDE,
                    NodeLabelPlacement.H_LEFT,
                    NodeLabelPlacement.V_CENTER),
            EnumSet.of(
                    NodeLabelPlacement.OUTSIDE,
                    NodeLabelPlacement.H_LEFT,
                    NodeLabelPlacement.V_CENTER,
                    NodeLabelPlacement.H_PRIORITY)),
    
    /** Outside left bottom. */
    OUT_L_B(HorizontalLabelAlignment.RIGHT,
            VerticalLabelAlignment.BOTTOM,
            ContainerArea.END,
            ContainerArea.BEGIN,
            EnumSet.of(
                    NodeLabelPlacement.OUTSIDE,
                    NodeLabelPlacement.H_LEFT,
                    NodeLabelPlacement.V_BOTTOM,
                    NodeLabelPlacement.H_PRIORITY)),
    
    /** Outside right top. */
    OUT_R_T(HorizontalLabelAlignment.LEFT,
            VerticalLabelAlignment.TOP,
            ContainerArea.BEGIN,
            ContainerArea.END,
            EnumSet.of(
                    NodeLabelPlacement.OUTSIDE,
                    NodeLabelPlacement.H_RIGHT,
                    NodeLabelPlacement.V_TOP,
                    NodeLabelPlacement.H_PRIORITY)),
    
    /** Outside right center. */
    OUT_R_C(HorizontalLabelAlignment.LEFT,
            VerticalLabelAlignment.CENTER,
            ContainerArea.CENTER,
            ContainerArea.END,
            EnumSet.of(
                    NodeLabelPlacement.OUTSIDE,
                    NodeLabelPlacement.H_RIGHT,
                    NodeLabelPlacement.V_CENTER),
            EnumSet.of(
                    NodeLabelPlacement.OUTSIDE,
                    NodeLabelPlacement.H_RIGHT,
                    NodeLabelPlacement.V_CENTER,
                    NodeLabelPlacement.H_PRIORITY)),
    
    /** Outside right bottom. */
    OUT_R_B(HorizontalLabelAlignment.LEFT,
            VerticalLabelAlignment.BOTTOM,
            ContainerArea.END,
            ContainerArea.END,
            EnumSet.of(
                    NodeLabelPlacement.OUTSIDE,
                    NodeLabelPlacement.H_RIGHT,
                    NodeLabelPlacement.V_BOTTOM,
                    NodeLabelPlacement.H_PRIORITY)),
    
    /** Inside top left. */
    IN_T_L(HorizontalLabelAlignment.LEFT,
            VerticalLabelAlignment.TOP,
            ContainerArea.BEGIN,
            ContainerArea.BEGIN,
            EnumSet.of(
                    NodeLabelPlacement.INSIDE,
                    NodeLabelPlacement.V_TOP,
                    NodeLabelPlacement.H_LEFT),
            EnumSet.of(
                    NodeLabelPlacement.INSIDE,
                    NodeLabelPlacement.V_TOP,
                    NodeLabelPlacement.H_LEFT,
                    NodeLabelPlacement.H_PRIORITY)),
    
    /** Inside top center. */
    IN_T_C(HorizontalLabelAlignment.CENTER,
            VerticalLabelAlignment.TOP,
            ContainerArea.BEGIN,
            ContainerArea.CENTER,
            EnumSet.of(
                    NodeLabelPlacement.INSIDE,
                    NodeLabelPlacement.V_TOP,
                    NodeLabelPlacement.H_CENTER),
            EnumSet.of(
                    NodeLabelPlacement.INSIDE,
                    NodeLabelPlacement.V_TOP,
                    NodeLabelPlacement.H_CENTER,
                    NodeLabelPlacement.H_PRIORITY)),
    
    /** Inside top right. */
    IN_T_R(HorizontalLabelAlignment.RIGHT,
            VerticalLabelAlignment.TOP,
            ContainerArea.BEGIN,
            ContainerArea.END,
            EnumSet.of(
                    NodeLabelPlacement.INSIDE,
                    NodeLabelPlacement.V_TOP,
                    NodeLabelPlacement.H_RIGHT),
            EnumSet.of(
                    NodeLabelPlacement.INSIDE,
                    NodeLabelPlacement.V_TOP,
                    NodeLabelPlacement.H_RIGHT,
                    NodeLabelPlacement.H_PRIORITY)),
    
    /** Inside center left. */
    IN_C_L(HorizontalLabelAlignment.LEFT,
            VerticalLabelAlignment.CENTER,
            ContainerArea.CENTER,
            ContainerArea.BEGIN,
            EnumSet.of(
                    NodeLabelPlacement.INSIDE,
                    NodeLabelPlacement.V_CENTER,
                    NodeLabelPlacement.H_LEFT),
            EnumSet.of(
                    NodeLabelPlacement.INSIDE,
                    NodeLabelPlacement.V_CENTER,
                    NodeLabelPlacement.H_LEFT,
                    NodeLabelPlacement.H_PRIORITY)),

    /** Inside center center. */
    IN_C_C(HorizontalLabelAlignment.CENTER,
            VerticalLabelAlignment.CENTER,
            ContainerArea.CENTER,
            ContainerArea.CENTER,
            EnumSet.of(
                    NodeLabelPlacement.INSIDE,
                    NodeLabelPlacement.V_CENTER,
                    NodeLabelPlacement.H_CENTER),
            EnumSet.of(
                    NodeLabelPlacement.INSIDE,
                    NodeLabelPlacement.V_CENTER,
                    NodeLabelPlacement.H_CENTER,
                    NodeLabelPlacement.H_PRIORITY)),

    /** Inside center right. */
    IN_C_R(HorizontalLabelAlignment.RIGHT,
            VerticalLabelAlignment.CENTER,
            ContainerArea.CENTER,
            ContainerArea.END,
            EnumSet.of(
                    NodeLabelPlacement.INSIDE,
                    NodeLabelPlacement.V_CENTER,
                    NodeLabelPlacement.H_RIGHT),
            EnumSet.of(
                    NodeLabelPlacement.INSIDE,
                    NodeLabelPlacement.V_CENTER,
                    NodeLabelPlacement.H_RIGHT,
                    NodeLabelPlacement.H_PRIORITY)),

    /** Inside bottom left. */
    IN_B_L(HorizontalLabelAlignment.LEFT,
            VerticalLabelAlignment.BOTTOM,
            ContainerArea.END,
            ContainerArea.BEGIN,
            EnumSet.of(
                    NodeLabelPlacement.INSIDE,
                    NodeLabelPlacement.V_BOTTOM,
                    NodeLabelPlacement.H_LEFT),
            EnumSet.of(
                    NodeLabelPlacement.INSIDE,
                    NodeLabelPlacement.V_BOTTOM,
                    NodeLabelPlacement.H_LEFT,
                    NodeLabelPlacement.H_PRIORITY)),

    /** Inside bottom center. */
    IN_B_C(HorizontalLabelAlignment.CENTER,
            VerticalLabelAlignment.BOTTOM,
            ContainerArea.END,
            ContainerArea.CENTER,
            EnumSet.of(
                    NodeLabelPlacement.INSIDE,
                    NodeLabelPlacement.V_BOTTOM,
                    NodeLabelPlacement.H_CENTER),
            EnumSet.of(
                    NodeLabelPlacement.INSIDE,
                    NodeLabelPlacement.V_BOTTOM,
                    NodeLabelPlacement.H_CENTER,
                    NodeLabelPlacement.H_PRIORITY)),

    /** Inside bottom right. */
    IN_B_R(HorizontalLabelAlignment.RIGHT,
            VerticalLabelAlignment.BOTTOM,
            ContainerArea.END,
            ContainerArea.END,
            EnumSet.of(
                    NodeLabelPlacement.INSIDE,
                    NodeLabelPlacement.V_BOTTOM,
                    NodeLabelPlacement.H_RIGHT),
            EnumSet.of(
                    NodeLabelPlacement.INSIDE,
                    NodeLabelPlacement.V_BOTTOM,
                    NodeLabelPlacement.H_RIGHT,
                    NodeLabelPlacement.H_PRIORITY)),
    
    /** Undefined or not decidable. */
    UNDEFINED(null, null, null, null);
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Fields

    /** The corresponding placements to this location. */
    private final List<? extends Set<NodeLabelPlacement>> assignedPlacements;
    /** The horizontal text alignment for this location. */
    private final HorizontalLabelAlignment horizontalAlignment;
    /** The vertical text alignment for this location. */
    private final VerticalLabelAlignment verticalAlignment;
    /** The appropriate label row in a container cell. */
    private final ContainerArea containerRow;
    /** The appropriate label column in a container cell. */
    private final ContainerArea containerColumn;
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Construction

    /**
     * Creates a new location with the {@link NodeLabelPlacement}s that lead to this location.
     *
     * @param horizontalAlignment
     *            the horizontal label alignment for this location.
     * @param verticalAlignment
     *            the vertical label alignment for this location.
     * @param row
     *            the appropriate container row in a container cell.
     * @param column
     *            the appropriate container column in a container cell.
     * @param assignedPlacements
     *            the valid {@link NodeLabelPlacement}s for this location.
     */
    NodeLabelLocation(final HorizontalLabelAlignment horizontalAlignment,
            final VerticalLabelAlignment verticalAlignment, final ContainerArea row, final ContainerArea column,
            final Set<NodeLabelPlacement>... assignedPlacements) {
        
        this.horizontalAlignment = horizontalAlignment;
        this.verticalAlignment = verticalAlignment;
        this.containerRow = row;
        this.containerColumn = column;
        this.assignedPlacements = Lists.newArrayList(assignedPlacements);
    }

    /**
     * Converts a set of {@link NodeLabelPlacement}s to a {@link NodeLabelLocation} if possible.
     * 
     * @param labelPlacement
     *            the set of placements to convert.
     * @return the corresponding location. If no valid combination is given,
     *         {@code LabelLocation.UNDEFINED} is returned.
     */
    public static NodeLabelLocation fromNodeLabelPlacement(final Set<NodeLabelPlacement> labelPlacement) {
        for (final NodeLabelLocation location : NodeLabelLocation.values()) {
            if (location.assignedPlacements.contains(labelPlacement)) {
                return location;
            }
        }
        return NodeLabelLocation.UNDEFINED;
    }
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Accessors

    /**
     * Returns the horizontal text alignment for this location.
     * 
     * @return the horizontalAlignment
     */
    public HorizontalLabelAlignment getHorizontalAlignment() {
        return horizontalAlignment;
    }

    /**
     * Returns the vertical text alignment for this location.
     * 
     * @return the verticalAlignment
     */
    public VerticalLabelAlignment getVerticalAlignment() {
        return verticalAlignment;
    }

    /**
     * Returns the appropriate row in a container cell.
     * 
     * @return the row.
     */
    public ContainerArea getContainerRow() {
        return containerRow;
    }

    /**
     * Returns the appropriate column in a container cell.
     * 
     * @return the column.
     */
    public ContainerArea getContainerColumn() {
        return containerColumn;
    }
    
    /**
     * Checks whether this location is inside the node or not.
     */
    public boolean isInsideLocation() {
        switch (this) {
        case IN_T_L:
        case IN_T_C:
        case IN_T_R:
        case IN_C_L:
        case IN_C_C:
        case IN_C_R:
        case IN_B_L:
        case IN_B_C:
        case IN_B_R:
            return true;
        default:
            return false;
        }
    }
    
    /**
     * Returns the side of the node an outside node label location corresponds to or {@link PortSide#UNDEFINED} if this
     * location is on the inside or {@link #UNDEFINED}.
     */
    public PortSide getOutsideSide() {
        switch (this) {
        case OUT_T_L:
        case OUT_T_C:
        case OUT_T_R:
            return PortSide.NORTH;
        
        case OUT_B_L:
        case OUT_B_C:
        case OUT_B_R:
            return PortSide.SOUTH;
            
        case OUT_L_T:
        case OUT_L_C:
        case OUT_L_B:
            return PortSide.WEST;
            
        case OUT_R_T:
        case OUT_R_C:
        case OUT_R_B:
            return PortSide.EAST;
        
        default:
            return PortSide.UNDEFINED;
        }
    }
    
}