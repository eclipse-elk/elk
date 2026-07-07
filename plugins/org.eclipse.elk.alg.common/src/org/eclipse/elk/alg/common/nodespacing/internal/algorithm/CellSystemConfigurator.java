/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
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
package org.eclipse.elk.alg.common.nodespacing.internal.algorithm;

import org.eclipse.elk.alg.common.nodespacing.cellsystem.AtomicCell;
import org.eclipse.elk.alg.common.nodespacing.cellsystem.LabelCell;
import org.eclipse.elk.alg.common.nodespacing.internal.NodeContext;
import org.eclipse.elk.alg.common.nodespacing.internal.NodeLabelLocation;
import org.eclipse.elk.core.options.PortConstraints;
import org.eclipse.elk.core.options.PortSide;
import org.eclipse.elk.core.options.SizeConstraint;
import org.eclipse.elk.core.options.SizeOptions;

/**
 * Configures constraints of the cell system such that the various cells contribute properly to the node size
 * calculation.
 */
public final class CellSystemConfigurator {

    /**
     * No instance required.
     */
    private CellSystemConfigurator() {
        
    }
    

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Size Contribution Configuration
    
    /**
     * Configures the cell system's constraints such that they work properly when calculating the required node space.
     */
    public static void configureCellSystemSizeContributions(final NodeContext nodeContext) {
        // If the node has a fixed size, we don't need to change anything because the cell system won't be used to
        // calculate the node's size
        if (nodeContext.sizeConstraints.isEmpty()) {
            return;
        }
        
        // Go through the different size constraint components
        if (nodeContext.sizeConstraints.contains(SizeConstraint.PORTS)) {
            // The northern and southern inside port label cells have the correct width for the node
            nodeContext.insidePortLabelCells.get(PortSide.NORTH).setContributesToMinimumWidth(true);
            nodeContext.insidePortLabelCells.get(PortSide.SOUTH).setContributesToMinimumWidth(true);
            
            // For the eastern and western cells, they only give a correct height if port placement is free instead of
            // constrained (the constrained case is handled separately in the node size calculator). This is due to the
            // fact that in the cell layout, the north and south cells are above and below the east and west cells,
            // which would cause the node to become too high
            boolean freePortPlacement = nodeContext.portConstraints != PortConstraints.FIXED_RATIO
                    && nodeContext.portConstraints != PortConstraints.FIXED_POS;
            
            nodeContext.insidePortLabelCells.get(PortSide.EAST).setContributesToMinimumHeight(freePortPlacement);
            nodeContext.insidePortLabelCells.get(PortSide.WEST).setContributesToMinimumHeight(freePortPlacement);
            
            // The main row needs to contribute height for the east and west port label cells to be able to contribute
            // their height
            nodeContext.nodeContainerMiddleRow.setContributesToMinimumHeight(freePortPlacement);
            
            // Port labels only contribute their size if ports are accounted for as well
            if (nodeContext.sizeConstraints.contains(SizeConstraint.PORT_LABELS)) {
                // The port label cells contribute the space they need for inside port label placement
                nodeContext.insidePortLabelCells.get(PortSide.NORTH).setContributesToMinimumHeight(true);
                nodeContext.insidePortLabelCells.get(PortSide.SOUTH).setContributesToMinimumHeight(true);
                nodeContext.insidePortLabelCells.get(PortSide.EAST).setContributesToMinimumWidth(true);
                nodeContext.insidePortLabelCells.get(PortSide.WEST).setContributesToMinimumWidth(true);
                
                // The main row needs to contribute Width for the east and west port label cells to be able to
                // contribute their width
                nodeContext.nodeContainerMiddleRow.setContributesToMinimumWidth(true);
            }
        }
        
        if (nodeContext.sizeConstraints.contains(SizeConstraint.NODE_LABELS)) {
            // The inside node label cell needs to contribute both width and height, as needs the middle row
            nodeContext.insideNodeLabelContainer.setContributesToMinimumHeight(true);
            nodeContext.insideNodeLabelContainer.setContributesToMinimumWidth(true);
            
            nodeContext.nodeContainerMiddleRow.setContributesToMinimumHeight(true);
            nodeContext.nodeContainerMiddleRow.setContributesToMinimumWidth(true);
            
            // All node label cells need to contribute height and width, but outside node labels only do so unless they
            // are configured to overhang
            boolean overhang = nodeContext.sizeOptions.contains(SizeOptions.OUTSIDE_NODE_LABELS_OVERHANG);
            for (NodeLabelLocation location : NodeLabelLocation.values()) {
                LabelCell labelCell = nodeContext.nodeLabelCells.get(location);
                if (labelCell != null) {
                    if (location.isInsideLocation()) {
                        labelCell.setContributesToMinimumHeight(true);
                        labelCell.setContributesToMinimumWidth(true);
                    } else {
                        labelCell.setContributesToMinimumHeight(!overhang);
                        labelCell.setContributesToMinimumWidth(!overhang);
                    }
                }
            }
        }
        
        // If the middle cell contributes to the node size, we need to set that up as well
        if (nodeContext.sizeConstraints.contains(SizeConstraint.MINIMUM_SIZE)
                && nodeContext.sizeOptions.contains(SizeOptions.MINIMUM_SIZE_ACCOUNTS_FOR_PADDING)) {
            
            // The middle row now needs to contribute width and height, and the center cell of the inside node label
            // container needs to contribute width and height as well
            nodeContext.nodeContainerMiddleRow.setContributesToMinimumHeight(true);
            nodeContext.nodeContainerMiddleRow.setContributesToMinimumHeight(true);
            
            // If the inside node label container is not already contributing to the minimum height and width, node
            // labels are not to be regarded. In that case, turn size contributions on, but limit them to the node
            // label container's center cell
            if (!nodeContext.insideNodeLabelContainer.isContributingToMinimumHeight()) {
                nodeContext.insideNodeLabelContainer.setContributesToMinimumHeight(true);
                nodeContext.insideNodeLabelContainer.setContributesToMinimumWidth(true);
                nodeContext.insideNodeLabelContainer.setOnlyCenterCellContributesToMinimumSize(true);
            }
        }
    }
    

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Update East and West Inside Port Label Cells
    
    /**
     * The padding of east and west inside port label cells was originally set to the surrounding port margins. If the
     * port constraints were fixed, the paddings were reset accordingly. In the free case, the paddings must be such
     * that the port placement will start below the northern and southern inside port label space, but still respects
     * the surrounding port margins.
     */
    public static void updateVerticalInsidePortLabelCellPadding(final NodeContext nodeContext) {
        // We only care for the free port placement case
        if (nodeContext.portConstraints == PortConstraints.FIXED_RATIO
                || nodeContext.portConstraints == PortConstraints.FIXED_POS) {
            
            return;
        }
        
        // Calculate where the east and west port cells will end up
        double topBorderOffset = nodeContext.nodeContainer.getPadding().top
                + nodeContext.insidePortLabelCells.get(PortSide.NORTH).getMinimumHeight()
                + nodeContext.labelCellSpacing;
        double bottomBorderOffset = nodeContext.nodeContainer.getPadding().bottom
                + nodeContext.insidePortLabelCells.get(PortSide.SOUTH).getMinimumHeight()
                + nodeContext.labelCellSpacing;
        
        AtomicCell eastCell = nodeContext.insidePortLabelCells.get(PortSide.EAST);
        AtomicCell westCell = nodeContext.insidePortLabelCells.get(PortSide.WEST);
        
        // Calculate how much top padding we actually need
        double topPadding = Math.max(0, eastCell.getPadding().top - topBorderOffset);
        topPadding = Math.max(topPadding, westCell.getPadding().top - topBorderOffset);
        double bottomPadding = Math.max(0, eastCell.getPadding().bottom - bottomBorderOffset);
        bottomPadding = Math.max(bottomPadding, westCell.getPadding().bottom - bottomBorderOffset);

        // Update paddings
        eastCell.getPadding().top = topPadding;
        westCell.getPadding().top = topPadding;
        eastCell.getPadding().bottom = bottomPadding;
        westCell.getPadding().bottom = bottomPadding;
    }
    
}
