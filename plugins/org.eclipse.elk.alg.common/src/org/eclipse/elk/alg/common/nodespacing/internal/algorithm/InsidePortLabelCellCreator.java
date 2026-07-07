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
import org.eclipse.elk.alg.common.nodespacing.cellsystem.ContainerArea;
import org.eclipse.elk.alg.common.nodespacing.cellsystem.StripContainerCell;
import org.eclipse.elk.alg.common.nodespacing.internal.NodeContext;
import org.eclipse.elk.alg.common.nodespacing.internal.PortContext;
import org.eclipse.elk.core.math.ElkPadding;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.options.PortLabelPlacement;
import org.eclipse.elk.core.options.PortSide;

/**
 * Sets up the inside port label cells. These are set up even when there are no inside port labels since they also
 * determine how much space we need to place ports along the node borders.
 */
public final class InsidePortLabelCellCreator {
    
    /**
     * No instance required.
     */
    private InsidePortLabelCellCreator() {
        
    }
    

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Inside Port Label Cells
    
    /**
     * Creates all the inside port label cells. Setting them up is left to the port and port label placement code.
     */
    public static void createInsidePortLabelCells(final NodeContext nodeContext) {
        // Create all inside port label cells
        createInsidePortLabelCell(nodeContext, nodeContext.nodeContainer, ContainerArea.BEGIN, PortSide.NORTH);
        createInsidePortLabelCell(nodeContext, nodeContext.nodeContainer, ContainerArea.END, PortSide.SOUTH);

        createInsidePortLabelCell(nodeContext, nodeContext.nodeContainerMiddleRow, ContainerArea.BEGIN, PortSide.WEST);
        createInsidePortLabelCell(nodeContext, nodeContext.nodeContainerMiddleRow, ContainerArea.END, PortSide.EAST);
        
        setupNorthOrSouthPortLabelCell(nodeContext, PortSide.NORTH);
        setupNorthOrSouthPortLabelCell(nodeContext, PortSide.SOUTH);
        setupEastOrWestPortLabelCell(nodeContext, PortSide.EAST);
        setupEastOrWestPortLabelCell(nodeContext, PortSide.WEST);
    }
    
    /**
     * Creates a new {@link AtomicCell} for inside port labels. The cell is placed in the given area of the given
     * container cell. Also, the cell is registered to be responsible for ports of the given side.
     */
    private static void createInsidePortLabelCell(final NodeContext nodeContext, final StripContainerCell container,
            final ContainerArea containerArea, final PortSide portSide) {
        
        AtomicCell portLabelCell = new AtomicCell();
        container.setCell(containerArea, portLabelCell);
        nodeContext.insidePortLabelCells.put(portSide, portLabelCell);
    }
    

    /////////////////////////////////////////////////////////////////
    // North or South
    
    /**
     * Sets up all paddings of this cell. We cannot calculate a height yet because that will only become known once
     * inside node label placement has finished.
     */
    private static void setupNorthOrSouthPortLabelCell(final NodeContext nodeContext, final PortSide portSide) {
        ElkPadding padding = nodeContext.insidePortLabelCells.get(portSide).getPadding();
        
        switch (portSide) {
        case NORTH:
            // In case of negative port spacing, do not use it as a padding since this would increase the node size.
            // Needed to ensure that negative label port spacing behaves the same as positive.
            if (nodeContext.portLabelSpacingVertical >= 0) {
                padding.top = nodeContext.portLabelSpacingVertical;
            }
            break;
            
        case SOUTH:
            // In case of negative port spacing, do not use it as a padding since this would increase the node size.
            // Needed to ensure that negative label port spacing behaves the same as positive.
            if (nodeContext.portLabelSpacingVertical >= 0) {
                padding.bottom = nodeContext.portLabelSpacingVertical;
            }
            break;
        }
        
        if (nodeContext.surroundingPortMargins != null) {
            padding.left = nodeContext.surroundingPortMargins.left;
            padding.right = nodeContext.surroundingPortMargins.right;
        }
    }
    

    /////////////////////////////////////////////////////////////////
    // East or West
    
    /**
     * Does the work for {@link #calculatePortLabelWidth(NodeContext)} for the given port side.
     */
    private static void setupEastOrWestPortLabelCell(final NodeContext nodeContext, final PortSide portSide) {
        if (nodeContext.portLabelsPlacement.contains(PortLabelPlacement.INSIDE)) {
            calculateWidthDueToLabels(nodeContext, portSide);   
        }
        setupTopAndBottomPadding(nodeContext, portSide);
    }

    /**
     * Calculates the cell's width to fit the largest label in there. If there actually is a label, also setup the
     * cell's top or bottom padding to ensure enough space between ports and labels.
     */
    private static void calculateWidthDueToLabels(final NodeContext nodeContext, final PortSide portSide) {
        // Retrieve the appropriate cell
        AtomicCell theAppropriateCell = nodeContext.insidePortLabelCells.get(portSide);
        KVector minCellSize = theAppropriateCell.getMinimumContentAreaSize();
        
        for (PortContext portContext : nodeContext.portContexts.get(portSide)) {
            // Update the maximum label width
            if (portContext.portLabelCell != null) {
                minCellSize.x = Math.max(minCellSize.x, portContext.portLabelCell.getMinimumWidth());
            }
        }
        
        // If the cell has a minimum width by now, that means we actually have labels in there. Which, in turn, means
        // that we need to add a padding to the cell to ensure enough space between ports and their inside labels
        if (minCellSize.x > 0) {
            switch (portSide) {
            case EAST:
                theAppropriateCell.getPadding().right = nodeContext.portLabelSpacingHorizontal;
                break;
                
            case WEST:
                theAppropriateCell.getPadding().left = nodeContext.portLabelSpacingHorizontal;
                break;
            }
        }
    }

    /**
     * Sets up the cell's top and bottom padding to match the surrounding port margins.
     */
    private static void setupTopAndBottomPadding(final NodeContext nodeContext, final PortSide portSide) {
        if (nodeContext.surroundingPortMargins != null) {
            ElkPadding padding = nodeContext.insidePortLabelCells.get(portSide).getPadding();
            padding.top = nodeContext.surroundingPortMargins.top;
            padding.bottom = nodeContext.surroundingPortMargins.bottom;
        }
    }
    
}
