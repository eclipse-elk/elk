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
package org.eclipse.elk.alg.common.nodespacing.cellsystem;

import java.util.Arrays;

import org.eclipse.elk.core.math.ElkPadding;
import org.eclipse.elk.core.math.ElkRectangle;

/**
 * A container cell that lays its children out along a strip. The strip can be horizontal (in which case the children
 * can be interpreted as columns) or vertical (in which case the children can be interpreted as rows). A strip
 * container will always make its outer cells the larger of their minimum sizes and left/right- or top/bottom-align
 * them. Whatever remains of its content area after that will be given to the center cell, but the center cell will
 * always be at least its minimum size, even if that causes it to overlap with the outer cells.
 */
public class StripContainerCell extends ContainerCell {
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Properties
    
    /** Whether we lay children out in rows or columns. */
    private final Strip containerMode;
    /** Whether the outer cells should be the same width or height. */
    private final boolean symmetrical;
    /** A container cell can include gaps between its children when calculating its preferred size. */
    private final double gap;
    /** A container cell consists of a number of cells that make up its content. */
    private final Cell[] cells = new Cell[ContainerArea.values().length];
    

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Constructors
    
    /**
     * Creates a new instance with the given settings.
     * 
     * @param mode
     *            whether to lay out children as rows or as columns.
     * @param symmetrical
     *            whether the outer columns should be the same width and the outer rows be the same height.
     * @param gap
     *            the gap inserted between each pair of consecutive cells.
     */
    public StripContainerCell(final Strip mode, final boolean symmetrical, final double gap) {
        this.containerMode = mode;
        this.symmetrical = symmetrical;
        this.gap = gap;
    }
    

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Getters / Setters

    /**
     * Whether this container lays out its children as rows or as columns.
     */
    public Strip getContainerMode() {
        return containerMode;
    }

    /**
     * Returns the gap inserted between each pair of consecutive children.
     */
    public double getGap() {
        return gap;
    }
    
    /**
     * Returns the cell placed at the container's given area.
     * 
     * @param area
     *            the area whose cell to retrieve.
     * @return the cell, which may by {@code null}.
     */
    public Cell getCell(final ContainerArea area) {
        return cells[area.ordinal()];
    }
    
    /**
     * Sets the cell in the container's given area.
     * 
     * @param area
     *            the container's area.
     * @param cell
     *            the new cell to be placed in the area.
     */
    public void setCell(final ContainerArea area, final Cell cell) {
        cells[area.ordinal()] = cell;
    }
    

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // (Container) Cell Methods

    /* (non-Javadoc)
     * @see org.eclipse.elk.core.util.nodespacing.internal.cells.Cell#getMinimumWidth()
     */
    @Override
    public double getMinimumWidth() {
        double width = 0;
        
        if (containerMode == Strip.VERTICAL) {
            // Take the maximum of the child cells
            width = Arrays.stream(cells)
                    .filter(cell -> cell != null && cell.isContributingToMinimumWidth())
                    .mapToDouble(cell -> cell.getMinimumWidth())
                    .max()
                    .orElse(0);
        } else {
            // Minimum widths of the different cells
            double[] cellWidths = minCellWidths(true);
            
            // Keep track of how many cells we have
            int activeCells = 0;
            for (double cellWidth : cellWidths) {
                if (cellWidth > 0) {
                    width += cellWidth;
                    activeCells++;
                }
            }
            
            // If there is more than a single cell, add necessary gaps
            if (activeCells > 1) {
                width += gap * (activeCells - 1);
            }
        }
        
        // If we don't have cells, we don't have width
        return width > 0
                ? width + getPadding().left + getPadding().right
                : 0;
    }
    
    /* (non-Javadoc)
     * @see org.eclipse.elk.core.util.nodespacing.internal.cells.Cell#getMinimumHeight()
     */
    @Override
    public double getMinimumHeight() {
        double height = 0;
        
        if (containerMode == Strip.VERTICAL) {
            // Minimum heights of the different cells
            double[] cellHeights = minCellHeights(true);
            
            // Keep track of how many cells we have
            int activeCells = 0;
            for (double cellHeight : cellHeights) {
                if (cellHeight > 0) {
                    height += cellHeight;
                    activeCells++;
                }
            }
            
            // If there is more than a single cell, add necessary gaps
            if (activeCells > 1) {
                height += gap * (activeCells - 1);
            }
        } else {
            // Take the maximum of the child cells
            height = Arrays.stream(cells)
                    .filter(cell -> cell != null && cell.isContributingToMinimumHeight())
                    .mapToDouble(cell -> cell.getMinimumHeight())
                    .max()
                    .orElse(0);
        }
        
        // If we don't have cells, we don't have height
        return height > 0
                ? height + getPadding().top + getPadding().bottom
                : 0;
    }
    
    /* (non-Javadoc)
     * @see org.eclipse.elk.core.util.nodespacing.internal.cells.ContainerCell#layoutChildrenHorizontally()
     */
    @Override
    public void layoutChildrenHorizontally() {
        ElkRectangle cellRectangle = getCellRectangle();
        ElkPadding cellPadding = getPadding();
        
        if (containerMode == Strip.VERTICAL) {
            // Each child cell begins at our left border (plus padding) and is as large as our content area
            double xPos = cellRectangle.x + cellPadding.left;
            double width = cellRectangle.width - cellPadding.left - cellPadding.right;
            
            for (Cell childCell : cells) {
                applyHorizontalLayout(childCell, xPos, width);
            }
        } else {
            double[] cellWidths = minCellWidths(false);
            
            // Left cell is left-aligned with our content area, right cell is right-aligned
            applyHorizontalLayout(cells[0],
                    cellRectangle.x + cellPadding.left,
                    cellWidths[0]);
            applyHorizontalLayout(cells[2],
                    cellRectangle.x + cellRectangle.width - cellPadding.right - cellWidths[2],
                    cellWidths[2]);
            
            // Size of the content area and size of the available space in the content area
            double freeContentAreaWidth = cellRectangle.width - cellPadding.left - cellPadding.right;
            
            if (cellWidths[0] > 0) {
                freeContentAreaWidth -= cellWidths[0] + gap;
                
                // We add the gap here because that will spare us to check if cellWidths[0] is zero later on
                cellWidths[0] += gap;
            }

            if (cellWidths[2] > 0) {
                freeContentAreaWidth -= cellWidths[2] + gap;
            }
            
            // If the available space is larger than the current size of the center cell, enlarge that thing
            cellWidths[1] = Math.max(cellWidths[1], freeContentAreaWidth);

            // Place the center cell, possibly enlarging it in the process
            applyHorizontalLayout(cells[1],
                    cellRectangle.x + cellPadding.left + cellWidths[0] - (cellWidths[1] - freeContentAreaWidth) / 2,
                    cellWidths[1]);
        }
        
        // Layout container cells recursively
        for (Cell childCell : cells) {
            if (childCell instanceof ContainerCell) {
                ((ContainerCell) childCell).layoutChildrenHorizontally();
            }
        }
    }
    
    /* (non-Javadoc)
     * @see org.eclipse.elk.core.util.nodespacing.internal.cells.ContainerCell#layoutChildrenVertically()
     */
    @Override
    public void layoutChildrenVertically() {
        ElkRectangle cellRectangle = getCellRectangle();
        ElkPadding cellPadding = getPadding();
        
        if (containerMode == Strip.VERTICAL) {
            double[] cellHeights = minCellHeights(false);
            
            // Top cell is top-aligned with our content area, bottom cell is bottom-aligned
            applyVerticalLayout(cells[0],
                    cellRectangle.y + cellPadding.top,
                    cellHeights[0]);
            applyVerticalLayout(cells[2],
                    cellRectangle.y + cellRectangle.height - cellPadding.bottom - cellHeights[2],
                    cellHeights[2]);
            
            // Size of the content area and size of the available space in the content area
            double contentAreaHeight = cellRectangle.height - cellPadding.top - cellPadding.bottom;
            double contentAreaFreeHeight = contentAreaHeight;
            
            if (cellHeights[0] > 0) {
                // We add the gap here because that will spare us to check if cellHeights[0] is zero later on
                cellHeights[0] += gap;
                contentAreaFreeHeight -= cellHeights[0];
            }

            if (cellHeights[2] > 0) {
                contentAreaFreeHeight -= cellHeights[2] + gap;
            }
            
            // If the available space is larger than the current size of the center cell, enlarge that thing
            cellHeights[1] = Math.max(cellHeights[1], contentAreaFreeHeight);

            // Place the center cell, possibly enlarging it in the process
            applyVerticalLayout(cells[1],
                    cellRectangle.y + cellPadding.top + cellHeights[0] - (cellHeights[1] - contentAreaFreeHeight) / 2,
                    cellHeights[1]);
        } else {
            // Each child cell begins at our top border (plus padding) and is as large as our content area
            double yPos = cellRectangle.y + cellPadding.top;
            double height = cellRectangle.height - cellPadding.top - cellPadding.bottom;
            
            for (Cell childCell : cells) {
                applyVerticalLayout(childCell, yPos, height);
            }
        }
        
        // Layout container cells recursively
        for (Cell childCell : cells) {
            if (childCell instanceof ContainerCell) {
                ((ContainerCell) childCell).layoutChildrenVertically();
            }
        }
    }
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Utilities
    
    /**
     * Returns an array containing the width of each cell, with symmetry applied if activated.
     */
    private double[] minCellWidths(final boolean respectContributionFlag) {
        double[] cellWidths = {
                minWidthOfCell(cells[0], respectContributionFlag),
                minWidthOfCell(cells[1], respectContributionFlag),
                minWidthOfCell(cells[2], respectContributionFlag)
        };
        
        // If we are to be symmetrical, the outer cells need to be the same size
        if (symmetrical) {
            cellWidths[0] = Math.max(cellWidths[0], cellWidths[2]);
            cellWidths[2] = cellWidths[0];
        }
        
        return cellWidths;
    }
    
    /**
     * Returns an array containing the height of each cell, with symmetry applied if activated.
     */
    private double[] minCellHeights(final boolean respectContributionFlag) {
        double[] cellHeights = {
                minHeightOfCell(cells[0], respectContributionFlag),
                minHeightOfCell(cells[1], respectContributionFlag),
                minHeightOfCell(cells[2], respectContributionFlag)
        };
        
        // If we are to be symmetrical, the outer cells need to be the same size
        if (symmetrical) {
            cellHeights[0] = Math.max(cellHeights[0], cellHeights[2]);
            cellHeights[2] = cellHeights[0];
        }
        
        return cellHeights;
    }

    /**
     * The mode of a {@link StripContainerCell} describes whether its child cells are to be regarded as rows or as
     * columns.
     */
    public static enum Strip {
        /** In a vertical strip, the container's child cells are its rows. */
        VERTICAL,
        /** In a horizontal strip, the container's child cells are its columns. */
        HORIZONTAL;
    }
    
}
