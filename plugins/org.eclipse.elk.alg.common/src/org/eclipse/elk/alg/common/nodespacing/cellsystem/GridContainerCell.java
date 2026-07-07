/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.common.nodespacing.cellsystem;

import org.eclipse.elk.core.math.ElkPadding;
import org.eclipse.elk.core.math.ElkRectangle;
import org.eclipse.elk.core.math.KVector;

/**
 * A container that lays out its child cells in three rows. It's basically not much more than an extension of
 * {@link StripContainerCell} into another dimension, which is as awesome as it sounds. Anyway, cells can be added to
 * the container as usual. The grid container can operate in two modes: tabular and non-tabular. Non-tabular mode works
 * similar to three {@link StripContainerCell} instances in that the columns of the three rows don't influence each
 * other in terms of their size. Tabular mode locks the columns of all rows, making them end up with the same widths.
 * 
 * <p>One special feature is that its center cell can be assigned a custom minimal size which will
 * be used regardless of the actual center cell's minimal size. This grid container can even be configured to only use
 * that custom size as its minimal size. As a final goodie, the grid allows you to retrieve the actual rectangle in the
 * center that is not overlapped by any of the outer grid cells. This may not be equal to the cell rectangle computed
 * for the center cell, which always has at least the cell's minimum size.</p>
 */
public class GridContainerCell extends ContainerCell {
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Constants
    
    /** Number of rows in the grid. */
    private static final int ROWS = ContainerArea.values().length;
    /** Number of columns in the grid. */
    private static final int COLUMNS = ROWS;
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Properties
    
    /** A container cell can include gaps between its children when calculating its preferred size. */
    private final double gap;
    /** Whether the container works in tabular mode or not. */
    private final boolean tabular;
    /** Whether the outer columns should be the same width and the outer rows be the same height. */
    private final boolean symmetrical;
    /** A container cell consists of a table of cells that make up its content. */
    private final Cell[][] cells = new Cell[ROWS][COLUMNS];
    /**
     * The center cell of the grid may have a custom minimum size (used to empose a mininum size on the node's
     * client area.
     */
    private KVector centerCellMinimumSize = null;
    /**
     * Whether only the center cell contributes to the minimum size of this cell. If this is true, we actually use
     * the custom center cell minimum size instead of the center cell's actual minimum size, if set.
     */
    private boolean onlyCenterCellContributesToMinimumSize = false;
    /** Rectangle that describes the part of the grid free of outer grid cells. May be empty. */
    private ElkRectangle centerCellRect = new ElkRectangle();
    

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Constructors
    
    /**
     * Creates a new instance with the given settings.
     * 
     * @param tabular
     *            whether the width of columns should be locked across the three rows.
     * @param symmetrical
     *            whether the outer columns should be the same width and the outer rows be the same height.
     * @param gap
     *            the gap inserted between each pair of consecutive cells.
     */
    public GridContainerCell(final boolean tabular, final boolean symmetrical, final double gap) {
        this.tabular = tabular;
        this.symmetrical = symmetrical;
        this.gap = gap;
    }
    

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Getters / Setters

    /**
     * Returns the gap inserted between each pair of consecutive children.
     */
    public double getGap() {
        return gap;
    }
    
    /**
     * Returns the cell placed at the container's given area.
     * 
     * @param row
     *            the row.
     * @param col
     *            the column.
     * @return the cell, which may by {@code null}.
     */
    public Cell getCell(final ContainerArea row, final ContainerArea col) {
        return cells[row.ordinal()][col.ordinal()];
    }
    
    /**
     * Sets the cell in the container's given area.
     * 
     * @param row
     *            the row.
     * @param col
     *            the column.
     * @param cell
     *            the new cell to be placed in the area.
     */
    public void setCell(final ContainerArea row, final ContainerArea col, final Cell cell) {
        cells[row.ordinal()][col.ordinal()] = cell;
    }
    
    /**
     * Sets the center cell's minimum size. Regardless of the actuall cell in the grid's center, this minimum size
     * imposes a lower bound on the center cell size.
     */
    public void setCenterCellMinimumSize(final KVector minimumSize) {
        this.centerCellMinimumSize = new KVector(minimumSize);
    }
    
    /**
     * Sets whether when this container is queried for its minimum size, only the center cell will be included in the
     * calculations. If the center cell has a minimum size set and this mode is active, that minimum size will be used
     * regardless of the actuall center cell.
     */
    public void setOnlyCenterCellContributesToMinimumSize(final boolean contribution) {
        this.onlyCenterCellContributesToMinimumSize = contribution;
    }
    
    /**
     * Returns a rectangle that describes the space left free by the outer rows and columns. May well be an empty
     * rectangle. The rectangle's information have only been computed once both {@link #layoutChildrenHorizontally()}
     * and {@link #layoutChildrenVertically()} have been called.
     */
    public ElkRectangle getCenterCellRectangle() {
        return new ElkRectangle(centerCellRect);
    }
    

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // (Container) Cell Methods

    /* (non-Javadoc)
     * @see org.eclipse.elk.core.util.nodespacing.internal.cells.Cell#getMinimumWidth()
     */
    @Override
    public double getMinimumWidth() {
        double width = 0;
        
        if (onlyCenterCellContributesToMinimumSize) {
            // Check if an explicit minimum size for the center cell was set
            if (centerCellMinimumSize != null) {
                width = centerCellMinimumSize.x;
            } else if (cells[1][1] != null) {
                width = cells[1][1].getMinimumWidth();
            }
            
        } else if (tabular) {
            // Use aggregated widths
            width = sumWithGaps(minColumnWidths(null, true));
            
        } else {
            // Use maximum width over each row
            for (ContainerArea area : ContainerArea.values()) {
                width = Math.max(width,
                        sumWithGaps(minColumnWidths(area, true)));
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
        
        if (onlyCenterCellContributesToMinimumSize) {
            // Check if an explicit minimum size for the center cell was set
            if (centerCellMinimumSize != null) {
                height = centerCellMinimumSize.y;
            } else if (cells[1][1] != null) {
                height = cells[1][1].getMinimumHeight();
            }
            
        } else {
            // Minimum height of the different rows (this is independent of whether we are in tabular mode or not)
            height = sumWithGaps(minRowHeights(true));
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
        // How we're going to do this depends on whether we're in tabular mode or not. If so, the column widths across
        // all rows are locked
        if (tabular) {
            double[] colWidths = minColumnWidths(null, false);
            for (ContainerArea row : ContainerArea.values()) {
                applyWidthsToRow(row, colWidths);
            }
            
        } else {
            for (ContainerArea row : ContainerArea.values()) {
                double[] colWidths = minColumnWidths(row, false);
                applyWidthsToRow(row, colWidths);
            }
        }
    }
    
    /* (non-Javadoc)
     * @see org.eclipse.elk.core.util.nodespacing.internal.cells.ContainerCell#layoutChildrenVertically()
     */
    @Override
    public void layoutChildrenVertically() {
        // Note that the code in this method does not depend on whether we're in tabular mode or not, that's why it
        // looks significantly different from its horizontal counterpart
        ElkRectangle cellRectangle = getCellRectangle();
        ElkPadding cellPadding = getPadding();
        
        double[] rowHeights = minRowHeights(false);
        
        // Top row is top-aligned with our content area, bottom row is bottom-aligned
        applyHeightToRow(ContainerArea.BEGIN,
                cellRectangle.y + cellPadding.top,
                rowHeights);
        applyHeightToRow(ContainerArea.END,
                cellRectangle.y + cellRectangle.height - cellPadding.bottom - rowHeights[2],
                rowHeights);
        
        // Size of the content area and size of the available space in the content area
        double freeContentAreaHeight = cellRectangle.height - cellPadding.top - cellPadding.bottom;
        
        if (rowHeights[0] > 0) {
            rowHeights[0] += gap;
            freeContentAreaHeight -= rowHeights[0];
        }
        
        if (rowHeights[2] > 0) {
            rowHeights[2] += gap;
            freeContentAreaHeight -= rowHeights[2];
        }
        
        // Compute the center cell rectangle
        centerCellRect.height = Math.max(0, freeContentAreaHeight);
        centerCellRect.y = cellRectangle.y + cellPadding.top + (centerCellRect.height - freeContentAreaHeight) / 2;
        
        // If the available space is larger than the current size of the center cell, enlarge that thing
        rowHeights[1] = Math.max(rowHeights[1], freeContentAreaHeight);

        // Place the center cell, possibly enlarging it in the process
        applyHeightToRow(ContainerArea.CENTER,
                cellRectangle.y + cellPadding.top + rowHeights[0] - (rowHeights[1] - freeContentAreaHeight) / 2,
                rowHeights);
    }
    

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Width and Height Calculations
    
    /**
     * Returns an array containing the width of each column for either the whole table or only a specific row. In the
     * former case, the maximum width of the respective column in each row is used. Symmetry is applied if activated.
     * 
     * @param row
     *            only return the column widths of the given row. Return the aggregated width for the whole table if
     *            {@code row == null}.
     */
    private double[] minColumnWidths(final ContainerArea row, final boolean respectContributionFlag) {
        double[] colWidths = {
                minWidthOfColumn(ContainerArea.BEGIN, row, respectContributionFlag),
                minWidthOfColumn(ContainerArea.CENTER, row, respectContributionFlag),
                minWidthOfColumn(ContainerArea.END, row, respectContributionFlag)
        };
        
        // If we are to be symmetrical, the outer cells need to be the same size
        if (symmetrical) {
            colWidths[0] = Math.max(colWidths[0], colWidths[2]);
            colWidths[2] = colWidths[0];
        }
        
        return colWidths;
    }
    
    /**
     * Returns the minimum width of the given column. If {@code row != null}, the maximum of the minimum width of that
     * column across all rows is returned. Otherwise, we concentrate only on the column of the given row.
     */
    private double minWidthOfColumn(final ContainerArea column, final ContainerArea row,
            final boolean respectContributionFlag) {
        
        double maxMinWidth = 0;
        
        if (row == null) {
            // Aggregate values for all rows
            for (int rowIndex = 0; rowIndex < ROWS; rowIndex++) {
                maxMinWidth = Math.max(
                        maxMinWidth,
                        minWidthOfCell(cells[rowIndex][column.ordinal()], respectContributionFlag));
            }
            
        } else {
            // Only concentrate on the specified row
            maxMinWidth = minWidthOfCell(cells[row.ordinal()][column.ordinal()], respectContributionFlag);
        }
        
        // If this is the center column, we might have an explicit minimal width for that
        if (column == ContainerArea.CENTER && centerCellMinimumSize != null) {
            maxMinWidth = Math.max(maxMinWidth, centerCellMinimumSize.x);
        }
        
        return maxMinWidth;
    }
    
    /**
     * Returns an array containing the height of each row. Symmetry is applied if activated.
     */
    private double[] minRowHeights(final boolean respectContributionFlag) {
        double[] rowHeights = {
                minHeightOfRow(ContainerArea.BEGIN, respectContributionFlag),
                minHeightOfRow(ContainerArea.CENTER, respectContributionFlag),
                minHeightOfRow(ContainerArea.END, respectContributionFlag)
        };
        
        // If we are to be symmetrical, the outer cells need to be the same size
        if (symmetrical) {
            rowHeights[0] = Math.max(rowHeights[0], rowHeights[2]);
            rowHeights[2] = rowHeights[0];
        }
        
        return rowHeights;
    }
    
    /**
     * Returns the minimum height of the given row.
     */
    private double minHeightOfRow(final ContainerArea row, final boolean respectContributionFlag) {
        double maxMinHeight = 0;
        for (int column = 0; column < COLUMNS; column++) {
            maxMinHeight = Math.max(
                    maxMinHeight,
                    minHeightOfCell(cells[row.ordinal()][column], respectContributionFlag));
        }
        
        // If this is the center row, we might have an explicit minimal height for that
        if (row == ContainerArea.CENTER && centerCellMinimumSize != null) {
            maxMinHeight = Math.max(maxMinHeight, centerCellMinimumSize.y);
        }
        
        return maxMinHeight;
    }
    
    /**
     * Takes an array of values and sums them up, inserting gaps between each pair of values bigger than zero.
     */
    private double sumWithGaps(final double[] values) {
        double sum = 0;
        
        // Keep track of how many components > 0 we have
        int activeComponents = 0;
        for (double val : values) {
            if (val > 0) {
                sum += val;
                activeComponents++;
            }
        }
        
        // If there is more than a single value, add necessary gaps
        if (activeComponents > 1) {
            sum += gap * (activeComponents - 1);
        }
        
        return sum;
    }
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Layout Application
    
    /**
     * Computes horizontal coordinates for all cells in the given row, based on the given array of column widths. The
     * array can either be one computed for all rows (tabular mode) or for each row individually. If the row is the
     * center row, this method also updates the X coordinate and the width of {@link #centerCellRect}.
     */
    private void applyWidthsToRow(final ContainerArea row, final double[] colWidths) {
        ElkRectangle cellRectangle = getCellRectangle();
        ElkPadding cellPadding = getPadding();

        // Left column is left-aligned with our content area, right column is right-aligned
        applyWidthToColumn(ContainerArea.BEGIN,
                cellRectangle.x + cellPadding.left,
                colWidths);
        applyWidthToColumn(ContainerArea.END,
                cellRectangle.x + cellRectangle.width - cellPadding.right - colWidths[2],
                colWidths);
        
        // Size of the content area and size of the available space in the content area
        double freeContentAreaWidth = cellRectangle.width - cellPadding.left - cellPadding.right;
        
        if (colWidths[0] > 0) {
            colWidths[0] += gap;
            freeContentAreaWidth -= colWidths[0];
        }
        
        if (colWidths[2] > 0) {
            colWidths[2] += gap;
            freeContentAreaWidth -= colWidths[2];
        }
        
        // Compute how wide the center cell can be
        double centerWidth = Math.max(0, freeContentAreaWidth);
        
        // If the available space is larger than the current size of the center cell, enlarge that thing
        colWidths[1] = Math.max(colWidths[1], freeContentAreaWidth);

        // Place the center cell, possibly enlarging it in the process
        applyWidthToColumn(ContainerArea.CENTER,
                cellRectangle.x + cellPadding.left + colWidths[0] - (colWidths[1] - freeContentAreaWidth) / 2,
                colWidths);
        
        // If this is the center row, remember the center cell's data for the center cell rectangle
        if (row == ContainerArea.CENTER) {
            centerCellRect.width = centerWidth;
            centerCellRect.x = cellRectangle.x + cellPadding.left + (centerWidth - freeContentAreaWidth) / 2;
        }
    }
    
    /**
     * Calls {@link ContainerCell#applyHorizontalLayout(Cell, double, double)} with the given information on all cells
     * in the given column.
     */
    private void applyWidthToColumn(final ContainerArea column, final double x, final double[] colWidths) {
        for (int row = 0; row < ROWS; row++) {
            applyHorizontalLayout(cells[row][column.ordinal()], x, colWidths[column.ordinal()]);
        }
    }
    
    /**
     * Calls {@link ContainerCell#applyVerticalLayout(Cell, double, double)} with the given information on all cells
     * in the given row.
     */
    private void applyHeightToRow(final ContainerArea row, final double y, final double[] rowHeights) {
        for (int column = 0; column < COLUMNS; column++) {
            applyVerticalLayout(cells[row.ordinal()][column], y, rowHeights[row.ordinal()]);
        }
    }
    
}
