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

import org.eclipse.elk.core.math.ElkRectangle;

/**
 * A container cell contains other cells. How it contains them depends on the actual container cell. This abstract class
 * simply adds methods that allow clients to tell the container to lay out its children.
 * 
 * @see StripContainerCell
 * @see GridContainerCell
 */
public abstract class ContainerCell extends Cell {

    /**
     * Compute x coordinates and widths of children.
     */
    public abstract void layoutChildrenHorizontally();

    /**
     * Compute y coordinates and heights of children.
     */
    public abstract void layoutChildrenVertically();

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Utility Methods

    /**
     * Returns the minimum width of the given cell.
     * 
     * @param cell
     *            the cell whose minimum width to return.
     * @param respectContributionFlag
     *            if {@code true}, zero is returned for cells that don't have the width contribution flag set.
     * @return minimum width.
     */
    protected static double minWidthOfCell(final Cell cell, final boolean respectContributionFlag) {
        // If there's no cell, there's no minimum with
        if (cell == null) {
            return 0;
        }

        // If the cell doesn't have its contribution flag activated, there's no minimum width
        if (respectContributionFlag && !cell.isContributingToMinimumWidth()) {
            return 0;
        }

        // If the cell is an atomic cell with a content area of no width, there's no minimum width
        if (cell instanceof AtomicCell) {
            if (((AtomicCell) cell).getMinimumContentAreaSize().x == 0) {
                return 0;
            }
        }

        return cell.getMinimumWidth();
    }

    /**
     * Returns the minimum height of the given cell.
     * 
     * @param cell
     *            the cell whose minimum height to return.
     * @param respectContributionFlag
     *            if {@code true}, zero is returned for cells that don't have the width contribution flag set.
     * @return minimum height.
     */
    protected static double minHeightOfCell(final Cell cell, final boolean respectContributionFlag) {
        // If there's no cell, there's no minimum height
        if (cell == null) {
            return 0;
        }

        // If the cell doesn't have its contribution flag activated, there's no minimum height
        if (respectContributionFlag && !cell.isContributingToMinimumHeight()) {
            return 0;
        }

        // If the cell is an atomic cell with a content area of no height, there's no minimum height
        if (cell instanceof AtomicCell) {
            if (((AtomicCell) cell).getMinimumContentAreaSize().y == 0) {
                return 0;
            }
        }

        return cell.getMinimumHeight();
    }
    
    /**
     * Applies the given horizontal layout information to the given cell if it's not {@code null}.
     * 
     * @param cell
     *            the cell to apply the layout information to.
     * @param x
     *            the cell's new x coordinate.
     * @param width
     *            the cell's new width.
     */
    protected void applyHorizontalLayout(final Cell cell, final double x, final double width) {
        if (cell != null) {
            ElkRectangle cellRect = cell.getCellRectangle();
            cellRect.x = x;
            cellRect.width = width;
        }
    }
    
    /**
     * Applies the given vertical layout information to the given cell if it's not {@code null}.
     * 
     * @param cell
     *            the cell to apply the layout information to.
     * @param y
     *            the cell's new y coordinate.
     * @param height
     *            the cell's new height.
     */
    protected void applyVerticalLayout(final Cell cell, final double y, final double height) {
        if (cell != null) {
            ElkRectangle cellRect = cell.getCellRectangle();
            cellRect.y = y;
            cellRect.height = height;
        }
    }

}
