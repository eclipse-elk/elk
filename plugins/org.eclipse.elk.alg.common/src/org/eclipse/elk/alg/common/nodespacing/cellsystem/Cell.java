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

/**
 * A cell is the basic component of the cell system. Each cell has a padding, which determines the amount of space
 * between its content area and its border. It also has a minimum width and height, which is the minimum size it
 * would like to be in the final layout. {@link ContainerCell Container cells} will use that information to compute
 * their own minimum size. Whether or not a cell contributes to that is controlled through its flags, for width and
 * height separately. Finally, a cell has a rectangle which desribes its actual position and size. While these
 * information can be set manually, they will often be computed by container cells.
 */
public abstract class Cell {
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Properties
    
    /** A cell has a padding. */
    private ElkPadding padding = new ElkPadding();
    /** The actual size and position of the cell. Includes the padding. */
    private ElkRectangle cellRectangle = new ElkRectangle();
    /** Whether the cell contributes to the minimum width calculation of a container cell or not. */
    private boolean contributesToMinimumWidth = false;
    /** Whether the cell contributes to the minimum height calculation of a container cell or not. */
    private boolean contributesToMinimumHeight = false;
    

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Getters / Setters

    /**
     * Returns this cell's padding, to be modified by the caller.
     */
    public ElkPadding getPadding() {
        return padding;
    }
    
    /**
     * Returns a rectangle that describes the cell's size and position, including padding, to be modified by the
     * caller.
     */
    public ElkRectangle getCellRectangle() {
        return cellRectangle;
    }
    
    /**
     * Checks whether this cell should be included when calculating a container cell's minimum width.
     */
    public boolean isContributingToMinimumWidth() {
        return contributesToMinimumWidth;
    }
    
    /**
     * Sets whether this cell should be included when calculating a container cell's minimum width.
     */
    public void setContributesToMinimumWidth(final boolean contributesToMinimumWidth) {
        this.contributesToMinimumWidth = contributesToMinimumWidth;
    }
    
    /**
     * Checks whether this cell should be included when calculating a container cell's minimum height.
     */
    public boolean isContributingToMinimumHeight() {
        return contributesToMinimumHeight;
    }
    
    /**
     * Sets whether this cell should be included when calculating a container cell's minimum height.
     */
    public void setContributesToMinimumHeight(final boolean contributesToMinimumHeight) {
        this.contributesToMinimumHeight = contributesToMinimumHeight;
    }
    

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Abstract Methods

    /**
     * Returns this cell's minimum width, including padding.
     */
    public abstract double getMinimumWidth();

    /**
     * Returns this cell's minimum height, including padding.
     */
    public abstract double getMinimumHeight();
    
}
