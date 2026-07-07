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
import org.eclipse.elk.core.math.KVector;

/**
 * An atomic cell is a simple cell that simply holds a size. Its minimum size can be directly manipulated. This is
 * basically a placeholder for other things to be placed in the cell system later.
 */
public class AtomicCell extends Cell {
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Properties

    /** The minimum size of a cell's content area (that is, this excludes the padding). */
    private KVector minimumContentAreaSize = new KVector();
    

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Getters / Setters
    
    /**
     * Returns the minimum content area size, to be modified by clients. Note that this size does not include any
     * padding.
     */
    public KVector getMinimumContentAreaSize() {
        return minimumContentAreaSize;
    }
    
    /**
     * Sets the cell's minimum size.
     * 
     * @param newMinimumContentAreaSize
     *            the new minimum content area size.
     * @param includesPadding
     *            if {@code true}, the new size includes padding that needs to be subtracted.
     */
    public void setMinimumContentAreaSize(final KVector newMinimumContentAreaSize, final boolean includesPadding) {
        if (includesPadding) {
            ElkPadding padding = getPadding();
            this.minimumContentAreaSize.x = newMinimumContentAreaSize.x - padding.left - padding.right;
            this.minimumContentAreaSize.y = newMinimumContentAreaSize.y - padding.top - padding.bottom;
        } else {
            this.minimumContentAreaSize.set(newMinimumContentAreaSize);
        }
    }
    

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Cell
    
    /* (non-Javadoc)
     * @see org.eclipse.elk.core.util.nodespacing.internal.cells.Cell#getMinimumWidth()
     */
    @Override
    public double getMinimumWidth() {
        ElkPadding padding = getPadding();
        return minimumContentAreaSize.x + padding.left + padding.right;
    }

    /* (non-Javadoc)
     * @see org.eclipse.elk.core.util.nodespacing.internal.cells.Cell#getMinimumHeight()
     */
    @Override
    public double getMinimumHeight() {
        ElkPadding padding = getPadding();
        return minimumContentAreaSize.y + padding.top + padding.bottom;
    }
    
}
