/*******************************************************************************
 * Copyright (c) 2017, 2019 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.common.nodespacing.cellsystem;

import java.util.List;
import java.util.Objects;

import org.eclipse.elk.alg.common.nodespacing.internal.NodeLabelLocation;
import org.eclipse.elk.core.math.ElkPadding;
import org.eclipse.elk.core.math.ElkRectangle;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.util.adapters.GraphAdapters.LabelAdapter;

import com.google.common.collect.Lists;

/**
 * A cell which manages the size and placement of labels. It inserts a customizable gap between each of the labels and
 * calculates its minimum size to be the area required to place all of its labels. The horizontal and vertical alignment
 * control what happens if the cell's actual size is larger than required by its labels. Call
 * {@link #applyLabelLayout()} to cause the cell to go ahead and assign positions to all of its labels.
 * 
 * <p>
 * The cell can operate in horizontal and vertical layout mode (with horizontal being the default). Horizontal layout
 * mode is used for horizontal layout directions and assumes that labels are horizontal lines of text that need to be
 * stacked atop one another. Vertical layout mode assumes that labels are vertical columns of text that need to be
 * stacked next to each other. If in doubt, use horizontal layout mode.
 * </p>
 */
public class LabelCell extends Cell {
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Properties

    /** Whether we operate in horizontal or vertical layout mode. */
    private final boolean horizontalLayoutMode;
    /** Horizontal alignment of labels. */
    private HorizontalLabelAlignment horizontalAlignment = HorizontalLabelAlignment.CENTER;
    /** Vertical alignment of labels. */
    private VerticalLabelAlignment verticalAlignment = VerticalLabelAlignment.CENTER;
    /** The gap inserted between two consecutive labels. */
    private final double gap;
    /** The labels in our cell. */
    private List<LabelAdapter<?>> labels = Lists.newArrayListWithCapacity(2);
    /** Our minimum size. This is basically the space we need at minimum to place our labels. */
    private KVector minimumContentAreaSize = new KVector();
    

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Constructors
    
    /**
     * Constructs a new instance with the given properties.
     * 
     * @param gap
     *            gap between labels.
     */
    public LabelCell(final double gap) {
        this.gap = gap;
        this.horizontalLayoutMode = true;
    }
    
    /**
     * Constructs a new instance with the given properties.
     * 
     * @param gap
     *            gap between labels.
     * @param horizontalLayoutMode
     *            whether the cell should operate in horizontal or vertical layout mode.
     */
    public LabelCell(final double gap, final boolean horizontalLayoutMode) {
        this.gap = gap;
        this.horizontalLayoutMode = horizontalLayoutMode;
    }
    
    
    /**
     * Constructs a new instance with the given properties.
     * 
     * @param gap
     *            gap between labels.
     * @param nodeLabelLocation
     *            the label location represented by this cell. This determines things like alignment.
     */
    public LabelCell(final double gap, final NodeLabelLocation nodeLabelLocation) {
        this(gap, nodeLabelLocation, true);
    }
    
    /**
     * Constructs a new instance with the given properties.
     * 
     * @param gap
     *            gap between labels.
     * @param nodeLabelLocation
     *            the label location represented by this cell. This determines things like alignment.
     * @param horizontalLayoutMode 
     *            Whether we operate in horizontal or vertical layout mode.
     */
    public LabelCell(final double gap, final NodeLabelLocation nodeLabelLocation, final boolean horizontalLayoutMode) {
        this.gap = gap;
        this.horizontalLayoutMode = horizontalLayoutMode;
        this.horizontalAlignment = nodeLabelLocation.getHorizontalAlignment();
        this.verticalAlignment = nodeLabelLocation.getVerticalAlignment();
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Getters / Setters
    
    /**
     * Returns the horizontal alignment of labels.
     */
    public HorizontalLabelAlignment getHorizontalAlignment() {
        return horizontalAlignment;
    }

    /**
     * Sets the horizontal alignment of labels. Can be chained with further method calls.
     */
    public LabelCell setHorizontalAlignment(final HorizontalLabelAlignment newHorizontalAlignment) {
        Objects.requireNonNull(newHorizontalAlignment, "Horizontal alignment cannot be null");
        this.horizontalAlignment = newHorizontalAlignment;
        return this;
    }

    /**
     * Returns the vertical alignment of labels.
     */
    public VerticalLabelAlignment getVerticalAlignment() {
        return verticalAlignment;
    }

    /**
     * Sets the vertical alignment of labels. Can be chained with further method calls.
     */
    public LabelCell setVerticalAlignment(final VerticalLabelAlignment newVerticalAlignment) {
        Objects.requireNonNull(newVerticalAlignment, "Vertical alignment cannot be null");
        this.verticalAlignment = newVerticalAlignment;
        return this;
    }

    /**
     * Returns the list of labels in this cell. The list should not be modified, except to reorder the labels. To
     * add labels to the cell, call {@link #addLabel(LabelAdapter)} instead.
     * 
     * @return unmodifiable list of labels.
     */
    public List<LabelAdapter<?>> getLabels() {
        return labels;
    }
    

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Cell
    
    @Override
    public double getMinimumWidth() {
        ElkPadding padding = getPadding();
        return minimumContentAreaSize.x + padding.left + padding.right;
    }

    @Override
    public double getMinimumHeight() {
        ElkPadding padding = getPadding();
        return minimumContentAreaSize.y + padding.top + padding.bottom;
    }
    

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Adding Labels
    
    /**
     * Adds a label to this cell.
     * 
     * @param label
     *            the label.
     */
    public void addLabel(final LabelAdapter<?> label) {
        labels.add(label);
        
        // Update our minimum size (how this works depends on the layout mode)
        KVector labelSize = label.getSize();
        
        if (horizontalLayoutMode) {
            minimumContentAreaSize.x = Math.max(minimumContentAreaSize.x, labelSize.x);
            minimumContentAreaSize.y += labelSize.y;
            
            // If this is not our first label, insert a gap
            if (labels.size() > 1) {
                minimumContentAreaSize.y += gap;
            }
        } else {
            minimumContentAreaSize.x += labelSize.x;
            minimumContentAreaSize.y = Math.max(minimumContentAreaSize.y, labelSize.y);
            
            // If this is not our first label, insert a gap
            if (labels.size() > 1) {
                minimumContentAreaSize.x += gap;
            }
        }
    }
    
    /**
     * Returns whether there are any labels in this cell.
     * 
     * @return {@code true} if this cell contains labels, {@code false} otherwise.
     */
    public boolean hasLabels() {
        return labels.size() > 0;
    }
    

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Label Layout
    
    /**
     * Actually applies coordinates to the labels based on this cell's content area rectangle.
     */
    public void applyLabelLayout() {
        // Delegate to the methods specialized for our two layout modes
        if (horizontalLayoutMode) {
            applyHorizontalModeLabelLayout();
        } else {
            applyVerticalModeLabelLayout();
        }
    }

    /**
     * Implementation of {@link #applyLabelLayout()} for horizontal layout mode.
     */
    private void applyHorizontalModeLabelLayout() {
        ElkRectangle cellRect = getCellRectangle();
        ElkPadding cellPadding = getPadding();
        
        // Calculate our starting y coordinate
        double yPos = cellRect.y;
        
        if (verticalAlignment == VerticalLabelAlignment.CENTER) {
            yPos += (cellRect.height - minimumContentAreaSize.y) / 2;
        } else if (verticalAlignment == VerticalLabelAlignment.BOTTOM) {
            yPos += cellRect.height - minimumContentAreaSize.y;
        }
        
        // Place them labels, I say!
        for (LabelAdapter<?> label : labels) {
            KVector labelSize = label.getSize();
            KVector labelPos = new KVector();
            
            // Y coordinate
            labelPos.y = yPos;
            yPos += labelSize.y + gap;
            
            // X coordinate
            switch (horizontalAlignment) {
            case LEFT:
                labelPos.x = cellRect.x + cellPadding.left;
                break;
                
            case CENTER:
                labelPos.x = cellRect.x + cellPadding.left + (cellRect.width - labelSize.x) / 2;
                break;
                
            case RIGHT:
                labelPos.x = cellRect.x + cellRect.width - cellPadding.right - labelSize.x;
                break;
            }
            
            // Apply position
            label.setPosition(labelPos);
        }
    }

    /**
     * Implementation of {@link #applyLabelLayout()} for vertical layout mode.
     */
    private void applyVerticalModeLabelLayout() {
        ElkRectangle cellRect = getCellRectangle();
        ElkPadding cellPadding = getPadding();
        
        // Calculate our starting y coordinate
        double xPos = cellRect.x;
        
        if (horizontalAlignment == HorizontalLabelAlignment.CENTER) {
            xPos += (cellRect.width - minimumContentAreaSize.x) / 2;
        } else if (horizontalAlignment == HorizontalLabelAlignment.RIGHT) {
            xPos += cellRect.width - minimumContentAreaSize.x;
        }
        
        // Place them labels, I say!
        for (LabelAdapter<?> label : labels) {
            KVector labelSize = label.getSize();
            KVector labelPos = new KVector();
            
            // X coordinate
            labelPos.x = xPos;
            xPos += labelSize.x + gap;
            
            // Y coordinate
            switch (verticalAlignment) {
            case TOP:
                labelPos.y = cellRect.y + cellPadding.top;
                break;
                
            case CENTER:
                labelPos.y = cellRect.y + cellPadding.top + (cellRect.height - labelSize.y) / 2;
                break;
                
            case BOTTOM:
                labelPos.y = cellRect.y + cellRect.height - cellPadding.bottom - labelSize.y;
                break;
            }
            
            // Apply position
            label.setPosition(labelPos);
        }
    }
    
}
