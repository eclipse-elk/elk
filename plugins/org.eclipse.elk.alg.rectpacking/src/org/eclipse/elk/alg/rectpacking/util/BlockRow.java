/*******************************************************************************
 * Copyright (c) 2020 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.rectpacking.util;

import java.util.ArrayList;
import java.util.List;

import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.util.ElkUtil;
import org.eclipse.elk.graph.ElkNode;

/**
 * This class corresponds to a row in a block of rectangles.
 * It introduces one abstraction layer for better handling of rows in a block.
 *
 */
public class BlockRow {
    //////////////////////////////////////////////////////////////////
    // Fields.
    
    /**
     * y-coordinate of the row.
     */
    private double x;
    /**
     * y-coordinate of the row.
     */
    private double y;
    /**
     * Width of the row + nodeNodeSpacing
     */
    private double width = 0;
    /**
     * Height of the row.
     * Corresponds to the height of the biggest node + nodeNodeSpacing
     */
    private double height = 0;
    /**
     * Spacing between two nodes.
     */
    private final double nodeNodeSpacing;
    /**
     * Rectangles assigned to the row.
     */
    private List<ElkNode> rects = new ArrayList<>();

    //////////////////////////////////////////////////////////////////
    // Constructors.
    
    /**
     * Create a new BlockRow.
     * @param y The y-coordinate of the row.
     * @param x The x-coordinate of the row.
     * @param nodeNodeSpacing Spacing between two rectangles.
     */
    public BlockRow(final double x, final double y, final double nodeNodeSpacing) {
        this.x = x;
        this.y = y;
        this.nodeNodeSpacing = nodeNodeSpacing;
    }

    //////////////////////////////////////////////////////////////////
    // Public methods.

    /**
     * Adds a rectangle to the row.
     * @param rect The rectangle
     * @return true if the rectangle fits in this row.
     */
    public boolean addRectangle(final ElkNode rect) {
        rect.setX(x + width + (this.rects.isEmpty() ? 0 : nodeNodeSpacing));
        rect.setY(y);
        this.height = Math.max(this.height, rect.getHeight());
        this.width += rect.getWidth() + (this.rects.isEmpty() ? 0 : nodeNodeSpacing);
        this.rects.add(rect);
        return true;
    }
    
    /**
     * Removes a rectangle from the row.
     * Updates the height and width of the row.
     * @param rect The rectangle.
     * @param update If enabled the x and y coordinate of the containing rectangles
     *  and the width and height of the row are updated.
     */
    public void removeRectangle(final ElkNode rect, final boolean update) {
        rects.remove(rect);
        if (update) {
            updateRow();
        }
    }
    
    /**
     * Updates the coordinates of all rectangles and the height and width of the row.
     */
    public void updateRow() {
        double width = 0;
        double height = 0;
        for (ElkNode rect : rects) {
            rect.setX(x + width);
            rect.setY(this.y);
            width += rect.getWidth() + nodeNodeSpacing;
            height = Math.max(height, rect.getHeight() + nodeNodeSpacing);
        }
        this.width = width - nodeNodeSpacing;
        this.height = height - nodeNodeSpacing;
    }

    /**
     * Expands the row to the given width and by the given additional height.
     * @param widthForRow The desired width for this row.
     * @param additionalHeightForRow The additional height for this row.
     * @param index Index of the row.
     */
    public void expand(final double widthForRow, final double additionalHeightForRow, final int index) {
        double additionalWidthForRect = (widthForRow - this.width) / this.rects.size();
        int i = 0;
        this.height += additionalHeightForRow;
        this.width = widthForRow;
        for (ElkNode rect : rects) {
            double oldWidth = rect.getWidth();
            double oldHeight = rect.getHeight();
            rect.setX(rect.getX() + i * additionalWidthForRect);
            rect.setY(rect.getY() + index * additionalHeightForRow);
            rect.setWidth(rect.getWidth() + additionalWidthForRect);
            rect.setHeight(this.height);
            i++;
            double newWidth = rect.getWidth();
            double newHeight = rect.getHeight();
            ElkUtil.translate(rect, new KVector(newWidth, newHeight), new KVector(oldWidth, oldHeight));
        }
    }
    
    //////////////////////////////////////////////////////////////////
    // Getters and setters.

    /**
     * @return the y
     */
    public double getY() {
        return y;
    }

    /**
     * @param y the y to set
     */
    public void setY(double y) {
        this.y = y;
    }

    /**
     * @return the width
     */
    public double getWidth() {
        return width;
    }

    /**
     * @return the height
     */
    public double getHeight() {
        return height;
    }

    /**
     * @return the nodes
     */
    public List<ElkNode> getNodes() {
        return rects;
    }

    /**
     * @return the x
     */
    public double getX() {
        return x;
    }

    /**
     * @param x the x to set
     */
    public void setX(double x) {
        this.x = x;
    }
}
