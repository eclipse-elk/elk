/*******************************************************************************
 * Copyright (c) 2020 Kiel University and others.
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
package org.eclipse.elk.alg.rectpacking.util;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

import org.eclipse.elk.core.math.ElkRectangle;
import org.eclipse.elk.graph.ElkNode;

/**
 * A block is a ordered set of rectangles that shall be layouted left to right. The height of these rectangles is "not that different".
 * The minimum/maximum width and height are defined by their rectangles and are kept up to date on add and remove of rectangles.
 * 
 * Depending on the width and height these rectangles can be layouted.
 * For a given height/width one can calculate the corresponding minimum height/width (if possible).
 * 
 */
public class Block {
    //////////////////////////////////////////////////////////////////
    // Fields.
    private double smallestRectWidth = Double.POSITIVE_INFINITY;
    /** Minimal width + spacing. All rectangles are in one column. */
    private double minWidth;
    /** Current width + spacing. This is corresponds to the maxWidth if no other value is set by the algorithm. */
    private double width;
    /** Minimal height + spacing. All rectangles are in one row. */
    private double minHeight;
    /** Smallest rect height + spacing */
    private double smallestRectHeight = Double.POSITIVE_INFINITY;
    /** Average block height */
    private double averageHeight;
    /** Maximum height + spacing. All rectangles are in one column. */
    private double maxHeight;
    /** Current height + spacing. This is corresponds to the minHeight if no other value is set by the algorithm. */
    private double height;
    /** Rectangles contained in this block. */
    private final List<ElkNode> children = new ArrayList<ElkNode>();
    /** The rectangles in children assigned to rows. */
    private List<BlockRow> rows = new ArrayList<BlockRow>();
    /** X coordinate of this block. */
    private double x;
    /** Y coordinate of this block. */
    private double y;
    /** The row this block is assigned to. */
    private RectRow parentRow;
    /** The stack the block belongs to. */
    private BlockStack stack;
    /** Spacing between two nodes */
    private double nodeNodeSpacing;
    /** Whether the current block has a fixed position and size. */
    private boolean fixed;
    /** Whether the current block has a fixed position. */
    private boolean positionFixed;

    //////////////////////////////////////////////////////////////////
    // Constructor
    
    /**
     * Creates a new block.
     * 
     * @param first
     *            first rectangle of the stack.
     * @param xCoord
     *            x-coordinate of the stack.
     * @param yCoord
     *            y-coordinate of the stack.
     * @param parentRow
     *            row this stack is assigned to.
     * @param nodeNodeSpacing
     *            The spacing between two nodes.
     */
    public Block(final double xCoord, final double yCoord, final RectRow parentRow,
            final double nodeNodeSpacing) {
        this.nodeNodeSpacing = nodeNodeSpacing;
        this.parentRow = parentRow;
        this.x = xCoord;
        this.y = yCoord;

        this.width = 0;
        this.height = 0;
    }

    //////////////////////////////////////////////////////////////////
    // Special methods.
    /**
     * Adds given rectangle to this stack. Adjusts height and width of the stack accordingly. Does not set the
     * coordinates of the ElkNode. Informs parent about size change.
     * 
     * @param rect
     *            rectangle that is added to this stack.
     */
    public void addChild(final ElkNode rect) {
        if (rows.isEmpty()) {
            rows.add(new BlockRow(this.x, this.y, nodeNodeSpacing));
        }
        this.children.add(rect);
        this.rows.get(rows.size() - 1).addRectangle(rect);
        adjustSizeAdd(rect);
    }
    
    public void addChildInNewRow(final ElkNode rect) {
        this.children.add(rect);
        BlockRow lastRow = getLastRow();
        rows.add(new BlockRow(this.x, lastRow.getY() + lastRow.getHeight() + nodeNodeSpacing, nodeNodeSpacing));
        getLastRow().addRectangle(rect);
        adjustSizeAdd(rect);
        
    }

    /**
     * Removes element from bottom of the abstract stack. Sets height and width accordingly. Informs parent about size
     * change.
     * 
     * @param rect
     *            rectangle that is removed from this stack.
     */
    public void removeChild(final ElkNode rect) {
        this.children.remove(rect);
        for (BlockRow row : rows) {
            if (row.getNodes().contains(rect)) {
                row.removeRectangle(rect, true);
                if (row.getNodes().isEmpty()) {
                    rows.remove(row);
                }
                break;
            }
        }
        adjustSizeAfterRemove();
    }

    /**
     * Relocates this block by adjusting its x- and y-coordinate. Passes this change to the children.
     * 
     * @param xCoord
     *            new x-coordinate of this block.
     * @param yCoord
     *            new y-coordinate of this block.
     */
    public void setLocation(final double xCoord, final double yCoord) {
        adjustChildrensXandY(xCoord - this.x, yCoord - this.y);
        for (BlockRow row : rows) {
            row.setX(row.getX() + xCoord - this.x);
            row.setY(row.getY() + yCoord - this.y);
        }
        this.x = xCoord;
        this.y = yCoord;
    }

    /**
     * Shifts rectangles by the given vector.
     * @param xChange Change in x-coordinate.
     * @param yChange Change in y-coordinate.
     */
    private void adjustChildrensXandY(final double xChange, final double yChange) {
        for (ElkNode rect : this.children) {
            rect.setLocation(rect.getX() + xChange, rect.getY() + yChange);
        }
    }

    /**
     * Adjusts size of block after the addition of a rectangle. Notifies parent row.
     * 
     * @param rect
     *            the rectangle that was assigned to this stack.
     */
    private void adjustSizeAdd(final ElkNode rect) {
        double widthOflastRow = getLastRow().getWidth();
        this.smallestRectWidth = Math.min(smallestRectWidth, rect.getWidth());
        this.width = Math.max(width, widthOflastRow);
        this.minWidth = Math.max(minWidth, rect.getWidth() + (this.children.size() == 1 ? 0 : nodeNodeSpacing));
        
        this.smallestRectHeight = Math.min(smallestRectHeight, rect.getHeight());
        this.maxHeight += rect.getHeight() + (this.children.size() == 1 ? 0 : nodeNodeSpacing);
        this.minHeight = Math.max(minHeight, rect.getHeight());
        double totalHeight = rows.size() > 0 ? (rows.size() - 1) * nodeNodeSpacing : 0;
        for (BlockRow row : rows) {
            totalHeight += row.getHeight();
        }
        this.height = totalHeight;
        this.averageHeight = maxHeight / this.children.size()
                - nodeNodeSpacing * ((double) (this.children.size() - 1) / ((double) this.children.size()));
        this.parentRow.notifyAboutNodeChange();
    }
    

    /**
     * Searches and returns the smallest possible width in which the rectangles of this block can be placed without exceeding the given height.
     * @param height The height.
     * @return The width of the rectangle block with specified height.
     */
    public double getWidthForTargetHeight(final double height) {
        // Check whether the block would just fit if all rectangles are drawn below each other.
        if (this.maxHeight <= height) {
            return this.minWidth;
        }
        
        // Check if the minimal width of the block is enough. This can be the case if some rows have more than one rectangle.
        if (placeRectsIn(this.minWidth, height, false)) {
            return this.minWidth;
        }
        
        // Binary search between minWidth and width.
        double upperBound = this.width;
        double lowerBound = this.minWidth;
        double viableWidth = this.width;
        double newWidth = (upperBound - lowerBound) / 2 + lowerBound;
        
        while (lowerBound + 1 < upperBound) {
            if (placeRectsIn(newWidth, height, false)) {
                viableWidth = newWidth;
                upperBound = newWidth;
            } else {
                lowerBound = newWidth;
            }
            newWidth = (upperBound - lowerBound) / 2 + lowerBound;
        }
        return viableWidth;
    }
    
    /**
     * Returns the height of this block if you would place it in the given width.
     * @param width The width
     * @return The height of this block defined by the width.
     */
    public double getHeightForTargetWidth(final double width) {
        ElkRectangle bounds = placeRectsIn(width, false);
        return bounds.height;
    }
    
    /**
     * Checks whether rectangles can fit in a given width.
     * Optionally the rectangles can be directly placed.
     * @param width The width
     * @param placeRects Whether the rectangles should be directly placed.
     * @return The bounds of this block if you would place all rectangles in the given width.
     */
    private ElkRectangle placeRectsIn(final double width, final boolean placeRects) {
        // Next x coordinate at which a rect shall be placed in a blockrow
        double currentX = 0;
        // Next y coordinate at which a rect shall be placed in a blockrow
        double currentY = this.y;
        double currentWidth = 0;
        double currentHeight = 0;
        double maxHeightInRow = 0;
        double widthInRow = 0;
        int row = 0;
        if (placeRects) {
            rows.clear();
            rows.add(new BlockRow(this.x, this.y, nodeNodeSpacing));
        }
        // Current index in row.
        int index = 0;
        for (ElkNode rect : children) {
            if (currentX + rect.getWidth() + (index > 0 ? nodeNodeSpacing : 0) > width
                    && maxHeightInRow > 0) {
                // Case new row
                currentX = 0;
                currentY += maxHeightInRow + nodeNodeSpacing;
                currentWidth = Math.max(currentWidth, widthInRow);
                currentHeight += maxHeightInRow + nodeNodeSpacing;
                maxHeightInRow = 0;
                widthInRow = 0;
                if (placeRects) {
                    row++;
                    rows.add(new BlockRow(this.x, currentY, nodeNodeSpacing));
                }
                // Reset current index in row
                index = 0;
            }
            widthInRow += rect.getWidth() + (index > 0 ? nodeNodeSpacing : 0);
            maxHeightInRow = Math.max(maxHeightInRow, rect.getHeight());
            if (placeRects) {
                rows.get(row).addRectangle(rect);
            }
            currentX += rect.getWidth() + (index > 0 ? nodeNodeSpacing : 0);
            index++;
        }
        currentWidth = Math.max(currentWidth, widthInRow);
        currentHeight += maxHeightInRow;
        if (placeRects) {
            this.width = currentWidth;
            this.height = currentHeight;
            this.parentRow.notifyAboutNodeChange();
        }
        return new ElkRectangle(x, y, currentWidth, currentHeight);
    }
    
    /**
     * Checks whether the rectangles of this block can be placed width a given width and height.
     * Optionally places the rectangles.
     * @param width The width.
     * @param height The height.
     * @param placeRects Whether the rectangles should be directly placed.
     * @return true if the rectangles could be placed within the given width and height.
     */
    private boolean placeRectsIn(final double width, final double height, final boolean placeRects) {
        ElkRectangle bounds = placeRectsIn(width, placeRects);
        return bounds.width <= width && bounds.height <= height;
    }
    
    /**
     * Places the rectangles of this block in a given width and height.
     * @param width The width
     * @param height The height
     * @return true if placing the rectangles exceeded the given width or height.
     */
    public boolean placeRectsIn(final double width, final double height) {
        return placeRectsIn(width, height, true);
    }
    

    /**
     * Places the rectangles of this block in a given width.
     * @param width The width
     */
    public boolean placeRectsIn(final double width) {
        double oldWidth = this.width;
        double oldHeight = this.height;
        ElkRectangle bounds =  placeRectsIn(width, true);
        return bounds.width != oldWidth || bounds.height != oldHeight;
    }

    /**
     * Adjusts size of stack after the removal of a rectangle. Notifies parent.
     */
    private void adjustSizeAfterRemove() {
        // Calculate new width and height
        double newWidth = 0;
        double newHeight = 0;
        // Delete empty rows and calculate new width and height.
        List<BlockRow> rowsToDelete = new LinkedList<>();
        int index = 0;
        for (BlockRow row : this.rows) {
            if (row.getNodes().isEmpty()) {
                rowsToDelete.add(row);
            } else {
                newWidth = Math.max(newWidth, row.getWidth());
                newHeight += row.getHeight() + (index > 0 ? nodeNodeSpacing : 0);
            }
            index++;
        }
        rows.removeAll(rowsToDelete);
        this.height = newHeight;
        this.width = newWidth;
        
        // Calculate new smallest/biggest rectangle, minimum/maximum width/height.
        this.minWidth = 0;
        this.minHeight = 0;
        this.maxHeight = 0;
        this.smallestRectHeight = Double.POSITIVE_INFINITY;
        this.smallestRectWidth = Double.POSITIVE_INFINITY;
        for (ElkNode rect : children) {
            this.smallestRectWidth = Math.min(smallestRectWidth, rect.getWidth());
            this.minWidth = Math.max(minWidth, rect.getWidth());
            this.minHeight = Math.max(minHeight, rect.getHeight());
            this.smallestRectHeight = Math.min(smallestRectHeight, rect.getHeight());
            this.maxHeight += rect.getHeight() + nodeNodeSpacing;
        }
        this.averageHeight = maxHeight / this.children.size()
                - nodeNodeSpacing * ((double) (this.children.size() - 1) / ((double) this.children.size()));
        
        this.parentRow.notifyAboutNodeChange();
    }

    /**
     * Stuff the empty space by distributing it to all rows if resizes the block by the given additional width and height.
     * @param additionalWidthPerBlock The additional width.
     * @param additionalHeightForBlock The additional height.
     */
    public void expand(final double additionalWidthPerBlock, final double additionalHeightForBlock) {
        double widthForRow = this.width + additionalWidthPerBlock;
        this.width += additionalWidthPerBlock;
        this.height += additionalHeightForBlock;
        double additionalHeightForRow = additionalHeightForBlock / this.rows.size();
        int index = 0;
        for (BlockRow row : rows) {
            row.expand(widthForRow, additionalHeightForRow, index);
            index++;
        }
    }

    //////////////////////////////////////////////////////////////////
    // Getters and setters.
    /**
     * Gets width of the stack.
     */
    public double getWidth() {
        return width;
    }

    /**
     * Sets width of the stack.
     * 
     * @param width
     *            new width of this stack.
     */
    public void setWidth(final double width) {
        this.width = width;
    }

    /**
     * Gets height of the stack.
     */
    public double getHeight() {
        return height;
    }

    /**
     * Sets height of the stack.
     * 
     * @param height
     *            new height of this stack.
     */
    public void setHeight(final double height) {
        this.height = height;
    }

    /**
     * Gets this stack's children.
     */
    public List<ElkNode> getChildren() {
        return children;
    }

    /**
     * Gets this stack's x-coordinate.
     */
    public double getX() {
        return x;
    }

    /**
     * Gets this stack's y/coordinate.
     */
    public double getY() {
        return y;
    }

    /**
     * Gets the row this stack is assigned to.
     */
    public RectRow getParentRow() {
        return parentRow;
    }

    /**
     * Sets the row this stack belongs to.
     * 
     * @param parentRow
     *            row this stack is assigned to.
     */
    public void setParentRow(final RectRow parentRow) {
        this.parentRow = parentRow;
    }

    /**
     * @return the minWidth
     */
    public double getMinWidth() {
        return minWidth;
    }

    /**
     * @return the minHeight
     */
    public double getMinHeight() {
        return minHeight;
    }

    /**
     * @return the maxHeight
     */
    public double getMaxHeight() {
        return maxHeight;
    }

    /**
     * @return the smallestRectHeight
     */
    public double getSmallestRectHeight() {
        return smallestRectHeight;
    }

    /**
     * @return the averageHeight
     */
    public double getAverageHeight() {
        return averageHeight;
    }

    /**
     * @return the fixed
     */
    public boolean isFixed() {
        return fixed;
    }

    /**
     * @param fixed the fixed to set
     */
    public void setFixed(boolean fixed) {
        this.fixed = fixed;
    }

    /**
     * @return the positionFixed
     */
    public boolean isPositionFixed() {
        return positionFixed;
    }

    /**
     * @param positionFixed the positionFixed to set
     */
    public void setPositionFixed(boolean positionFixed) {
        this.positionFixed = positionFixed;
    }

    /**
     * @return the rows
     */
    public List<BlockRow> getRows() {
        return rows;
    }

    /**
     * @return The x position at which the next rect would be added to the last row.
     */
    public double getLastRowNewX() {
        BlockRow lastRow = getLastRow();
        return lastRow.getX() + lastRow.getWidth();
    }

    /**
     * @return The y position at which the rectangles of the last row are placed.
     */
    public double getLastRowY() {
        BlockRow lastRow = getLastRow();
        return lastRow.getY();
    }
    
    public BlockRow getLastRow() {
        return this.rows.get(rows.size() - 1);
    }

    /**
     * @return the stack
     */
    public BlockStack getStack() {
        return stack;
    }

    /**
     * @param stack the stack to set
     */
    public void setStack(BlockStack stack) {
        this.stack = stack;
    }

    /**
     * @return the smallestRectWidth
     */
    public double getSmallestRectWidth() {
        return smallestRectWidth;
    }

    /**
     * @param smallestRectWidth the smallestRectWidth to set
     */
    public void setSmallestRectWidth(double smallestRectWidth) {
        this.smallestRectWidth = smallestRectWidth;
    }

    /**
     * Recalculates the size of this block.
     */
    public void resetBlock() {
        adjustSizeAfterRemove();
    }
}
