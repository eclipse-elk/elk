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

/**
 * Inside of a row blocks are stacked in the compaction step.
 * This class represents a stack of blocks. All blocks in this stack are placed below each other.
 */
public class BlockStack {
    //////////////////////////////////////////////////////////////////
    // Fields.
    
    /**
     * The blocks assigned to this stack.
     */
    private List<Block> blocks;
    /**
     * x-coordinate of this stack.
     */
    private double x;
    /**
     * y-coordinate of this stack.
     */
    private double y;
    /**
     * Width of this stack.
     */
    private double width;
    /**
     * Height of this stack.
     */
    private double height;
    /**
     * Spacing between two nodes.
     */
    private double nodeNodeSpacing;

    //////////////////////////////////////////////////////////////////
    // Constructors.
    
    /**
     * Creates a new stack at a specific location.
     * The y-coordinate should correspond to the y-coordinate of the row.
     * @param x The x-coordinate.
     * @param y The y-coordinate.
     */
    public BlockStack(final double x, final double y, final double nodeNodeSpacing) {
        blocks = new ArrayList<>();
        this.x = x;
        this.y = y;
        this.nodeNodeSpacing = nodeNodeSpacing;
    }

    //////////////////////////////////////////////////////////////////
    // Public methods.
    
    /**
     * Adds a block to this stack.
     * @param block The block to add.
     */
    public void addBlock(final Block block) {
        block.setStack(this);
        this.width = Math.max(this.width, block.getWidth());
        this.height += block.getHeight() + (this.blocks.isEmpty() ? 0 : nodeNodeSpacing);
        blocks.add(block);
    }
    
    /**
     * Recalculates the dimension (width and height) of the stack.
     */
    public void updateDimension() {
        double height = 0;
        double width = 0;
        int index = 0;
        for (Block block : blocks) {
            width = Math.max(width, block.getWidth());
            height += block.getHeight() + (index > 0 ? nodeNodeSpacing : 0);
            index++;
        }
        this.height = height;
        this.width = width;
    }
    
    /**
     * Updates the location of this stack including all of its blocks.
     * @param x The x-coordinate.
     * @param y The y-coordinate.
     */
    public void setLocation(final double x, final double y) {
        double xDiff = x - this.x;
        double yDiff = y - this.y;
        for (Block block : blocks) {
            block.setLocation(block.getX() + xDiff, block.getY() + yDiff);
        }
        this.x = x;
        this.y = y;
    }

    /**
     * Returns the minimal width for this stack that would still fit the given height.
     * @param height The maximum height for this stack.
     * @return The width of the stack such that all blocks fit the height.
     */
    public double getWidthForFixedHeight(final double height) {
        // One element special case.
        if (blocks.size() == 1) {
            return blocks.get(0).getWidthForTargetHeight(height);
        }
        
        // Try some kind of binary search between the widest block length and the minWidth.
        double minWidth = getMinimumWidth();
        double totalHeight = 0;
        double upperBound = this.width;
        double lowerBound = minWidth;
        double viableWidth = this.width;
        double newWidth = (upperBound - lowerBound) / 2 + lowerBound;
        while (lowerBound + 1 < upperBound) {
            totalHeight = 0;
            for (Block block : blocks) {
                totalHeight += block.getHeightForTargetWidth(newWidth);
            }
            if (totalHeight < height) {
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
     * Places all blocks in this stack in the specified target width.
     * Also updates the width and height of the stack.
     * @param targetWidth
     */
    public void placeRectsIn(final double targetWidth) {
        double currentY = this.getY();
        double currentHeight = 0;
        double currentWidth = 0;
        for (Block block : blocks) {
            block.setLocation(this.x, currentY);
            block.placeRectsIn(targetWidth);
            currentWidth = Math.max(currentWidth, block.getWidth());
            currentY += block.getHeight() + nodeNodeSpacing;
            currentHeight = currentY;
        }
        this.width = currentWidth;
        this.height = currentHeight;
        
    }
    
    /**
     * Expands this stack by the given additional width and height.
     * The components of this stack are enlarged rather than recalculating the layout.
     * @param additionalWidth The additional width.
     * @param additionalHeight The additional height.
     */
    public void expand(final double additionalWidth, final double additionalHeight) {
        int index = 0;
        double additionalHeightPerBlock = additionalHeight / blocks.size();
        for (Block block : blocks) {
            block.setLocation(block.getX(),  block.getY() + index * additionalHeightPerBlock);
            block.expand(this.getWidth() - block.getWidth() + additionalWidth, additionalHeightPerBlock);
            index++;
        }
    }

    //////////////////////////////////////////////////////////////////
    // Helper methods.
    
    /**
     * Returns the minimum width of this stack.
     * @return The minimum width.
     */
    private double getMinimumWidth() {
        double minWidth = 0;
        for (Block block : blocks) {
            minWidth = Math.max(minWidth, block.getMinWidth());
        }
        return minWidth;
    }

    //////////////////////////////////////////////////////////////////
    // Getters and setters.

    /**
     * @return the blocks
     */
    public List<Block> getBlocks() {
        return blocks;
    }

    /**
     * @return the x
     */
    public double getX() {
        return x;
    }

    /**
     * @return the y
     */
    public double getY() {
        return y;
    }

    /**
     * @return the width
     */
    public double getWidth() {
        return width;
    }

    /**
     * @param width the width to set
     */
    public void setWidth(double width) {
        this.width = width;
    }

    /**
     * @return the height
     */
    public double getHeight() {
        return height;
    }

    /**
     * @param height the height to set
     */
    public void setHeight(double height) {
        this.height = height;
    }
}
