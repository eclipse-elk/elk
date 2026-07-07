/*******************************************************************************
 * Copyright (c) 2018, 2020 Kiel University and others.
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
import java.util.List;

/**
 * A row of rectangles inside a given bounding box. Tracks the blocks inside the row and the row's width and height.
 * A might have stacks of blocks for ordering them after the compaction step.
 */
public class RectRow {
    //////////////////////////////////////////////////////////////////
    // Fields.
    /** Height of row, given by the highest stack. */
    private double height = 0;
    /** Sum of this row's stack's widths. */
    private double width = 0;
    /** Y-coordinate of this row. */
    private double y;
    /** Spacing between two nodes (used for block and stack spacing) */
    private double nodeNodeSpacing;
    /** This row's blocks. */
    private final List<Block> children = new ArrayList<Block>();
    /**
     * This row's stacks of blocks. 
     * Used during the compaction step for better handling and during the stuffing (expansion) process.
     */
    private List<BlockStack> stacks = new ArrayList<BlockStack>();

    /**
     * The width that was missing to add the first block of the next row into this one.
     */
    private double potentialAdditionalWidthToGetLastBlock = 0;

    //////////////////////////////////////////////////////////////////
    // Constructors.

    /**
     * Creates a row.
     * 
     * @param y The y-coordinate of the row.
     */
    public RectRow(final double y, final double nodeNodeSpacing) {
        this.y = y;
        this.nodeNodeSpacing = nodeNodeSpacing;
    }

    //////////////////////////////////////////////////////////////////
    // Public methods.

    /**
     * Called by one of its assigned blocks when a change was made to respective block like removing or adding a
     * rectangle. By removing or adding a rectangle, the blocks dimensions change which might affect this rows
     * dimensions. This method recalculates its dimensions.
     */
    protected void notifyAboutNodeChange() {
        double totalStackWidth = 0;
        double newMaxHeight = Double.NEGATIVE_INFINITY;
        int index = 0;
        for (Block child : this.children) {
            totalStackWidth += child.getWidth() + (index > 0 ? nodeNodeSpacing : 0);
            newMaxHeight = Math.max(newMaxHeight, child.getHeight());
            index++;
        }

        this.width = totalStackWidth;
        this.height = newMaxHeight;
    }
    
    public void expand(double width, double additionalHeight) {
        double additionalWidth = width - this.width;
        double additionalWidthPerStack = additionalWidth / stacks.size();
        int index = 0;
        for (BlockStack stack : stacks) {
            double additionalHeightForStack = this.getHeight() - stack.getHeight() + additionalHeight;
            stack.setLocation(stack.getX() + index * additionalWidthPerStack, stack.getY());
            stack.expand(additionalWidthPerStack, additionalHeightForStack);
            index++;
        }
    }

    /**
     * Puts blocks with the same x coordinate in stacks for better expansion calculation.
     */
    public void calculateBlockStacks() {
        double currentX = -1;
        for (Block block : children) {
            if (block.getX() != currentX) {
                stacks.add(new BlockStack(block.getX(), this.y, nodeNodeSpacing));
                stacks.get(stacks.size() - 1).addBlock(block);
                currentX = block.getX();
            } else {
                stacks.get(stacks.size() - 1).addBlock(block);
            }
        }
        
    }

    /**
     * Gets the first block of this row.
     */
    public Block getFirstBlock() {
        return this.children.get(0);
    }

    /**
     * Gets the last block of this row.
     */
    public Block getLastBlock() {
        return this.children.get(this.children.size() - 1);
    }

    /**
     * Amount of stacks assigned to this row.
     */
    public int getNumberOfAssignedBlocks() {
        return this.children.size();
    }

    //////////////////////////////////////////////////////////////////
    // Add and remove methods.
    /**
     * Assigns stack to the row. Adjusts height and width of the row. Does not adjust the stack's coordinates.
     * 
     * @param stack The stack to add.
     */
    public void addBlock(final Block block) {

        this.height = Math.max(this.height, block.getHeight());
        this.width += block.getWidth() + (this.children.isEmpty() ? 0 : nodeNodeSpacing);
        this.children.add(block);
    }

    /**
     * Removes specified stack from this row. Adjusts width and height.
     * @param stack The stack to remove.
     */
    public void removeBlock(final Block block) {
        this.children.remove(block);
        this.width -= block.getWidth() + (this.children.isEmpty() ? 0 : nodeNodeSpacing);

        double newMaxHeight = Double.MIN_VALUE;
        for (Block child : this.children) {
            newMaxHeight = Math.max(newMaxHeight, child.getHeight());
        }
        this.height = newMaxHeight;
    }

    //////////////////////////////////////////////////////////////////
    // Getters and setters.
    /**
     * Gets height of this row.
     */
    public double getHeight() {
        return height;
    }

    /**
     * Sets the height of this row.
     * 
     * @param height
     *            new height of this row.
     */
    public void setHeight(final double height) {
        this.height = height;
    }

    /**
     * Gets width of this row.
     */
    public double getWidth() {
        return width;
    }

    /**
     * Sets width of this row.
     * 
     * @param width
     *            width of this row.
     */
    public void setWidth(final double width) {
        this.width = width;
    }

    /**
     * Gets y-coordinate of this row.
     */
    public double getY() {
        return y;
    }

    /**
     * Sets y-coordinate of this row.
     * 
     * @param y The new y-coordinate of this row.
     */
    public void setY(final double y) {
        double yChange = y - this.y;
        for (BlockStack stack : stacks) {
            stack.setLocation(stack.getX(), stack.getY() + yChange);
        }
        this.y = y;
    }

    /**
     * Gets stacks assigned to this row.
     */
    public List<Block> getChildren() {
        return children;
    }

    /**
     * @return the stacks
     */
    public List<BlockStack> getStacks() {
        return stacks;
    }

    /**
     * @return the potentialAdditionalWidthToGetLastBlock
     */
    public double getPotentialAdditionalWidthToGetLastBlock() {
        return potentialAdditionalWidthToGetLastBlock;
    }

    /**
     * @param potentialAdditionalWidthToGetLastBlock the potentialAdditionalWidthToGetLastBlock to set
     */
    public void setPotentialAdditionalWidthToGetLastBlock(double potentialAdditionalWidthToGetLastBlock) {
        this.potentialAdditionalWidthToGetLastBlock = potentialAdditionalWidthToGetLastBlock;
    }

    /**
     * Removes all created stacks. Should only be used when reevaluating the packing in a row.
     */
    public void resetStacks() {
        this.stacks = new ArrayList<BlockStack>();

    }
}
