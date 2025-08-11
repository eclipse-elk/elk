/*******************************************************************************
 * Copyright (c) 2018, 2020 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.alg.rectpacking.p2packing;

import java.util.ArrayList;
import java.util.List;

import org.eclipse.elk.alg.rectpacking.options.RectPackingOptions;
import org.eclipse.elk.alg.rectpacking.util.Block;
import org.eclipse.elk.alg.rectpacking.util.RectRow;
import org.eclipse.elk.graph.ElkNode;

/**
 * Class that offers methods that calculate the initial placement of rectangles in {@link RowFillingAndCompaction}.
 * 
 * @see RowFillingAndCompaction
 */
public final class InitialPlacement {

    //////////////////////////////////////////////////////////////////
    // Private Constructor.
    private InitialPlacement() {
    }

    //////////////////////////////////////////////////////////////////
    // Static methods.
    /**
     * Simply places the rectangles as {@link RectRow}s onto the drawing area, bounded by the calculated bounding box
     * width.
     * 
     * @param rectangles The rectangles to be placed.
     * @param boundingWidthThe width of the bounding box.
     * @param nodeNodeSpacing The spacing between two nodes.
     * @return returns the rows in which the rectangles were placed.
     */
    protected static List<RectRow> place(final List<ElkNode> rectangles, final double boundingWidth, final double nodeNodeSpacing) {
        List<RectRow> rows = new ArrayList<RectRow>();
        RectRow row = new RectRow(0, nodeNodeSpacing);
        double drawingHeight = 0;
        row.addBlock(new Block(0, 0, row, nodeNodeSpacing));
        double currentWidth = 0;
        
        for (ElkNode rect : rectangles) {
            // Check whether current rectangle can be added to the last block
            Block block = row.getLastBlock();
            double potentialRowWidth = currentWidth + rect.getWidth() + 
                    (row.getChildren().get(0).getChildren().isEmpty() ? 0 : nodeNodeSpacing);
            if (potentialRowWidth > boundingWidth || rect.getProperty(RectPackingOptions.IN_NEW_ROW)) {
                // Add rect in new block in new row.
                currentWidth = 0;
                drawingHeight += row.getHeight() + nodeNodeSpacing;
                rows.add(row);
                row = new RectRow(drawingHeight, nodeNodeSpacing);
                block = new Block(0, row.getY(), row, nodeNodeSpacing);
                row.addBlock(block);
                currentWidth = 0;
            }
            
            if (block.getChildren().isEmpty()
                    || !rect.getParent().getProperty(RectPackingOptions.PACKING_COMPACTION_ROW_HEIGHT_REEVALUATION)
                        && isSimilarHeight(block, rect, nodeNodeSpacing)) {
                // Every rect is in its own block before comapction.
                block.addChild(rect);
            } else {
                // Case rect does not fit in block. Add new block to the right of it.
                Block newBlock = new Block(block.getX() + block.getWidth() + nodeNodeSpacing, row.getY(), row,
                        nodeNodeSpacing);
                row.addBlock(newBlock);
                newBlock.addChild(rect);
            }
            currentWidth = rect.getX() + rect.getWidth();
        }
        rows.add(row);
        return rows;
    }
    
    /**
     * Check whether a rectangle has a similar height as the other rectangles in the block and add it too it depending on
     * the available space.
     * @param row The row.
     * @param block The block.
     * @param rect The rectangle.
     * @param boundingWidth The bounding width of the row.
     * @param nodeNodeSpacing The spacing between two nodes.
     * @return true, if the rectangle was successfully added.
     */
    public static boolean placeRectInBlock(final RectRow row, final Block block, final ElkNode rect,
            final double boundingWidth, final double nodeNodeSpacing) {
        if (isSimilarHeight(block, rect, nodeNodeSpacing)) {
            if (block.getLastRowNewX() + rect.getWidth() + nodeNodeSpacing <= boundingWidth &&
                    (block.getLastRowY() - row.getY() + rect.getHeight() <= row.getHeight() || row.getChildren().size() == 1)) {
                // Case it fits in a row in the same block
                block.addChild(rect);
                return true;
            } else if (block.getX() + rect.getWidth() <= boundingWidth &&
                    (block.getY() + block.getHeight() + rect.getHeight() + nodeNodeSpacing <= row.getY() + row.getHeight())) {
                // Case a new row in the block can be opened
                block.addChildInNewRow(rect);
                return true;
            }
        }
        return false;
        
    }
    
    /**
     * Checks whether a rectangle has a similar height as the other rectangles in a block.
     * @param block The block.
     * @param rect The rectangle.
     * @param nodeNodeSpacing The spacing between two nodes.
     * @return true if the height of the rectangle is similar to the rectangles in the block.
     */
    public static boolean isSimilarHeight(final Block block, final ElkNode rect, final double nodeNodeSpacing) {
        if (rect.getHeight() >= block.getSmallestRectHeight() && rect.getHeight() <= block.getMinHeight()) {
            return true;
        } else {
            return block.getAverageHeight() * 0.5 <= rect.getHeight() &&
                    block.getAverageHeight() * 1.5 >= rect.getHeight();
        }
    }
}
