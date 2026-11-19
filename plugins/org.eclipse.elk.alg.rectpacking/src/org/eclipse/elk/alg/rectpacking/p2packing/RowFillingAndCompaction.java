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
package org.eclipse.elk.alg.rectpacking.p2packing;

import java.util.List;

import org.eclipse.elk.alg.rectpacking.options.InternalProperties;
import org.eclipse.elk.alg.rectpacking.options.RectPackingOptions;
import org.eclipse.elk.alg.rectpacking.util.Block;
import org.eclipse.elk.alg.rectpacking.util.BlockRow;
import org.eclipse.elk.alg.rectpacking.util.BlockStack;
import org.eclipse.elk.alg.rectpacking.util.DrawingData;
import org.eclipse.elk.alg.rectpacking.util.DrawingDataDescriptor;
import org.eclipse.elk.alg.rectpacking.util.DrawingUtil;
import org.eclipse.elk.alg.rectpacking.util.RectRow;
import org.eclipse.elk.core.math.ElkPadding;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.util.IElkProgressMonitor;
import org.eclipse.elk.core.util.Pair;
import org.eclipse.elk.graph.ElkNode;

/**
 * Second iteration of the algorithm. Actual placement of the boxes inside the approximated bounding box. Rectangles are
 * placed in rows, which are compacted and then filled with elements from the subsequent row.
 */
public class RowFillingAndCompaction {
    //////////////////////////////////////////////////////////////////
    // Fields.
    /** Desired aspect ratio. */
    private double aspectRatio;
    /** Spacing between two nodes. */
    private double nodeNodeSpacing;

    /**
     * The minimum width that has to be removed from the target width not allow the last element of a row in it.
     */
    public double potentialRowWidthDecreaseMin = Double.POSITIVE_INFINITY;
    /**
     * The maximum width that has to be removed from the target width not allow the last element of a row in it.
     */
    public double potentialRowWidthDecreaseMax = 0;
    /**
     * The minimum width that has to be added to the target width to allow an element from the next row to be placed
     * the previous one.
     */
    public double potentialRowWidthIncreaseMin = Double.POSITIVE_INFINITY;
    /**
     * The maximum width that has to be added to the target width to allow an element from the next row to be placed
     * the previous one.
     */
    public double potentialRowWidthIncreaseMax = 0;

    //////////////////////////////////////////////////////////////////
    // Constructors.
    /**
     * Creates an {@link RowFillingAndCompaction} object to execute the second iteration on.
     * 
     * @param aspectRatio The desired aspect ratio.
     * @param nodeNodeSpacing The spacing between two nodes.
     */
    public RowFillingAndCompaction(final double aspectRatio, final double nodeNodeSpacing) {
        this.aspectRatio = aspectRatio;
        this.nodeNodeSpacing = nodeNodeSpacing;
    }

    //////////////////////////////////////////////////////////////////
    // Starting method.
    /**
     * Placement of the rectangles given by {@link ElkNode} inside the given bounding box.
     * @param rectangles The set of rectangles to be placed inside the bounding box.
     * @param minParentSize The minimal size of the parent.
     * @param nodeNodeSpacing The spacing between two nodes.
     * @return Drawing data for a produced drawing.
     */
    public DrawingData start(final List<ElkNode> rectangles,
            final IElkProgressMonitor progressMonitor, final ElkNode layoutGraph, final ElkPadding padding) {
        double targetWidth = layoutGraph.getProperty(InternalProperties.TARGET_WIDTH);
        double minWidth = layoutGraph.getProperty(InternalProperties.MIN_WIDTH);
        double minHeight = layoutGraph.getProperty(InternalProperties.MIN_HEIGHT);
        // Reset coordinates potentially set by width approximation.
        DrawingUtil.resetCoordinates(layoutGraph.getChildren());

        // Initial placement for rectangles in blocks in each row.
        List<RectRow> rows = InitialPlacement.place(layoutGraph.getChildren(), targetWidth, nodeNodeSpacing);
        
        // Compaction of blocks.
        for (int rowIdx = 0; rowIdx < rows.size(); rowIdx++) {
            RectRow currentRow = rows.get(rowIdx);
            if (rowIdx != 0) {
                RectRow previousRow = rows.get(rowIdx - 1);
                currentRow.setY(previousRow.getY() + previousRow.getHeight() + nodeNodeSpacing);
            }
            Pair<Boolean, Boolean> result = Compaction.compact(rowIdx, rows, targetWidth, nodeNodeSpacing,
                    layoutGraph.getProperty(RectPackingOptions.PACKING_COMPACTION_ROW_HEIGHT_REEVALUATION));
            if (result.getSecond()) {
                // Reset the row such that stacks are removed, blocks are not fixed.
                for (Block block : currentRow.getChildren()) {
                    block.setFixed(false);
                    block.setPositionFixed(false);
                    // Reset precalculated min/max width/heights.
                    block.resetBlock();
                }
                currentRow.resetStacks();
                currentRow.setWidth(targetWidth);
                rowIdx--;
            } else {
                adjustWidthAndHeight(currentRow);
                // Check how much space would be needed in the current row to add the first block from the next one.
                // And how much space could be removed from each row by removing the last block.
                if (rowIdx + 1 < rows.size()) {

                    potentialRowWidthIncreaseMax =
                            Math.max(currentRow.getWidth() + nodeNodeSpacing
                                            + rows.get(rowIdx + 1).getFirstBlock().getWidth() - targetWidth,
                                    potentialRowWidthDecreaseMax);
                    potentialRowWidthIncreaseMin =
                            Math.min(currentRow.getWidth() + nodeNodeSpacing
                                            + rows.get(rowIdx + 1).getFirstBlock().getWidth() - targetWidth,
                                    potentialRowWidthDecreaseMin);
                    if (currentRow.getStacks().size() != 0) {
                        potentialRowWidthDecreaseMax = Math.max(potentialRowWidthDecreaseMax,
                                currentRow.getStacks().get(currentRow.getStacks().size() - 1).getWidth()
                                    + (currentRow.getStacks().size() <= 1 ? 0 : nodeNodeSpacing));
                        potentialRowWidthDecreaseMin = Math.min(potentialRowWidthDecreaseMax,
                                currentRow.getStacks().get(currentRow.getStacks().size() - 1).getWidth()
                                    + (currentRow.getStacks().size() <= 1 ? 0 : nodeNodeSpacing));
                    }
                }
                // Special case the graph has only one row with one block with several subrows
                if (rows.size() == 1) {
                    BlockStack lastStack = currentRow.getStacks().get(currentRow.getStacks().size() - 1);
                    Block lastBlock = lastStack.getBlocks().get(lastStack.getBlocks().size() -1);
                    for (BlockRow blockRow : lastBlock.getRows()) {
                        potentialRowWidthDecreaseMax =
                                Math.max(potentialRowWidthDecreaseMax, lastBlock.getWidth() - blockRow.getWidth());
                        potentialRowWidthDecreaseMin =
                                Math.min(potentialRowWidthDecreaseMin, lastBlock.getWidth() - blockRow.getWidth());
                        potentialRowWidthIncreaseMax =
                                Math.max(potentialRowWidthIncreaseMax, blockRow.getWidth() + nodeNodeSpacing);
                        potentialRowWidthIncreaseMin =
                                Math.min(potentialRowWidthIncreaseMin, blockRow.getWidth() + nodeNodeSpacing);
                    }
                }
                // Log graph after first row compaction.
                // elkjs-exclude-start
                if (progressMonitor.isLoggingEnabled()) {
                    progressMonitor.logGraph(layoutGraph, "Compacted row " + rowIdx);
                }
                // elkjs-exclude-end
            }
        }
         
        KVector size = DrawingUtil.calculateDimensions(rows, nodeNodeSpacing);
        
        double totalWidth = Math.max(size.x, minWidth - padding.getHorizontal());
        double height = Math.max(size.y, minHeight - padding.getVertical());
        double additionalHeight = height - size.y;
        layoutGraph.setProperty(InternalProperties.ADDITIONAL_HEIGHT, additionalHeight);
        layoutGraph.setProperty(InternalProperties.ROWS, rows);

        return new DrawingData(this.aspectRatio, totalWidth, size.y + additionalHeight, DrawingDataDescriptor.WHOLE_DRAWING);
    }

    //////////////////////////////////////////////////////////////////
    // Helper method.

    /**
     * Recalculates the width and height of the row.
     * @param row The row.
     */
    private void adjustWidthAndHeight(final RectRow row) {
        double maxHeight = 0;
        double maxWidth = 0;
        int index = 0;
        for (BlockStack stack : row.getStacks()) {
            stack.updateDimension();
            maxHeight = Math.max(maxHeight, stack.getHeight());
            maxWidth += stack.getWidth() + (index > 0 ? nodeNodeSpacing : 0);
            index++;
        }
        row.setHeight(maxHeight);
        row.setWidth(maxWidth);
    }
}
