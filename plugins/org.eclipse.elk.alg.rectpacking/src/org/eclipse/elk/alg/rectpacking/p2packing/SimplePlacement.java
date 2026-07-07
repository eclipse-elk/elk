/*******************************************************************************
 * Copyright (c) 2022 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 * 
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later 
 *******************************************************************************/
package org.eclipse.elk.alg.rectpacking.p2packing;

import java.util.List;

import org.eclipse.elk.alg.rectpacking.RectPackingLayoutPhases;
import org.eclipse.elk.alg.rectpacking.options.InternalProperties;
import org.eclipse.elk.alg.rectpacking.options.RectPackingOptions;
import org.eclipse.elk.alg.rectpacking.util.Block;
import org.eclipse.elk.alg.rectpacking.util.BlockStack;
import org.eclipse.elk.alg.rectpacking.util.DrawingUtil;
import org.eclipse.elk.alg.rectpacking.util.RectRow;
import org.eclipse.elk.core.alg.ILayoutPhase;
import org.eclipse.elk.core.alg.LayoutProcessorConfiguration;
import org.eclipse.elk.core.math.ElkPadding;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.util.IElkProgressMonitor;
import org.eclipse.elk.graph.ElkNode;

/**
 * Places the given rectangles next to each other to form rows.
 * 
 * <dl>
 *   <dt>Precondition:</dt>
 *     <dd>{@link InternalProperties#TARGET_WIDTH} is set.</dd>
 *   <dt>Postcondition:</dt>
 *     <dd>All child nodes have coordinates such that:</dd>
 *     <dd>the next node is either right of the current one</dd>
 *     <dd>or in a new row.</dd>
 *     <dd>{@link InternalProperties#DRAWING_WIDTH},</dd>
 *     <dd>{@link InternalProperties#DRAWING_HEIGHT},</dd>
 *     <dd>{@link InternalProperties#ROWS},</dd>
 *     <dd>and {@link InternalProperties#ADDITIONAL_HEIGHT} are set.</dd>
 * </dl>
 */
public class SimplePlacement implements ILayoutPhase<RectPackingLayoutPhases, ElkNode> {
    /*
     * (non-Javadoc)
     * 
     * @see org.eclipse.elk.core.alg.ILayoutProcessor#process(java.lang.Object,
     * org.eclipse.elk.core.util.IElkProgressMonitor)
     */
    @Override
    public void process(ElkNode graph, IElkProgressMonitor progressMonitor) {
        progressMonitor.begin("No Compaction", 1);
        double targetWidth = graph.getProperty(InternalProperties.TARGET_WIDTH);
        double nodeNodeSpacing = graph.getProperty(RectPackingOptions.SPACING_NODE_NODE);
        ElkPadding padding = graph.getProperty(RectPackingOptions.PADDING);
        // Reset coordinates potentially set by width approximation.
        DrawingUtil.resetCoordinates(graph.getChildren());

        // Initial placement for rectangles in blocks in each row.
        List<RectRow> rows = InitialPlacement.place(graph.getChildren(), targetWidth, nodeNodeSpacing);
        KVector size;
        List<ElkNode> rectangles = graph.getChildren();
        if (rows == null) {
            size = DrawingUtil.calculateDimensions(rectangles);
        } else {
            // Put every block in its own block stack.
            for (RectRow row : rows) {
                for (Block block : row.getChildren()) {
                    BlockStack stack = new BlockStack(block.getX(), block.getY(),
                            graph.getProperty(RectPackingOptions.SPACING_NODE_NODE));
                    stack.addBlock(block);
                    row.getStacks().add(stack);
                }
            }
            size = DrawingUtil.calculateDimensions(rows, nodeNodeSpacing);
        }
        double width = Math.max(size.x, graph.getProperty(InternalProperties.MIN_WIDTH) - padding.getHorizontal());
        double height = Math.max(size.y, graph.getProperty(InternalProperties.MIN_HEIGHT) - padding.getVertical());
        double additionalHeight = height - size.y;
        graph.setProperty(InternalProperties.ADDITIONAL_HEIGHT, additionalHeight);
        graph.setProperty(InternalProperties.DRAWING_WIDTH, width);
        graph.setProperty(InternalProperties.DRAWING_HEIGHT, height + additionalHeight);
        graph.setProperty(InternalProperties.ROWS, rows);
        progressMonitor.done();
    }

    /*
     * (non-Javadoc)
     * 
     * @see org.eclipse.elk.core.alg.ILayoutPhase#getLayoutProcessorConfiguration(java.lang.Object)
     */
    @Override
    public LayoutProcessorConfiguration<RectPackingLayoutPhases, ElkNode> getLayoutProcessorConfiguration(
            ElkNode graph) {
        return null;
    }
}
