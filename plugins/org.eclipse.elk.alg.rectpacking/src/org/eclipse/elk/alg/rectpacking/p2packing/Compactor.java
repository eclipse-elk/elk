/*******************************************************************************
 * Copyright (c) 2022 Kiel University and others.
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

import org.eclipse.elk.alg.rectpacking.RectPackingLayoutPhases;
import org.eclipse.elk.alg.rectpacking.options.InternalProperties;
import org.eclipse.elk.alg.rectpacking.options.RectPackingOptions;
import org.eclipse.elk.alg.rectpacking.util.DrawingData;
import org.eclipse.elk.core.alg.ILayoutPhase;
import org.eclipse.elk.core.alg.LayoutProcessorConfiguration;
import org.eclipse.elk.core.math.ElkPadding;
import org.eclipse.elk.core.util.IElkProgressMonitor;
import org.eclipse.elk.graph.ElkGraphFactory;
import org.eclipse.elk.graph.ElkNode;

/**
 * Places and compact the given rectangles by forming rows of stacks of blocks of subrows to maintain a common
 * reading direction inside a target width.
 * 
 * <dl>
 *   <dt>Precondition:</dt>
 *     <dd>{@link InternalProperties#TARGET_WIDTH} is set.</dd>
 *   <dt>Postcondition:</dt>
 *     <dd>All child nodes have coordinates such that:</dd>
 *     <dd>the next node is either in the same subrows/row,</dd>
 *     <dd>a new subrows/row,</dd>
 *     <dd>or in a new stack in the same row.</dd>
 *     <dd>{@link InternalProperties#DRAWING_WIDTH},</dd>
 *     <dd>{@link InternalProperties#DRAWING_HEIGHT},</dd>
 *     <dd>{@link InternalProperties#ROWS},</dd>
 *     <dd>and {@link InternalProperties#ADDITIONAL_HEIGHT} are set.</dd>
 * </dl>
 */
public class Compactor implements ILayoutPhase<RectPackingLayoutPhases, ElkNode> {

    /* (non-Javadoc)
     * @see org.eclipse.elk.core.alg.ILayoutProcessor#process(java.lang.Object, org.eclipse.elk.core.util.IElkProgressMonitor)
     */
    @Override
    public void process(ElkNode graph, IElkProgressMonitor progressMonitor) {
        progressMonitor.begin("Compaction", 1);
        List<ElkNode> rectangles = graph.getChildren();
        double aspectRatio = graph.getProperty(RectPackingOptions.ASPECT_RATIO);
        double nodeNodeSpacing = graph.getProperty(RectPackingOptions.SPACING_NODE_NODE);
        ElkPadding padding = graph.getProperty(RectPackingOptions.PADDING);
        
        RowFillingAndCompaction secondIt = new RowFillingAndCompaction(aspectRatio, nodeNodeSpacing);
        DrawingData drawing = secondIt.start(rectangles, progressMonitor, graph, padding);
        // elkjs-exclude-start
        if (progressMonitor.isLoggingEnabled()) {
            progressMonitor.logGraph(graph, "Compacted");
        }
        // elkjs-exclude-end
        // Begin possible iterations to improve rectpacking by setting a new target width and repeating the compaction.
        copyRowWidthChangeValues(graph, secondIt);
        
        // Begin more compaction iterations if more than one iteration is specified.
        int iterations = graph.getProperty(RectPackingOptions.PACKING_COMPACTION_ITERATIONS);
        while (iterations > 1) {
            // Create a shallow clone based on properties and sizes of children (not grandchildren).
            ElkNode clone = clone(graph);
            double oldSM = drawing.getScaleMeasure();
            // Calculate new target width and configure clone.
            configureSecondIteration(graph, clone, drawing);
            // Run additional compaction step.
            secondIt = new RowFillingAndCompaction(aspectRatio, nodeNodeSpacing);
            DrawingData newDrawing = secondIt.start(rectangles, progressMonitor, clone, padding);

            // elkjs-exclude-start
            if (progressMonitor.isLoggingEnabled()) {
                progressMonitor.logGraph(clone, "Layouted clone " + iterations);
            }
            // elkjs-exclude-end
            // Compare scale measure and choose the best packing.
            double newSM = newDrawing.getScaleMeasure();

            if (newSM >= oldSM && newSM == (double) newSM) {
                // If the new packing is better apply packing to original graph.
                for (int i = 0; i < clone.getChildren().size(); i++) {
                    copyPosition(clone.getChildren().get(i), graph.getChildren().get(i));
                }
                copyRowWidthChangeValues(graph, secondIt);
                drawing.setDrawingWidth(newDrawing.getDrawingWidth());
                drawing.setDrawingHeight(newDrawing.getDrawingHeight());

            }
            iterations--;
        }
        
        graph.setProperty(InternalProperties.DRAWING_HEIGHT, drawing.getDrawingHeight());
        graph.setProperty(InternalProperties.DRAWING_WIDTH, drawing.getDrawingWidth());
        progressMonitor.done();
    }

    /**
     * Copies the row width increase and decrease to the graph.
     * 
     * @param graph The graph
     * @param compaction Compaction data
     */
    private void copyRowWidthChangeValues(ElkNode graph, RowFillingAndCompaction compaction) {
        graph.setProperty(InternalProperties.MIN_ROW_INCREASE, compaction.potentialRowWidthIncreaseMin);
        graph.setProperty(InternalProperties.MAX_ROW_INCREASE, compaction.potentialRowWidthIncreaseMax);
        graph.setProperty(InternalProperties.MIN_ROW_DECREASE, compaction.potentialRowWidthDecreaseMin);
        graph.setProperty(InternalProperties.MAX_ROW_DECREASE, compaction.potentialRowWidthDecreaseMax);
        
    }

    /* (non-Javadoc)
     * @see org.eclipse.elk.core.alg.ILayoutPhase#getLayoutProcessorConfiguration(java.lang.Object)
     */
    @Override
    public LayoutProcessorConfiguration<RectPackingLayoutPhases, ElkNode> getLayoutProcessorConfiguration(
            ElkNode graph) {
        return null;
    }
    
    /**
     * Set new target width on clone.
     * 
     * @param layoutGraph The original graph.
     * @param drawing The initial drawing.
     */
    private void configureSecondIteration(ElkNode layoutGraph, ElkNode clone, DrawingData drawing) {
        ElkPadding padding = layoutGraph.getProperty(RectPackingOptions.PADDING);
        double aspectRatio = layoutGraph.getProperty(RectPackingOptions.ASPECT_RATIO);
        // Try to layout again if the aspect ratio seems to be bad
        if (layoutGraph.getChildren().size() > 1
                && layoutGraph.getProperty(InternalProperties.MIN_ROW_INCREASE) != Double.POSITIVE_INFINITY
                && (drawing.getDrawingWidth() + padding.getHorizontal())
                        / (drawing.getDrawingHeight() + padding.getVertical()) < aspectRatio) {
            // The drawing is too high, this means the approximated target width is too low
            // The new target width will be set to the next higher value that would change something.
            clone.setProperty(InternalProperties.TARGET_WIDTH, layoutGraph.getProperty(InternalProperties.TARGET_WIDTH)
                    + layoutGraph.getProperty(InternalProperties.MIN_ROW_INCREASE));
        } else if (layoutGraph.getChildren().size() > 1
                && layoutGraph.getProperty(InternalProperties.MIN_ROW_DECREASE) != Double.POSITIVE_INFINITY
                && (drawing.getDrawingWidth() + padding.getHorizontal())
                        / (drawing.getDrawingHeight() + padding.getVertical()) > aspectRatio) {
            // The drawing is too high, this means the approximated target width is too high.
            // The new target width will be set to the next smaller value that would change something.
            clone.setProperty(InternalProperties.TARGET_WIDTH,
                    Math.max(layoutGraph.getProperty(InternalProperties.MIN_WIDTH),
                    clone.getProperty(InternalProperties.TARGET_WIDTH)
                    - layoutGraph.getProperty(InternalProperties.MIN_ROW_DECREASE)));
        }
    }
    
    /**
     * Clones a node including all properties and its children with their properties.
     * 
     * @param node The node to clone
     * @return The clones node
     */
    private ElkNode clone(ElkNode node) {
        ElkNode clone = ElkGraphFactory.eINSTANCE.createElkNode();
        clone.copyProperties(node);
        for (ElkNode child : node.getChildren()) {
            ElkNode newChild = ElkGraphFactory.eINSTANCE.createElkNode();
            newChild.setParent(clone);
            newChild.setDimensions(child.getWidth(), child.getHeight());
            newChild.setIdentifier(child.getIdentifier());
            newChild.setLocation(child.getX(), child.getY());
            clone.getChildren().add(newChild);
            newChild.copyProperties(child);
        }
        return clone;
    }

    /**
     * Copy the position of dimension of a node to the other and does the same for its children.
     * 
     * @param node The node to copy from.
     * @param other The node to copy to.
     */
    private void copyPosition(ElkNode node, ElkNode other) {
        other.setDimensions(node.getWidth(), node.getHeight());
        other.setLocation(node.getX(), node.getY());
        for (int i = 0; i < node.getChildren().size(); i++) {
            copyPosition(node.getChildren().get(i), other.getChildren().get(i));
        }
    }

}
