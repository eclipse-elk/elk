/*******************************************************************************
 * Copyright (c) 2022 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 * 
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later 
 *******************************************************************************/
package org.eclipse.elk.alg.rectpacking.p1widthapproximation;

import java.util.List;

import org.eclipse.elk.alg.rectpacking.RectPackingLayoutPhases;
import org.eclipse.elk.alg.rectpacking.options.InternalProperties;
import org.eclipse.elk.alg.rectpacking.options.OptimizationGoal;
import org.eclipse.elk.alg.rectpacking.options.RectPackingOptions;
import org.eclipse.elk.alg.rectpacking.util.DrawingData;
import org.eclipse.elk.alg.rectpacking.util.DrawingUtil;
import org.eclipse.elk.core.alg.ILayoutPhase;
import org.eclipse.elk.core.alg.LayoutProcessorConfiguration;
import org.eclipse.elk.core.math.ElkPadding;
import org.eclipse.elk.core.util.IElkProgressMonitor;
import org.eclipse.elk.graph.ElkNode;

/**
 * Approximates the target width based on a greedy placement strategy, which can be configured to use the scale measure,
 * area, or aspect ratio as its {@link OptimizationGoal}, to transform the rectangle packing problem into a strip
 * packing problem.
 * 
 * <dl>
 *   <dt>Precondition:</dt>
 *   <dt>Postcondition:</dt>
 *     <dd>{@link InternalProperties#TARGET_WIDTH} is set on the graph.</dd>
 * </dl>
 */
public class GreedyWidthApproximator implements ILayoutPhase<RectPackingLayoutPhases, ElkNode> {

    /* (non-Javadoc)
     * @see org.eclipse.elk.core.alg.ILayoutProcessor#process(java.lang.Object, org.eclipse.elk.core.util.IElkProgressMonitor)
     */
    @Override
    public void process(ElkNode graph, IElkProgressMonitor progressMonitor) {
        progressMonitor.begin("Greedy Width Approximator", 1);
        // The desired aspect ratio.
        double aspectRatio = graph.getProperty(RectPackingOptions.ASPECT_RATIO);
        // The padding surrounding the drawing.
        ElkPadding padding = graph.getProperty(RectPackingOptions.PADDING);
        // The strategy for the initial width approximation.
        OptimizationGoal goal = graph.getProperty(RectPackingOptions.WIDTH_APPROXIMATION_OPTIMIZATION_GOAL);
        // Option for better width approximation.
        boolean lastPlaceShift = graph.getProperty(RectPackingOptions.WIDTH_APPROXIMATION_LAST_PLACE_SHIFT);
        //  The spacing between two nodes.
        double nodeNodeSpacing = graph.getProperty(RectPackingOptions.SPACING_NODE_NODE);
        
        List<ElkNode> rectangles = graph.getChildren();
        DrawingUtil.resetCoordinates(rectangles);
        DrawingData drawing;
        
        // Initial width approximation.
        AreaApproximation firstIt = new AreaApproximation(aspectRatio, goal, lastPlaceShift);
        drawing = firstIt.approxBoundingBox(rectangles, nodeNodeSpacing, padding);
              
        graph.setProperty(InternalProperties.TARGET_WIDTH, drawing.getDrawingWidth());
        progressMonitor.done();
    }

    /* (non-Javadoc)
     * @see org.eclipse.elk.core.alg.ILayoutPhase#getLayoutProcessorConfiguration(java.lang.Object)
     */
    @Override
    public LayoutProcessorConfiguration<RectPackingLayoutPhases, ElkNode> getLayoutProcessorConfiguration(
            ElkNode graph) {
        return null;
    }

}
