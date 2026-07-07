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
import org.eclipse.elk.alg.rectpacking.util.DrawingUtil;
import org.eclipse.elk.core.alg.ILayoutPhase;
import org.eclipse.elk.core.alg.LayoutProcessorConfiguration;
import org.eclipse.elk.core.math.ElkPadding;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.util.IElkProgressMonitor;
import org.eclipse.elk.graph.ElkNode;

/**
 * Strategy that just does nothing with the input and just calculates the drawing size based on the coordinates of the
 * children. Whitespace elimination cannot be done with this packing strategy since it might not succeed.
 * 
 * <dl>
 *   <dt>Precondition:</dt>
 *     <dd>{@link InternalProperties#TARGET_WIDTH} is set.</dd>
 *     <dd>Child nodes already have sensible coordinates.</dd>
 *   <dt>Postcondition:</dt>
 *     <dd>{@link InternalProperties#DRAWING_WIDTH},</dd>
 *     <dd>{@link InternalProperties#DRAWING_HEIGHT},</dd>
 *     <dd>and {@link InternalProperties#ADDITIONAL_HEIGHT} are set.</dd>
 * </dl>
 */
public class NoPlacement implements ILayoutPhase<RectPackingLayoutPhases, ElkNode> {
    /*
     * (non-Javadoc)
     * 
     * @see org.eclipse.elk.core.alg.ILayoutProcessor#process(java.lang.Object,
     * org.eclipse.elk.core.util.IElkProgressMonitor)
     */
    @Override
    public void process(ElkNode graph, IElkProgressMonitor progressMonitor) {
        progressMonitor.begin("No Compaction", 1);
        ElkPadding padding = graph.getProperty(RectPackingOptions.PADDING);
        List<ElkNode> rectangles = graph.getChildren();
        
        KVector size = DrawingUtil.calculateDimensions(rectangles);
        double width = Math.max(size.x, graph.getProperty(InternalProperties.MIN_WIDTH) - padding.getHorizontal());
        double height = Math.max(size.y, graph.getProperty(InternalProperties.MIN_HEIGHT) - padding.getVertical());
        double additionalHeight = height - size.y;
        graph.setProperty(InternalProperties.ADDITIONAL_HEIGHT, additionalHeight);
        graph.setProperty(InternalProperties.DRAWING_WIDTH, width);
        graph.setProperty(InternalProperties.DRAWING_HEIGHT, height + additionalHeight);
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
