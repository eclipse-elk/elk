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

import org.eclipse.elk.alg.rectpacking.RectPackingLayoutPhases;
import org.eclipse.elk.alg.rectpacking.options.InternalProperties;
import org.eclipse.elk.alg.rectpacking.options.RectPackingOptions;
import org.eclipse.elk.core.UnsupportedConfigurationException;
import org.eclipse.elk.core.alg.ILayoutPhase;
import org.eclipse.elk.core.alg.LayoutProcessorConfiguration;
import org.eclipse.elk.core.util.IElkProgressMonitor;
import org.eclipse.elk.graph.ElkNode;

/**
 * Sets a {@link InternalProperties#TARGET_WIDTH} based on the given target width.
 * 
 * <dl>
 *   <dt>Precondition: {@link RectPackingOptions#WIDTH_APPROXIMATION_TARGET_WIDTH}</dt>
 *   <dt>Postcondition:</dt>
 *     <dd>{@link InternalProperties#TARGET_WIDTH} is set on the graph.</dd>
 * </dl>
 */
public class TargetWidthWidthApproximator implements ILayoutPhase<RectPackingLayoutPhases, ElkNode> {

    /* (non-Javadoc)
     * @see org.eclipse.elk.core.alg.ILayoutProcessor#process(java.lang.Object, org.eclipse.elk.core.util.IElkProgressMonitor)
     */
    @Override
    public void process(ElkNode graph, IElkProgressMonitor progressMonitor) {
        progressMonitor.begin("Target Width Setter", 1);
        if (graph.hasProperty(RectPackingOptions.WIDTH_APPROXIMATION_TARGET_WIDTH)) {
            graph.setProperty(InternalProperties.TARGET_WIDTH,
                    graph.getProperty(RectPackingOptions.WIDTH_APPROXIMATION_TARGET_WIDTH));
        } else {
            throw new UnsupportedConfigurationException(
                    "A target width has to be set if the TargetWidthWidthApproximator should be used.");
        }
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
