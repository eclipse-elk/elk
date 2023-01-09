/*******************************************************************************
 * Copyright (c) 2022 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 * 
 * SPDX-License-Identifier: EPL-2.0 
 *******************************************************************************/
package org.eclipse.elk.alg.rectpacking.intermediate;

import org.eclipse.elk.alg.rectpacking.options.InternalProperties;
import org.eclipse.elk.core.alg.ILayoutProcessor;
import org.eclipse.elk.core.util.IElkProgressMonitor;
import org.eclipse.elk.graph.ElkNode;

/**
 * Sorts all child nodes by their desired position.
 * 
 * <dl>
 *   <dt>Precondition:</dt>
 *   <dt>Postcondition:</dt>
 *     <dd>Sets the {@link InternalProperties#TARGET_WIDTH} to a size at least as high as the minimum width.</dd>
 *   <dt>Slots:</dt>
 *     <dd>After phase 1.</dd>
 *   <dt>Same-slot dependencies:</dt>
 * </dl>
 */
public class MinSizePostProcessor implements ILayoutProcessor<ElkNode> {

    /* (non-Javadoc)
     * @see org.eclipse.elk.core.alg.ILayoutProcessor#process(java.lang.Object, org.eclipse.elk.core.util.IElkProgressMonitor)
     */
    @Override
    public void process(ElkNode graph, IElkProgressMonitor progressMonitor) {
        progressMonitor.begin("Min Size Postprocessing", 1);
        // Remove padding to get the space the algorithm can use.
        graph.setProperty(InternalProperties.TARGET_WIDTH, Math.max(graph.getProperty(InternalProperties.TARGET_WIDTH),
                graph.getProperty(InternalProperties.MIN_WIDTH)));
        progressMonitor.done();
    }

}
