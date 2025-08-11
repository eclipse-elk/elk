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

import org.eclipse.elk.alg.common.nodespacing.NodeLabelAndSizeCalculator;
import org.eclipse.elk.alg.rectpacking.options.InternalProperties;
import org.eclipse.elk.core.alg.ILayoutProcessor;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.util.ElkUtil;
import org.eclipse.elk.core.util.IElkProgressMonitor;
import org.eclipse.elk.core.util.adapters.ElkGraphAdapters;
import org.eclipse.elk.core.util.adapters.GraphAdapters.GraphAdapter;
import org.eclipse.elk.core.util.adapters.GraphAdapters.NodeAdapter;
import org.eclipse.elk.graph.ElkNode;

/**
 * Sets the minimum width and height of the parent since this influences decisions based on the aspect ratio.
 * 
 * <dl>
 *   <dt>Precondition:</dt>
 *   <dt>Postcondition:</dt>
 *     <dd>The minimum width and height on the graph.</dd>
 *   <dt>Slots:</dt>
 *     <dd>Before phase 1.</dd>
 *   <dt>Same-slot dependencies:</dt>
 * </dl>
 */
public class MinSizePreProcessor implements ILayoutProcessor<ElkNode> {

    /* (non-Javadoc)
     * @see org.eclipse.elk.core.alg.ILayoutProcessor#process(java.lang.Object, org.eclipse.elk.core.util.IElkProgressMonitor)
     */
    @Override
    public void process(ElkNode graph, IElkProgressMonitor progressMonitor) {
        progressMonitor.begin("Min Size Preprocessing", 1);
        // Get minimum size based on children.
        KVector minSize = ElkUtil.effectiveMinSizeConstraintFor(graph);
        
        // Get minimum size based on labels.
        if (graph.getParent() != null) {
            // only possible if a parent exists.
            GraphAdapter<?> graphAdapter = ElkGraphAdapters.adapt(graph.getParent());
            NodeAdapter<?> nodeAdapter = ElkGraphAdapters.adaptSingleNode(graph);
            
            KVector minSize2 = NodeLabelAndSizeCalculator.process(graphAdapter, nodeAdapter, false, true);
            minSize.x = Math.max(minSize.x, minSize2.x);
            minSize.y = Math.max(minSize.y, minSize2.y);
        }
        graph.setProperty(InternalProperties.MIN_WIDTH, minSize.x);
        graph.setProperty(InternalProperties.MIN_HEIGHT, minSize.y);
        progressMonitor.done();
    }

}
