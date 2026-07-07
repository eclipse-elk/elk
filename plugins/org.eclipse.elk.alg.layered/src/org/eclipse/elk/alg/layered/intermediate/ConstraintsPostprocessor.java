/*******************************************************************************
 * Copyright (c) 2019, 2020 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.layered.intermediate;

import org.eclipse.elk.alg.layered.DebugUtil;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LNode.NodeType;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.core.alg.ILayoutProcessor;
import org.eclipse.elk.core.util.IElkProgressMonitor;

/**
 * Adds to each LNode the layerID and positionID that has been computed by ELK Layered.
 * <dl>
 *   <dt>Precondition:</dt>
 *      <dd>none</dd>
 *   <dt>Postcondition:</dt>
 *      <dd>Nodes have a layer and position id based on the layout.</dd>
 *   <dt>Slots:</dt>
 *      <dd>After phase 5.</dd>
 * </dl>
 */
public final class ConstraintsPostprocessor implements ILayoutProcessor<LGraph> {

    @Override
    public void process(final LGraph graph, final IElkProgressMonitor progressMonitor) {

        progressMonitor.begin("Constraints Postprocessor", 1);
        // elkjs-exclude-start
        if (progressMonitor.isLoggingEnabled()) {
            DebugUtil.logDebugGraph(progressMonitor, graph, 0, "ConstraintPostProcessor");
        }
        // elkjs-exclude-end
        // Uses for each loops since they should be more efficient because LinkedLists are accessed.
        // It still creates counter variables so that no getter calls on layer and node are needed.
        int layerIndex = 0;

        for (Layer layer : graph.getLayers()) {
            int posIndex = 0;
            
            boolean nodeLayer = false;
            for (LNode currentNode : layer.getNodes()) {
                if (currentNode.getType() == NodeType.NORMAL) {
                    nodeLayer = true;
                    currentNode.setProperty(LayeredOptions.LAYERING_LAYER_ID, layerIndex);
                    currentNode.setProperty(LayeredOptions.CROSSING_MINIMIZATION_POSITION_ID, posIndex);
                    posIndex++;
                }
            }
            // layers with no nodes in it should not increase the layer id
            if (nodeLayer) {
                layerIndex++;
            }
        }
        // elkjs-exclude-start
        if (progressMonitor.isLoggingEnabled()) {
            DebugUtil.logDebugGraph(progressMonitor, graph, 0, "ConstraintPostProcessor");
        }
        // elkjs-exclude-end
        progressMonitor.done();

    }

}
