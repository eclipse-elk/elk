/*******************************************************************************
 * Copyright (c) 2026 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 * 
 * SPDX-License-Identifier: EPL-2.0 
 *******************************************************************************/
package org.eclipse.elk.alg.layered.intermediate;

import org.eclipse.elk.alg.common.nodespacing.NodeLabelAndSizeCalculator;
import org.eclipse.elk.alg.layered.LayeredPhases;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LGraphAdapters;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LNode.NodeType;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.core.alg.ILayoutPhase;
import org.eclipse.elk.core.alg.ILayoutProcessor;
import org.eclipse.elk.core.util.IElkProgressMonitor;
import org.eclipse.elk.core.util.adapters.GraphAdapters.GraphAdapter;

/**
 * Post-processor for the node flexibility option of network simplex. Fixes the node sizes according to the output of
 * node flexibility, then disables the option and computes node placement again with the set node placement algorithm.
 *
 */
public class NodePlacementRepeater implements ILayoutProcessor<LGraph> {

    /* (non-Javadoc)
     * @see org.eclipse.elk.core.alg.ILayoutProcessor#process(java.lang.Object, org.eclipse.elk.core.util.IElkProgressMonitor)
     */
    @Override
    public void process(LGraph graph, IElkProgressMonitor progressMonitor) {

        progressMonitor.begin("Node placement repeater", 1);

        ILayoutPhase<LayeredPhases, LGraph> nodePlacer = graph.getProperty(LayeredOptions.NODE_PLACEMENT_NETWORK_SIMPLEX_NODE_FLEXIBILITY_RECOMPUTE_NODE_PLACEMENT).create();
        // disable node flexibility for second node placement run
        graph.setProperty(LayeredOptions.NODE_PLACEMENT_NETWORK_SIMPLEX_NODE_FLEXIBILITY_DEFAULT, null);
        for (Layer layer : graph) {
            for (LNode node : layer.getNodes()) {
                // make sure node-level flexibility settings don't affect the second run
                node.setProperty(LayeredOptions.NODE_PLACEMENT_NETWORK_SIMPLEX_NODE_FLEXIBILITY, null);
                // default node size constraints == fixed size
                node.setProperty(LayeredOptions.NODE_SIZE_CONSTRAINTS, null);
                // reset node positions
                node.getPosition().y = 0;
            }
        }
        GraphAdapter<LGraph> adapterGraph = LGraphAdapters.adapt(
                graph,
                true,
                true,
                node -> node.getType() == NodeType.NORMAL);
        // recompute port positions based on now fixed node sizes
        NodeLabelAndSizeCalculator.process(adapterGraph);
        // second node placement run
        nodePlacer.process(graph, progressMonitor);

        progressMonitor.done();
    }
}
