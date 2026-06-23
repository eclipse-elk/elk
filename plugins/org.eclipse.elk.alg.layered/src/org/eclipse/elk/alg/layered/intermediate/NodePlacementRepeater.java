/*******************************************************************************
 * Copyright (c) 2026 sdo and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 * 
 * SPDX-License-Identifier: EPL-2.0 
 *******************************************************************************/
package org.eclipse.elk.alg.layered.intermediate;

import org.eclipse.elk.alg.common.nodespacing.internal.NodeContext;
import org.eclipse.elk.alg.common.nodespacing.internal.algorithm.CellSystemConfigurator;
import org.eclipse.elk.alg.common.nodespacing.internal.algorithm.HorizontalPortPlacementSizeCalculator;
import org.eclipse.elk.alg.common.nodespacing.internal.algorithm.InsidePortLabelCellCreator;
import org.eclipse.elk.alg.common.nodespacing.internal.algorithm.NodeLabelAndSizeUtilities;
import org.eclipse.elk.alg.common.nodespacing.internal.algorithm.NodeSizeCalculator;
import org.eclipse.elk.alg.common.nodespacing.internal.algorithm.PortContextCreator;
import org.eclipse.elk.alg.common.nodespacing.internal.algorithm.PortLabelPlacementCalculator;
import org.eclipse.elk.alg.common.nodespacing.internal.algorithm.PortPlacementCalculator;
import org.eclipse.elk.alg.common.nodespacing.internal.algorithm.VerticalPortPlacementSizeCalculator;
import org.eclipse.elk.alg.layered.DebugUtil;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LGraphAdapters;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LNode.NodeType;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.options.FixedAlignment;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.layered.p4nodes.bk.BKNodePlacer;
import org.eclipse.elk.core.alg.ILayoutProcessor;
import org.eclipse.elk.core.util.IElkProgressMonitor;
import org.eclipse.elk.core.util.adapters.GraphAdapters.GraphAdapter;
import org.eclipse.elk.core.util.adapters.GraphAdapters.NodeAdapter;

/**
 * @author sdo
 *
 */
public class NodePlacementRepeater implements ILayoutProcessor<LGraph> {

    /* (non-Javadoc)
     * @see org.eclipse.elk.core.alg.ILayoutProcessor#process(java.lang.Object, org.eclipse.elk.core.util.IElkProgressMonitor)
     */
    @Override
    public void process(LGraph graph, IElkProgressMonitor progressMonitor) {
//        NetworkSimplexPlacer test = new NetworkSimplexPlacer();
//        SimpleNodePlacer test = new SimpleNodePlacer();
        BKNodePlacer test = new BKNodePlacer();
        graph.setProperty(LayeredOptions.NODE_PLACEMENT_BK_FIXED_ALIGNMENT, FixedAlignment.BALANCED);
        graph.setProperty(LayeredOptions.NODE_PLACEMENT_NETWORK_SIMPLEX_NODE_FLEXIBILITY, null);
        graph.setProperty(LayeredOptions.NODE_PLACEMENT_NETWORK_SIMPLEX_NODE_FLEXIBILITY_DEFAULT, null);
        for (Layer layer : graph) {
            for (LNode node : layer.getNodes()) {
                node.setProperty(LayeredOptions.NODE_PLACEMENT_NETWORK_SIMPLEX_NODE_FLEXIBILITY, null);
                node.setProperty(LayeredOptions.NODE_PLACEMENT_NETWORK_SIMPLEX_NODE_FLEXIBILITY_DEFAULT, null);
                node.setProperty(LayeredOptions.NODE_SIZE_CONSTRAINTS, null);
                node.getPosition().y = 0;
//                node.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_ORDER);
                
//                NodeContext nodeContext = new NodeContext(graph, node);
//                PortContextCreator.createPortContexts(nodeContext, ignoreInsidePortLabels);
            }
        }
        GraphAdapter<LGraph> adapterGraph = LGraphAdapters.adapt(
                graph,
                true,
                true,
                node -> node.getType() == NodeType.NORMAL);
        for (NodeAdapter<?> node : adapterGraph.getNodes()) {
            NodeContext nodeContext = new NodeContext(adapterGraph, node);
            PortContextCreator.createPortContexts(nodeContext, true);
            
            InsidePortLabelCellCreator.createInsidePortLabelCells(nodeContext);
            
            NodeLabelAndSizeUtilities.setupMinimumClientAreaSize(nodeContext);
            NodeLabelAndSizeUtilities.setupNodePaddingForPortsWithOffset(nodeContext);
            
            HorizontalPortPlacementSizeCalculator.calculateHorizontalPortPlacementSize(nodeContext);
            VerticalPortPlacementSizeCalculator.calculateVerticalPortPlacementSize(nodeContext);
            
            CellSystemConfigurator.configureCellSystemSizeContributions(nodeContext);
            
            NodeSizeCalculator.setNodeWidth(nodeContext);
            
            PortPlacementCalculator.placeHorizontalPorts(nodeContext);
            PortLabelPlacementCalculator.placeHorizontalPortLabels(nodeContext);
            
            CellSystemConfigurator.updateVerticalInsidePortLabelCellPadding(nodeContext);
            
            NodeSizeCalculator.setNodeHeight(nodeContext);
            
            NodeLabelAndSizeUtilities.offsetSouthernPortsByNodeSize(nodeContext);

            PortPlacementCalculator.placeVerticalPorts(nodeContext);
            PortLabelPlacementCalculator.placeVerticalPortLabels(nodeContext);

            NodeLabelAndSizeUtilities.setNodePadding(nodeContext);
            NodeLabelAndSizeUtilities.applyStuff(nodeContext);
        }
//        graph.setProperty(InternalProperties.WAS_FLEXIBLE, true);

        // elkjs-exclude-start
        if (progressMonitor.isLoggingEnabled()) {
            DebugUtil.logDebugGraph(progressMonitor, graph, 0, "etst");
        }
        // elkjs-exclude-end
        test.process(graph, progressMonitor);

        // elkjs-exclude-start
        if (progressMonitor.isLoggingEnabled()) {
            DebugUtil.logDebugGraph(progressMonitor, graph, 0, "test1");
        }
        // elkjs-exclude-end

    }

}
