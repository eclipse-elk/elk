/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.layered.intermediate;

import java.util.List;

import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LNode.NodeType;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.core.alg.ILayoutProcessor;
import org.eclipse.elk.core.options.PortConstraints;
import org.eclipse.elk.core.options.PortSide;
import org.eclipse.elk.core.util.IElkProgressMonitor;

/**
 * Sets the y coordinate of external node dummies representing eastern or western hierarchical ports. Note that due to
 * additional space required to route edges connected to northern external ports, the y coordinate set here may become
 * invalid and may need to be fixed later. That fixing is part of what {@link HierarchicalPortOrthogonalEdgeRouter}
 * does.
 * 
 * <p>This processor is only necessary for node placers that do not respect the
 * {@link InternalProperties#PORT_RATIO_OR_POSITION} property themselves.</p>
 * 
 * <p>Note that this code involves a few subtleties with the graph's offset. The vertical offset is calculated by the
 * {@link LayerSizeAndGraphHeightCalculator}, which is executed before this processor has fixed the positions of the
 * external port dummies of eastern and western ports. However, those dummies were most likely not placed at a smaller
 * y coordinate than other nodes. This processor either moves dummies downwards, which has no influence on the graph's
 * offset, or upwards. If this had an influence on the offset, external ports could not be placed outside their border,
 * which is a perfectly valid configuration. We thus do not touch the offset again after moving the dummies around.</p>
 * 
 * <dl>
 *   <dt>Precondition:</dt>
 *     <dd>A layered graph with finished node placement.</dd>
 *     <dd>The graph's vertical coordinate offset is correct.</dd>
 *     <dd>External port dummies of eastern and western ports appear in their layer's node list in the order they
 *         will later appear on their node side in.</dd>
 *   <dt>Postcondition:</dt>
 *     <dd>External node dummies representing western or eastern ports have a correct y coordinate.</dd>
 *   <dt>Slots:</dt>
 *     <dd>Before phase 5.</dd>
 *   <dt>Same-slot dependencies:</dt>
 *     <dd>{@link LayerSizeAndGraphHeightCalculator}</dd>
 * </dl>
 * 
 * @see HierarchicalPortConstraintProcessor
 * @see HierarchicalPortDummySizeProcessor
 * @see HierarchicalPortOrthogonalEdgeRouter
 */
public final class HierarchicalPortPositionProcessor implements ILayoutProcessor<LGraph> {
    
    @Override
    public void process(final LGraph layeredGraph, final IElkProgressMonitor monitor) {
        monitor.begin("Hierarchical port position processing", 1);
        
        List<Layer> layers = layeredGraph.getLayers();
        
        // We're interested in EAST and WEST external port dummies only; since they can only be in
        // the first or last layer, only fix coordinates of nodes in those two layers
        if (layers.size() > 0) {
            fixCoordinates(layers.get(0), layeredGraph);
        }
        
        if (layers.size() > 1) {
            fixCoordinates(layers.get(layers.size() - 1), layeredGraph);
        }
        
        monitor.done();
    }
    
    /**
     * Fixes the y coordinates of external port dummies in the given layer.
     * 
     * @param layer the layer.
     * @param layeredGraph the layered graph.
     * @param portConstraints the port constraints that apply to external ports.
     * @param graphHeight height of the graph.
     */
    private void fixCoordinates(final Layer layer, final LGraph layeredGraph) {
        PortConstraints portConstraints = layeredGraph.getProperty(LayeredOptions.PORT_CONSTRAINTS);
        if (!(portConstraints.isRatioFixed() || portConstraints.isPosFixed())) {
            // If coordinates are free to be set, we're done
            return;
        }
        
        double graphHeight = layeredGraph.getActualSize().y;
        
        // Iterate over the layer's nodes
        for (LNode node : layer) {
            // We only care about external port dummies...
            if (node.getType() != NodeType.EXTERNAL_PORT) {
                continue;
            }
            
            // ...representing eastern or western ports.
            PortSide extPortSide = node.getProperty(InternalProperties.EXT_PORT_SIDE);
            if (extPortSide != PortSide.EAST && extPortSide != PortSide.WEST) {
                continue;
            }
            
            double finalYCoordinate = node.getProperty(InternalProperties.PORT_RATIO_OR_POSITION);
            
            if (portConstraints == PortConstraints.FIXED_RATIO) {
                // finalYCoordinate is a ratio that must be multiplied with the graph's height
                finalYCoordinate *= graphHeight;
            }

            // Apply the node's new Y coordinate
            node.getPosition().y = finalYCoordinate - node.getProperty(LayeredOptions.PORT_ANCHOR).y;
            node.borderToContentAreaCoordinates(false, true);
        }
    }
    
}
