/*******************************************************************************
 * Copyright (c) 2010, 2020 Kiel University and others.
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
package org.eclipse.elk.alg.layered.intermediate;

import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LMargin;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LNode.NodeType;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.core.alg.ILayoutProcessor;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.util.IElkProgressMonitor;

/**
 * Calculates the size of each layer, the height of the graph, and the graph's vertical offset. The
 * graph's width cannot be set at this point since that would require the edge routing phase to have
 * already finished.
 * 
 * <dl>
 *   <dt>Precondition:</dt>
 *     <dd>A layered graph with finished node placement.</dd>
 *     <dd>Node margins are calculated.</dd>
 *   <dt>Postcondition:</dt>
 *     <dd>Layer sizes are set</dd>
 *     <dd>Graph's height is set</dd>
 *     <dd>Graph's offset includes node positions</dd>
 *   <dt>Slots:</dt>
 *     <dd>Before phase 5.</dd>
 *   <dt>Same-slot dependencies:</dt>
 *     <dd>None.</dd>
 * </dl>
 */
public final class LayerSizeAndGraphHeightCalculator implements ILayoutProcessor<LGraph> {

    @Override
    public void process(final LGraph layeredGraph, final IElkProgressMonitor monitor) {
        monitor.begin("Layer size calculation", 1);
        
        // Remember the lowest and the highest y coordinates
        double minY = Double.POSITIVE_INFINITY;
        double maxY = Double.NEGATIVE_INFINITY;
        
        boolean foundNodes = false;
        for (Layer layer : layeredGraph) {
            // Retrieve and reset layer size (just to be sure...)
            KVector layerSize = layer.getSize();
            layerSize.x = 0.0;
            layerSize.y = 0.0;
            
            // If the layer is empty, skip the rest (should not usually happen)
            if (layer.getNodes().isEmpty()) {
                continue;
            }
            
            foundNodes = true;
            
            // Calculate the layer's width
            for (LNode node : layer) {
                KVector nodeSize = node.getSize();
                LMargin nodeMargin = node.getMargin();
                
                layerSize.x = Math.max(
                        layerSize.x,
                        nodeSize.x + nodeMargin.left + nodeMargin.right);
            }
            
            // Calculate the layer's height
            LNode firstNode = layer.getNodes().get(0);
            double top = firstNode.getPosition().y - firstNode.getMargin().top;
            if (firstNode.getType() == NodeType.EXTERNAL_PORT) {
                // Consider port spacing options for first node
                top -= layeredGraph.getProperty(LayeredOptions.SPACING_PORTS_SURROUNDING).getTop();
            }
            LNode lastNode = layer.getNodes().get(layer.getNodes().size() - 1);
            double bottom = lastNode.getPosition().y + lastNode.getSize().y
                    + lastNode.getMargin().bottom;
            if (lastNode.getType() == NodeType.EXTERNAL_PORT) {
                // Consider port spacing options for last node
                bottom += layeredGraph.getProperty(LayeredOptions.SPACING_PORTS_SURROUNDING).getBottom();
            }
            layerSize.y = bottom - top;
            
            // Update the lowest and highest encountered Y coordinate
            minY = Math.min(minY, top);
            maxY = Math.max(maxY, bottom);
        }
        
        // If we didn't find any nodes in the graph, be sure to set the values properly. This can happen if a connected
        // component only contains external ports.
        if (!foundNodes) {
            minY = 0;
            maxY = 0;
        }

        // Set the graph's height
        layeredGraph.getSize().y = maxY - minY;
        layeredGraph.getOffset().y -= minY;
        
        monitor.done();
    }

}
