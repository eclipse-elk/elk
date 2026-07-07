/*******************************************************************************
 * Copyright (c) 2010, 2015 Kiel University and others.
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

import java.util.List;

import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LNode.NodeType;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.core.alg.ILayoutProcessor;
import org.eclipse.elk.core.options.Alignment;
import org.eclipse.elk.core.options.PortSide;
import org.eclipse.elk.core.util.IElkProgressMonitor;

import com.google.common.collect.Lists;

/**
 * Sets the width of hierarchical port dummies and sets the layer alignment of North/South port dummies
 * to Center.
 * 
 * <p>To see why this is necessary, please refer to the processor's Wiki documentation.</p>
 * 
 * <dl>
 *   <dt>Precondition:</dt>
 *     <dd>A layered graph with finished node ordering</dd>
 *     <dd>Node order respects in-layer constraints.</dd>
 *   <dt>Postcondition:</dt>
 *     <dd>External port dummies are assigned a width.</dd>
 *     <dd>Layer alignment of North/South external port dummies is set to Center.</dd> 
 *   <dt>Slots:</dt>
 *     <dd>Before phase 4.</dd>
 *   <dt>Same-slot dependencies:</dt>
 *     <dd>None.</dd>
 * </dl>
 * 
 * @see HierarchicalPortConstraintProcessor
 * @see HierarchicalPortOrthogonalEdgeRouter
 * @see HierarchicalPortPositionProcessor
 * @author cds
 */
public final class HierarchicalPortDummySizeProcessor implements ILayoutProcessor<LGraph> {

    @Override
    public void process(final LGraph layeredGraph, final IElkProgressMonitor monitor) {
        monitor.begin("Hierarchical port dummy size processing", 1);

        List<LNode> northernDummies = Lists.newArrayList();
        List<LNode> southernDummies = Lists.newArrayList();
        
        // Calculate the width difference (this assumes CENTER node alignment)
        //  The idea behind this is that ports are stacked on top of each other at the center of the 
        //  node. By iteratively increasing their size by twice the edgeEdge spacing between layers, 
        //  vertical edge segments are spaced by that amount since the edge's anchor will either 
        //  be the very first or very last point of the port.
        double edgeSpacing = layeredGraph.getProperty(LayeredOptions.SPACING_EDGE_EDGE_BETWEEN_LAYERS);
        double delta = edgeSpacing * 2;
        
        // Iterate through the layers
        for (Layer layer : layeredGraph) {
            northernDummies.clear();
            southernDummies.clear();
            
            // Collect northern and southern hierarchical port dummies
            for (LNode node : layer) {
                if (node.getType() == NodeType.EXTERNAL_PORT) {
                    PortSide side = node.getProperty(InternalProperties.EXT_PORT_SIDE);
                    
                    if (side == PortSide.NORTH) {
                        northernDummies.add(node);
                    } else if (side == PortSide.SOUTH) {
                        southernDummies.add(node);
                    }
                }
            }
            
            // Set widths
            setWidths(northernDummies, true, delta);
            setWidths(southernDummies, false, delta);
        }
        
        monitor.done();
    }
    
    /**
     * Sets the widths of the given list of nodes and sets their layer alignment properly.
     * 
     * @param nodes the list of nodes.
     * @param topDown {@code true} if the nodes should widen with increasing index, {@code false}
     *                if it should be the other way round.
     * @param delta the width difference from one node to the next.
     */
    private void setWidths(final List<LNode> nodes, final boolean topDown, final double delta) {
        double currentWidth = 0.0;
        double step = delta;
        
        if (!topDown) {
            // Start with the widest node, decreasing node size
            currentWidth = delta * (nodes.size() - 1);
            step *= -1.0;
        }
        
        for (LNode node : nodes) {
            node.setProperty(LayeredOptions.ALIGNMENT, Alignment.CENTER);
            node.getSize().x = currentWidth;
            
            // Move eastern ports to the node's right border
            for (LPort port : node.getPorts(PortSide.EAST)) {
                port.getPosition().x = currentWidth;
            }
            
            currentWidth += step;
        }
    }
    
}
