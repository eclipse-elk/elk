/*******************************************************************************
 * Copyright (c) 2011, 2015 Kiel University and others.
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
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.core.alg.ILayoutProcessor;
import org.eclipse.elk.core.options.PortConstraints;
import org.eclipse.elk.core.options.PortSide;
import org.eclipse.elk.core.util.IElkProgressMonitor;

/**
 * Makes sure ports have at least a fixed side. If they don't, input ports are assigned
 * to the left and output ports to the right side. This processor can run either before the
 * first phase or before the third phase. In the first slot it assigns port sides before any
 * edges are reversed, hence edges that are reversed later are routed "around" the nodes
 * using inverted ports. This can only be handled correctly if the {@link InvertedPortProcessor}
 * is active. In the third slot, however, the port sides are assigned <em>after</em> edges
 * are reversed, so no inverted ports will occur. This behavior is controlled by the option
 * {@link org.eclipse.elk.alg.layered.options.LayeredOptions#FEEDBACK_EDGES FEEDBACK_EDGES}.
 * 
 * <dl>
 *   <dt>Precondition:</dt>
 *     <dd>a layered graph.</dd>
 *   <dt>Postcondition:</dt>
 *     <dd>all nodes have their ports distributed, with port constraints set to fixed sides at the least.</dd>
 *   <dt>Slots:</dt>
 *     <dd>Before phase 1 or before phase 3.</dd>
 *   <dt>Same-slot dependencies:</dt>
 *     <dd>None.</dd>
 * </dl>
 */
public final class PortSideProcessor implements ILayoutProcessor<LGraph> {

    @Override
    public void process(final LGraph layeredGraph, final IElkProgressMonitor monitor) {
        monitor.begin("Port side processing", 1);
        
        // IF USED BEFORE PHASE 1
        // Iterate through the unlayered nodes of the graph
        for (LNode node : layeredGraph.getLayerlessNodes()) {
            process(node);
        }
        
        // IF USED BEFORE PHASE 3
        // Iterate through the nodes of all layers
        for (Layer layer : layeredGraph) {
            for (LNode node : layer) {
                process(node);
            }
        }
        
        monitor.done();
    }
    
    /**
     * Process the ports of the given node.
     * 
     * @param node a node
     */
    private void process(final LNode node) {
        if (node.getProperty(LayeredOptions.PORT_CONSTRAINTS).isSideFixed()) {
            // Check whether there are ports with undefined side
            for (LPort port : node.getPorts()) {
                if (port.getSide() == PortSide.UNDEFINED) {
                    setPortSide(port);
                }
            }
        } else {
            // Distribute all ports and change constraints to fixed side
            for (LPort port : node.getPorts()) {
                setPortSide(port);
            }
            node.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        }
    }
    
    /**
     * Places input ports on the left side and output ports on the right side of nodes.
     * 
     * @param port the port to set side and anchor position
     */
    public static void setPortSide(final LPort port) {
        LNode portDummy = port.getProperty(InternalProperties.PORT_DUMMY);
        if (portDummy != null) {
            port.setSide(portDummy.getProperty(InternalProperties.EXT_PORT_SIDE));
        } else if (port.getNetFlow() < 0) {
            port.setSide(PortSide.EAST);
        } else {
            port.setSide(PortSide.WEST);
        }
    }

}
