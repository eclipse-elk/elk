/*******************************************************************************
 * Copyright (c) 2012, 2020 Kiel University and others.
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

import static org.junit.Assert.assertEquals;

import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.test.framework.LayoutTestRunner;
import org.eclipse.elk.alg.test.framework.annotations.Algorithm;
import org.eclipse.elk.alg.test.framework.annotations.DefaultConfiguration;
import org.eclipse.elk.alg.test.framework.annotations.GraphProvider;
import org.eclipse.elk.alg.test.framework.annotations.TestAfterProcessor;
import org.eclipse.elk.core.options.PortConstraints;
import org.eclipse.elk.core.options.PortSide;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.ElkPort;
import org.eclipse.elk.graph.util.ElkGraphUtil;
import org.junit.runner.RunWith;

/**
 * Tests whether the port offset is correctly applied to ports of all sides on nodes with all port
 * constraints. The test graph always consists of a node with one port on each side whose position is
 * somewhat away from that side and whose offset is set to zero. The result should always be that all
 * ports touch the border of the node.
 */
@RunWith(LayoutTestRunner.class)
@Algorithm(LayeredOptions.ALGORITHM_ID)
@DefaultConfiguration()
public class LabelAndNodeSizeProcessorOffsetTest {

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Sources

    /**
     * Creates a graph containing a node for each port constraints configuration and four ports each.
     */
    @GraphProvider
    public ElkNode produceGraph() {
        ElkNode graph = ElkGraphUtil.createGraph();
        
        generateNode(graph, PortConstraints.FREE);
        generateNode(graph, PortConstraints.FIXED_SIDE);
        generateNode(graph, PortConstraints.FIXED_ORDER);
        generateNode(graph, PortConstraints.FIXED_RATIO);
        generateNode(graph, PortConstraints.FIXED_POS);
        
        return graph;
    }

    /**
     * Generates a node for the given graph with one port on each side and with its port constraints
     * set to the given value.
     */
    private void generateNode(final ElkNode graph, final PortConstraints portConstraints) {
        // We need to set coordinates and don't want to add constants for them all
        // CHECKSTYLEOFF MagicNumber
        
        // Create the node
        ElkNode node = ElkGraphUtil.createNode(graph);
        node.setProperty(LayeredOptions.PORT_CONSTRAINTS, portConstraints);
        node.setWidth(100);
        
        // Create ports
        ElkPort nPort = ElkGraphUtil.createPort(node);
        nPort.setProperty(LayeredOptions.PORT_BORDER_OFFSET, 0.0);
        nPort.setProperty(LayeredOptions.PORT_SIDE, PortSide.NORTH);
        nPort.setX(50);
        nPort.setY(-300);

        ElkPort sPort = ElkGraphUtil.createPort(node);
        sPort.setProperty(LayeredOptions.PORT_BORDER_OFFSET, 0.0);
        sPort.setProperty(LayeredOptions.PORT_SIDE, PortSide.SOUTH);
        sPort.setX(50);
        sPort.setY(300);

        ElkPort ePort = ElkGraphUtil.createPort(node);
        ePort.setProperty(LayeredOptions.PORT_BORDER_OFFSET, 0.0);
        ePort.setProperty(LayeredOptions.PORT_SIDE, PortSide.EAST);
        ePort.setX(300);
        ePort.setY(50);

        ElkPort wPort = ElkGraphUtil.createPort(node);
        wPort.setProperty(LayeredOptions.PORT_BORDER_OFFSET, 0.0);
        wPort.setProperty(LayeredOptions.PORT_SIDE, PortSide.WEST);
        wPort.setX(-300);
        wPort.setY(50);

        // CHECKSTYLEON MagicNumber
    }
    

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Tests
    
    /**
     * Tests if the ports of the given graph really touch the border of their node.
     * 
     * @param graph the graph to test.
     */
    @TestAfterProcessor(LabelAndNodeSizeProcessor.class)
    public void testPortPositions(final Object graph) {
        LGraph lGraph = (LGraph) graph;
        
        for (Layer layer : lGraph) {
            for (LNode node : layer) {
                for (LPort port : node.getPorts()) {
                    switch (port.getSide()) {
                    case NORTH:
                        assertEquals(
                                generateErrorMessage(port),
                                -port.getSize().y,
                                port.getPosition().y,
                                0.0);
                        break;
                        
                    case SOUTH:
                        assertEquals(
                                generateErrorMessage(port),
                                node.getSize().y,
                                port.getPosition().y,
                                0.0);
                        break;
                        
                    case EAST:
                        assertEquals(
                                generateErrorMessage(port),
                                node.getSize().x,
                                port.getPosition().x,
                                0.0);
                        break;
                        
                    case WEST:
                        assertEquals(
                                generateErrorMessage(port),
                                -port.getSize().x,
                                port.getPosition().x,
                                0.0);
                        break;
                    }
                }
            }
        }
    }
    
    /**
     * Generates a meaningful error message in case the given port doesn't touch its parent node's
     * border.
     * 
     * @param port the port.
     * @return a meaningful error message.
     */
    private String generateErrorMessage(final LPort port) {
        return "Port (" + port.getSide() + ") does not touch border of node with port constraints "
                + port.getNode().getProperty(LayeredOptions.PORT_CONSTRAINTS);
    }
    
}
