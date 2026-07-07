/*******************************************************************************
 * Copyright (c) 2023 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 * 
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later 
 *******************************************************************************/
package org.eclipse.elk.alg.topdown.layered.test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;

import java.util.EnumSet;

import org.eclipse.elk.alg.test.PlainJavaInitialization;
import org.eclipse.elk.core.LayoutConfigurator;
import org.eclipse.elk.core.RecursiveGraphLayoutEngine;
import org.eclipse.elk.core.UnsupportedGraphException;
import org.eclipse.elk.core.data.LayoutAlgorithmResolver;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.options.SizeConstraint;
import org.eclipse.elk.core.options.TopdownNodeTypes;
import org.eclipse.elk.core.util.BasicProgressMonitor;
import org.eclipse.elk.core.util.ElkUtil;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.ElkPort;
import org.eclipse.elk.graph.util.ElkGraphUtil;
import org.junit.Test;

/**
 * Test interaction between top-down layout and the node size calculation 
 * performed for port placement in layered layout.
 *
 */
public class PortPlacementCalculationTest {
    
    /**
     * Tests that during a bottom-up layout the port positions are calculated
     * according to the minimum size they require.
     */
    @Test
    public void testBottomUp() {
        PlainJavaInitialization.initializePlainJavaLayout();
        ElkNode graph = ElkGraphUtil.createGraph();
 
        final double portSpacing = 10.0;
        
        ElkNode node = ElkGraphUtil.createNode(graph);
        node.setProperty(CoreOptions.NODE_SIZE_CONSTRAINTS, EnumSet.of(SizeConstraint.PORT_LABELS, SizeConstraint.PORTS));
        node.setProperty(CoreOptions.SPACING_PORT_PORT, portSpacing);
        node.setDimensions(20, 60);
        
        // create three ports that would require less than 60 height
        ElkPort port1 = ElkGraphUtil.createPort(node);
        ElkPort port2 = ElkGraphUtil.createPort(node);
        ElkPort port3 = ElkGraphUtil.createPort(node);
        
        // prepare layout engine
        LayoutConfigurator config = new LayoutConfigurator();
        ElkUtil.applyVisitors(graph, config, new LayoutAlgorithmResolver());
        // call layout with layout engine
        try {
            new RecursiveGraphLayoutEngine().layout(graph, new BasicProgressMonitor());
        } catch (UnsupportedGraphException exception) {
            fail(exception.toString());
        }
        
        assertEquals(portSpacing, Math.abs(port1.getY() - port2.getY()), 0.0001);
        assertEquals(portSpacing, Math.abs(port2.getY() - port3.getY()), 0.0001);
    }
    
    /**
     * Tests that during a top-down layout the actually set node size is considered
     * beside the port placement requirements. If the node is already larger the 
     * ports will need to be placed further apart.
     */
    @Test
    public void testTopDown() {
        PlainJavaInitialization.initializePlainJavaLayout();
        ElkNode graph = ElkGraphUtil.createGraph();
        graph.setProperty(CoreOptions.TOPDOWN_LAYOUT, true);
        graph.setProperty(CoreOptions.TOPDOWN_NODE_TYPE, TopdownNodeTypes.ROOT_NODE);
        
        final double portSpacing = 10.0;
        final double fixedNodeHeight = 60.0;
        
        ElkNode node = ElkGraphUtil.createNode(graph);
        node.setProperty(CoreOptions.TOPDOWN_LAYOUT, true);
        node.setProperty(CoreOptions.TOPDOWN_NODE_TYPE, TopdownNodeTypes.HIERARCHICAL_NODE);
        node.setProperty(CoreOptions.NODE_SIZE_FIXED_GRAPH_SIZE, true);
        node.setProperty(CoreOptions.NODE_SIZE_CONSTRAINTS, EnumSet.of(SizeConstraint.PORT_LABELS, SizeConstraint.PORTS));
        node.setProperty(CoreOptions.SPACING_PORT_PORT, portSpacing);
        node.setDimensions(20, fixedNodeHeight);
        
        // create three ports that would require less than 60 height
        ElkPort port1 = ElkGraphUtil.createPort(node);
        ElkPort port2 = ElkGraphUtil.createPort(node);
        ElkPort port3 = ElkGraphUtil.createPort(node);
        
        // prepare layout engine
        LayoutConfigurator config = new LayoutConfigurator();
        ElkUtil.applyVisitors(graph, config, new LayoutAlgorithmResolver());
        // call layout with layout engine
        try {
            new RecursiveGraphLayoutEngine().layout(graph, new BasicProgressMonitor());
        } catch (UnsupportedGraphException exception) {
            fail(exception.toString());
        }
        
        // Since the node size is fixed, the port spacing should now be the height
        // divided by the number of ports plus one
        double expectedPortSpacing = fixedNodeHeight / 4;
        assertEquals(expectedPortSpacing, Math.abs(port1.getY() - port2.getY()), 0.0001);
        assertEquals(expectedPortSpacing, Math.abs(port2.getY() - port3.getY()), 0.0001);
    }

}
