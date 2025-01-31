/*******************************************************************************
 * Copyright (c) 2022-2024 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 * 
 * SPDX-License-Identifier: EPL-2.0 
 *******************************************************************************/
package org.eclipse.elk.alg.topdown.test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;

import org.eclipse.elk.alg.test.PlainJavaInitialization;
import org.eclipse.elk.core.LayoutConfigurator;
import org.eclipse.elk.core.RecursiveGraphLayoutEngine;
import org.eclipse.elk.core.UnsupportedGraphException;
import org.eclipse.elk.core.data.LayoutAlgorithmResolver;
import org.eclipse.elk.core.math.ElkPadding;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.options.TopdownNodeTypes;
import org.eclipse.elk.core.util.BasicProgressMonitor;
import org.eclipse.elk.core.util.ElkUtil;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.util.ElkGraphUtil;
import org.junit.Test;

/**
 * Tests for checking that topdown layout works correctly in different configurations.
 */
public class TopdownLayoutTest {
    
    // CHECKSTYLEOFF MagicNumber
    
    /**
     * Tests a topdown layout with two hierarchical levels. All nodes are hierarchical nodes. Size of the
     * toplevel node causes a horizontal scaling.
     */
    @Test
    public void testTwoLevelLayoutHorizontalScaling() {
        PlainJavaInitialization.initializePlainJavaLayout();
        ElkNode graph = ElkGraphUtil.createGraph();
        graph.setProperty(CoreOptions.TOPDOWN_LAYOUT, true);
        graph.setProperty(CoreOptions.TOPDOWN_NODE_TYPE, TopdownNodeTypes.ROOT_NODE);
        
        ElkNode toplevel = ElkGraphUtil.createNode(graph);
        toplevel.setProperty(CoreOptions.TOPDOWN_LAYOUT, true);
        toplevel.setProperty(CoreOptions.TOPDOWN_NODE_TYPE, TopdownNodeTypes.HIERARCHICAL_NODE);
        toplevel.setProperty(CoreOptions.NODE_SIZE_FIXED_GRAPH_SIZE, true);
        toplevel.setProperty(CoreOptions.ALGORITHM, "org.eclipse.elk.fixed");
        toplevel.setProperty(CoreOptions.PADDING, new ElkPadding());
        toplevel.setProperty(CoreOptions.SPACING_NODE_NODE, 0.0);
        // set size of node so that children will be scaled down 
        toplevel.setProperty(CoreOptions.TOPDOWN_HIERARCHICAL_NODE_WIDTH, 20.0);
        toplevel.setProperty(CoreOptions.TOPDOWN_HIERARCHICAL_NODE_ASPECT_RATIO, 0.4);
        
        ElkNode child1 = ElkGraphUtil.createNode(toplevel);
        child1.setProperty(CoreOptions.TOPDOWN_LAYOUT, true);
        child1.setProperty(CoreOptions.TOPDOWN_NODE_TYPE, TopdownNodeTypes.HIERARCHICAL_NODE);
        child1.setProperty(CoreOptions.NODE_SIZE_FIXED_GRAPH_SIZE, true);
        child1.setX(0);
        child1.setY(0);
        child1.setDimensions(30, 30);
        
        ElkNode child2 = ElkGraphUtil.createNode(toplevel);
        child2.setProperty(CoreOptions.TOPDOWN_LAYOUT, true);
        child2.setProperty(CoreOptions.TOPDOWN_NODE_TYPE, TopdownNodeTypes.HIERARCHICAL_NODE);
        child2.setProperty(CoreOptions.NODE_SIZE_FIXED_GRAPH_SIZE, true);
        child2.setX(0);
        child2.setY(40);
        child2.setDimensions(30, 30);
        
        
        // prepare layout engine
        LayoutConfigurator config = new LayoutConfigurator();
        ElkUtil.applyVisitors(graph, config, new LayoutAlgorithmResolver());
        // call layout with layout engine
        try {
            new RecursiveGraphLayoutEngine().layout(graph, new BasicProgressMonitor());
        } catch (UnsupportedGraphException exception) {
            fail(exception.toString());
        }
               
        // topdown scale factor of toplevel node
        assertEquals(20.0 / 30.0, toplevel.getProperty(CoreOptions.TOPDOWN_SCALE_FACTOR), 0.00001);
        
    }
    
    /**
     * Tests a topdown layout with two hierarchical levels. All nodes are hierarchical nodes. Size of the
     * toplevel node causes a vertical scaling.
     */
    @Test
    public void testTwoLevelLayoutVerticalScaling() {
        PlainJavaInitialization.initializePlainJavaLayout();
        ElkNode graph = ElkGraphUtil.createGraph();
        graph.setProperty(CoreOptions.TOPDOWN_LAYOUT, true);
        graph.setProperty(CoreOptions.TOPDOWN_NODE_TYPE, TopdownNodeTypes.ROOT_NODE);
        
        ElkNode toplevel = ElkGraphUtil.createNode(graph);
        toplevel.setProperty(CoreOptions.TOPDOWN_LAYOUT, true);
        toplevel.setProperty(CoreOptions.TOPDOWN_NODE_TYPE, TopdownNodeTypes.HIERARCHICAL_NODE);
        toplevel.setProperty(CoreOptions.NODE_SIZE_FIXED_GRAPH_SIZE, true);
        toplevel.setProperty(CoreOptions.ALGORITHM, "org.eclipse.elk.fixed");
        toplevel.setProperty(CoreOptions.PADDING, new ElkPadding());
        toplevel.setProperty(CoreOptions.SPACING_NODE_NODE, 0.0);
        // set size of node so that children will be scaled down
        toplevel.setProperty(CoreOptions.TOPDOWN_HIERARCHICAL_NODE_WIDTH, 40.0);
        toplevel.setProperty(CoreOptions.TOPDOWN_HIERARCHICAL_NODE_ASPECT_RATIO, 1.33333);
        
        ElkNode child1 = ElkGraphUtil.createNode(toplevel);
        child1.setProperty(CoreOptions.TOPDOWN_LAYOUT, true);
        child1.setProperty(CoreOptions.TOPDOWN_NODE_TYPE, TopdownNodeTypes.HIERARCHICAL_NODE);
        child1.setProperty(CoreOptions.NODE_SIZE_FIXED_GRAPH_SIZE, true);
        child1.setX(0);
        child1.setY(0);
        child1.setDimensions(30, 30);
        
        ElkNode child2 = ElkGraphUtil.createNode(toplevel);
        child2.setProperty(CoreOptions.TOPDOWN_LAYOUT, true);
        child2.setProperty(CoreOptions.TOPDOWN_NODE_TYPE, TopdownNodeTypes.HIERARCHICAL_NODE);
        child2.setProperty(CoreOptions.NODE_SIZE_FIXED_GRAPH_SIZE, true);
        child2.setX(0);
        child2.setY(40);
        child2.setDimensions(30, 30);
        
        
        // prepare layout engine
        LayoutConfigurator config = new LayoutConfigurator();
        ElkUtil.applyVisitors(graph, config, new LayoutAlgorithmResolver());
        // call layout with layout engine
        try {
            new RecursiveGraphLayoutEngine().layout(graph, new BasicProgressMonitor());
        } catch (UnsupportedGraphException exception) {
            fail(exception.toString());
        }
               
        // topdown scale factor of toplevel node
        assertEquals(30.0 / 70.0, toplevel.getProperty(CoreOptions.TOPDOWN_SCALE_FACTOR), 0.00001);
        
    }
    
    /**
     * Test scale cap property when the result scale is larger than default, but smaller than
     * the set scale cap.
     */
    @Test
    public void testScaleCap() {
        PlainJavaInitialization.initializePlainJavaLayout();
        ElkNode graph = ElkGraphUtil.createGraph();
        graph.setProperty(CoreOptions.TOPDOWN_LAYOUT, true);
        graph.setProperty(CoreOptions.TOPDOWN_NODE_TYPE, TopdownNodeTypes.ROOT_NODE);
        
        ElkNode toplevel = ElkGraphUtil.createNode(graph);
        toplevel.setProperty(CoreOptions.TOPDOWN_LAYOUT, true);
        toplevel.setProperty(CoreOptions.TOPDOWN_NODE_TYPE, TopdownNodeTypes.HIERARCHICAL_NODE);
        toplevel.setProperty(CoreOptions.NODE_SIZE_FIXED_GRAPH_SIZE, true);
        toplevel.setProperty(CoreOptions.TOPDOWN_SCALE_CAP, 100.0);
        toplevel.setProperty(CoreOptions.ALGORITHM, "org.eclipse.elk.fixed");
        toplevel.setProperty(CoreOptions.PADDING, new ElkPadding());
        toplevel.setProperty(CoreOptions.SPACING_NODE_NODE, 0.0);
        // set size of node so that children will be scaled down
        toplevel.setProperty(CoreOptions.TOPDOWN_HIERARCHICAL_NODE_WIDTH, 300.0);
        toplevel.setProperty(CoreOptions.TOPDOWN_HIERARCHICAL_NODE_ASPECT_RATIO, 1.0);
        
        ElkNode child1 = ElkGraphUtil.createNode(toplevel);
        child1.setProperty(CoreOptions.TOPDOWN_LAYOUT, true);
        child1.setProperty(CoreOptions.TOPDOWN_NODE_TYPE, TopdownNodeTypes.HIERARCHICAL_NODE);
        child1.setProperty(CoreOptions.NODE_SIZE_FIXED_GRAPH_SIZE, true);
        child1.setX(0);
        child1.setY(0);
        child1.setDimensions(30, 30);
        
        ElkNode child2 = ElkGraphUtil.createNode(toplevel);
        child2.setProperty(CoreOptions.TOPDOWN_LAYOUT, true);
        child2.setProperty(CoreOptions.TOPDOWN_NODE_TYPE, TopdownNodeTypes.HIERARCHICAL_NODE);
        child2.setProperty(CoreOptions.NODE_SIZE_FIXED_GRAPH_SIZE, true);
        child2.setX(0);
        child2.setY(40);
        child2.setDimensions(30, 30);
        
        
        // prepare layout engine
        LayoutConfigurator config = new LayoutConfigurator();
        ElkUtil.applyVisitors(graph, config, new LayoutAlgorithmResolver());
        // call layout with layout engine
        try {
            new RecursiveGraphLayoutEngine().layout(graph, new BasicProgressMonitor());
        } catch (UnsupportedGraphException exception) {
            fail(exception.toString());
        }
               
        // topdown scale factor of toplevel node
        assertEquals(300.0 / 70.0, toplevel.getProperty(CoreOptions.TOPDOWN_SCALE_FACTOR), 0.00001);
    }


    /**
     * Test scale cap property when the scale cap bounds the resulting topdown scale factor.
     */
    @Test
    public void testScaleCapBounded() {
        PlainJavaInitialization.initializePlainJavaLayout();
        ElkNode graph = ElkGraphUtil.createGraph();
        graph.setProperty(CoreOptions.TOPDOWN_LAYOUT, true);
        graph.setProperty(CoreOptions.TOPDOWN_NODE_TYPE, TopdownNodeTypes.ROOT_NODE);
        
        ElkNode toplevel = ElkGraphUtil.createNode(graph);
        toplevel.setProperty(CoreOptions.TOPDOWN_LAYOUT, true);
        toplevel.setProperty(CoreOptions.TOPDOWN_NODE_TYPE, TopdownNodeTypes.HIERARCHICAL_NODE);
        toplevel.setProperty(CoreOptions.NODE_SIZE_FIXED_GRAPH_SIZE, true);
        toplevel.setProperty(CoreOptions.TOPDOWN_SCALE_CAP, 3.0);
        toplevel.setProperty(CoreOptions.ALGORITHM, "org.eclipse.elk.fixed");
        toplevel.setProperty(CoreOptions.PADDING, new ElkPadding());
        toplevel.setProperty(CoreOptions.SPACING_NODE_NODE, 0.0);
        // set size of node so that children will be scaled down
        toplevel.setProperty(CoreOptions.TOPDOWN_HIERARCHICAL_NODE_WIDTH, 300.0);
        toplevel.setProperty(CoreOptions.TOPDOWN_HIERARCHICAL_NODE_ASPECT_RATIO, 1.0);
        
        ElkNode child1 = ElkGraphUtil.createNode(toplevel);
        child1.setProperty(CoreOptions.TOPDOWN_LAYOUT, true);
        child1.setProperty(CoreOptions.TOPDOWN_NODE_TYPE, TopdownNodeTypes.HIERARCHICAL_NODE);
        child1.setProperty(CoreOptions.NODE_SIZE_FIXED_GRAPH_SIZE, true);
        child1.setX(0);
        child1.setY(0);
        child1.setDimensions(30, 30);
        
        ElkNode child2 = ElkGraphUtil.createNode(toplevel);
        child2.setProperty(CoreOptions.TOPDOWN_LAYOUT, true);
        child2.setProperty(CoreOptions.TOPDOWN_NODE_TYPE, TopdownNodeTypes.HIERARCHICAL_NODE);
        child2.setProperty(CoreOptions.NODE_SIZE_FIXED_GRAPH_SIZE, true);
        child2.setX(0);
        child2.setY(40);
        child2.setDimensions(30, 30);
        
        
        // prepare layout engine
        LayoutConfigurator config = new LayoutConfigurator();
        ElkUtil.applyVisitors(graph, config, new LayoutAlgorithmResolver());
        // call layout with layout engine
        try {
            new RecursiveGraphLayoutEngine().layout(graph, new BasicProgressMonitor());
        } catch (UnsupportedGraphException exception) {
            fail(exception.toString());
        }
               
        // topdown scale factor of toplevel node
        assertEquals(3.0, toplevel.getProperty(CoreOptions.TOPDOWN_SCALE_FACTOR), 0.00001);
    }
}
