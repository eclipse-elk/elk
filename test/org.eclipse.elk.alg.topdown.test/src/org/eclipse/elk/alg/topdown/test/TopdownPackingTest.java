/*******************************************************************************
 * Copyright (c) 2022 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 * 
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later 
 *******************************************************************************/
package org.eclipse.elk.alg.topdown.test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;

import org.eclipse.elk.alg.test.PlainJavaInitialization;
import org.eclipse.elk.alg.topdownpacking.TopdownpackingLayoutProvider;
import org.eclipse.elk.core.math.ElkPadding;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.util.BasicProgressMonitor;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.util.ElkGraphUtil;
import org.junit.Test;

/**
 * Tests the topdown packing algorithm. The node sizes must be externally set and determine the layout results.
 * Unless otherwise specified, the tests expect the default values of properties to be used.
 */
public class TopdownPackingTest {

    // CHECKSTYLEOFF MagicNumber

    /**
     * Test topdown packing on an empty graph.
     */
    @Test
    public void testEmptyGraph() {
        PlainJavaInitialization.initializePlainJavaLayout();
        ElkNode graph = createEmptyGraph();
        TopdownpackingLayoutProvider layoutProvider = new TopdownpackingLayoutProvider();
        layoutProvider.layout(graph, new BasicProgressMonitor());
    }
    
    /**
     * Test topdown packing for a graph with two nodes. In this case the resulting graph
     * size should be vertically compacted.
     */
    @Test
    public void testTwoNodes() {
        PlainJavaInitialization.initializePlainJavaLayout();
        ElkNode graph = createGraph(2);
        // property default values
        double hierarchicalWidth = graph.getProperty(CoreOptions.TOPDOWN_HIERARCHICAL_NODE_WIDTH);
        double hierarchicalAspectRatio = graph.getProperty(CoreOptions.TOPDOWN_HIERARCHICAL_NODE_ASPECT_RATIO);
        ElkPadding padding = graph.getProperty(CoreOptions.PADDING);
        double nodeNodeSpacing = graph.getProperty(CoreOptions.SPACING_NODE_NODE);
        TopdownpackingLayoutProvider layoutProvider = new TopdownpackingLayoutProvider();

        // test size prediction
        KVector predictedSize = layoutProvider.getPredictedGraphSize(graph);
        double expectedWidth = padding.left + 2 * hierarchicalWidth + nodeNodeSpacing + padding.right;
        assertEquals(expectedWidth, predictedSize.x, 0.00001);

        double expectedHeight = padding.top + hierarchicalWidth / hierarchicalAspectRatio + padding.bottom;
        assertEquals(expectedHeight, predictedSize.y, 0.00001);


        // test actual layout
        layoutProvider.layout(graph, new BasicProgressMonitor());
        double childOneX = graph.getChildren().get(0).getX();
        double childOneY = graph.getChildren().get(0).getY();
        double childOneWidth = graph.getChildren().get(0).getWidth();
        double childOneHeight = graph.getChildren().get(0).getHeight();
        
        assertEquals(padding.left, childOneX, 0.00001);
        assertEquals(padding.top, childOneY, 0.00001);
        assertEquals(hierarchicalWidth, childOneWidth, 0.00001);
        assertEquals(hierarchicalWidth / hierarchicalAspectRatio, childOneHeight, 0.00001);
        
        double childTwoX = graph.getChildren().get(1).getX();
        double childTwoY = graph.getChildren().get(1).getY();
        double childTwoWidth = graph.getChildren().get(1).getWidth();
        double childTwoHeight = graph.getChildren().get(1).getHeight();
        
        assertEquals(padding.left + hierarchicalWidth + nodeNodeSpacing, childTwoX, 0.00001);
        assertEquals(padding.top, childTwoY, 0.00001);
        assertEquals(hierarchicalWidth, childTwoWidth, 0.00001);
        assertEquals(hierarchicalWidth / hierarchicalAspectRatio, childTwoHeight, 0.00001);
    }
    
    /**
     * Test topdown packing for a graph with three nodes. In this case the third node should
     * be expanded to fill the whitespace of the parent.
     */
    @Test
    public void testThreeNodes() {
        PlainJavaInitialization.initializePlainJavaLayout();
        ElkNode graph = createGraph(3);
        // property default values
        double hierarchicalWidth = graph.getProperty(CoreOptions.TOPDOWN_HIERARCHICAL_NODE_WIDTH);
        double hierarchicalAspectRatio = graph.getProperty(CoreOptions.TOPDOWN_HIERARCHICAL_NODE_ASPECT_RATIO);
        ElkPadding padding = graph.getProperty(CoreOptions.PADDING);
        double nodeNodeSpacing = graph.getProperty(CoreOptions.SPACING_NODE_NODE);
        TopdownpackingLayoutProvider layoutProvider = new TopdownpackingLayoutProvider();
        
        // test size prediction
        KVector predictedSize = layoutProvider.getPredictedGraphSize(graph);
        double expectedWidth = padding.left + 2 * hierarchicalWidth + nodeNodeSpacing + padding.right;
        assertEquals(expectedWidth, predictedSize.x, 0.00001);
        
        double expectedHeight = padding.top + 2 * (hierarchicalWidth / hierarchicalAspectRatio) 
                + nodeNodeSpacing + padding.bottom;
        assertEquals(expectedHeight, predictedSize.y, 0.00001);
        
        // test actual layout
        layoutProvider.layout(graph, new BasicProgressMonitor());
        double childOneX = graph.getChildren().get(0).getX();
        double childOneY = graph.getChildren().get(0).getY();
        double childOneWidth = graph.getChildren().get(0).getWidth();
        double childOneHeight = graph.getChildren().get(0).getHeight();
        
        assertEquals(padding.left, childOneX, 0.00001);
        assertEquals(padding.top, childOneY, 0.00001);
        assertEquals(hierarchicalWidth, childOneWidth, 0.00001);
        assertEquals(hierarchicalWidth / hierarchicalAspectRatio, childOneHeight, 0.00001);
        
        double childTwoX = graph.getChildren().get(1).getX();
        double childTwoY = graph.getChildren().get(1).getY();
        double childTwoWidth = graph.getChildren().get(1).getWidth();
        double childTwoHeight = graph.getChildren().get(1).getHeight();
        
        assertEquals(padding.left + hierarchicalWidth + nodeNodeSpacing, childTwoX, 0.00001);
        assertEquals(padding.top, childTwoY, 0.00001);
        assertEquals(hierarchicalWidth, childTwoWidth, 0.00001);
        assertEquals(hierarchicalWidth / hierarchicalAspectRatio, childTwoHeight, 0.00001);
        
        // test third node is expanded correctly
        double childThreeX = graph.getChildren().get(2).getX();
        double childThreeY = graph.getChildren().get(2).getY();
        double childThreeWidth = graph.getChildren().get(2).getWidth();
        double childThreeHeight = graph.getChildren().get(2).getHeight();
        
        assertEquals(padding.left, childThreeX, 0.00001);
        assertEquals(padding.top + (hierarchicalWidth / hierarchicalAspectRatio) 
                + nodeNodeSpacing, childThreeY, 0.00001);
        assertEquals(2 * hierarchicalWidth + nodeNodeSpacing, childThreeWidth, 0.00001);
        assertEquals(hierarchicalWidth / hierarchicalAspectRatio, childThreeHeight, 0.00001);
    }
    
    /**
     * Test topdown packing for a graph with five nodes. In this case there is both vertical
     * compaction and whitespace elimination in the second row (involving two nodes).
     */
    @Test
    public void testFiveNodes() {
        PlainJavaInitialization.initializePlainJavaLayout();
        ElkNode graph = createGraph(5);
        // property default values
        double hierarchicalWidth = graph.getProperty(CoreOptions.TOPDOWN_HIERARCHICAL_NODE_WIDTH);
        double hierarchicalAspectRatio = graph.getProperty(CoreOptions.TOPDOWN_HIERARCHICAL_NODE_ASPECT_RATIO);
        ElkPadding padding = graph.getProperty(CoreOptions.PADDING);
        double nodeNodeSpacing = graph.getProperty(CoreOptions.SPACING_NODE_NODE);
        TopdownpackingLayoutProvider layoutProvider = new TopdownpackingLayoutProvider();
        
        // test size prediction
        KVector predictedSize = layoutProvider.getPredictedGraphSize(graph);
        double expectedWidth = padding.left + 3 * hierarchicalWidth + 2 * nodeNodeSpacing + padding.right;
        assertEquals(expectedWidth, predictedSize.x, 0.00001);
        
        double expectedHeight = padding.top + 2 * (hierarchicalWidth / hierarchicalAspectRatio) 
                + nodeNodeSpacing + padding.bottom;
        assertEquals(expectedHeight, predictedSize.y, 0.00001);
        
        // test white space expansion in bottom row is performed correctly
        layoutProvider.layout(graph, new BasicProgressMonitor());
        double childFourX = graph.getChildren().get(3).getX();
        double childFourY = graph.getChildren().get(3).getY();
        double childFourWidth = graph.getChildren().get(3).getWidth();
        double childFourHeight = graph.getChildren().get(3).getHeight();
        
        // additional half of empty space
        double expandedWidth = hierarchicalWidth + 0.5 * (hierarchicalWidth + nodeNodeSpacing);
        
        assertEquals(padding.left, childFourX, 0.00001);
        assertEquals(padding.top + (hierarchicalWidth / hierarchicalAspectRatio) 
                + nodeNodeSpacing, childFourY, 0.00001);
        assertEquals(expandedWidth, childFourWidth, 0.00001);
        assertEquals(hierarchicalWidth / hierarchicalAspectRatio, childFourHeight, 0.00001);
        
        layoutProvider.layout(graph, new BasicProgressMonitor());
        double childFiveX = graph.getChildren().get(4).getX();
        double childFiveY = graph.getChildren().get(4).getY();
        double childFiveWidth = graph.getChildren().get(4).getWidth();
        double childFiveHeight = graph.getChildren().get(4).getHeight();
        
        assertEquals(padding.left + expandedWidth + nodeNodeSpacing, childFiveX, 0.00001);
        assertEquals(padding.top + (hierarchicalWidth / hierarchicalAspectRatio) 
                + nodeNodeSpacing, childFiveY, 0.00001);
        assertEquals(expandedWidth, childFiveWidth, 0.00001);
        assertEquals(hierarchicalWidth / hierarchicalAspectRatio, childFiveHeight, 0.00001);
        
    }
    
    private static ElkNode createEmptyGraph() {
        return ElkGraphUtil.createGraph();
    }
    
    private static ElkNode createGraph(final int numberOfNodes) {
        if (numberOfNodes < 0) {
            fail();
        }
        ElkNode graph = ElkGraphUtil.createGraph();
        for (int i = 0; i < numberOfNodes; i++) {
            ElkGraphUtil.createNode(graph);
        }
        return graph;
    }
}
