/*******************************************************************************
 * Copyright (c) 2024 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 * 
 * SPDX-License-Identifier: EPL-2.0 
 *******************************************************************************/
package org.eclipse.elk.alg.layered.intermediate;

import static org.junit.Assert.assertTrue;
import org.eclipse.elk.alg.layered.LayeredLayoutProvider;
import org.eclipse.elk.alg.layered.options.LayerUnzippingStrategy;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.layered.options.OrderingStrategy;
import org.eclipse.elk.alg.test.PlainJavaInitialization;
import org.eclipse.elk.core.util.NullElkProgressMonitor;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.util.ElkGraphUtil;
import org.junit.Before;
import org.junit.Test;

/**
 * Tests the general n-way alternating layer unzipper.
 *
 */
public class AlternatingLayerUnzipperTest {
    
    LayeredLayoutProvider layeredLayout;
    
    @Before
    public void setUp() {
        PlainJavaInitialization.initializePlainJavaLayout();
        layeredLayout = new LayeredLayoutProvider();
    }

    /**
     * Tests splitting a layer of three nodes into two layers in an A-B-A pattern.
     */
    @Test
    public void simpleTwoSplit() {
        ElkNode graph = ElkGraphUtil.createGraph();

        ElkNode node1 = ElkGraphUtil.createNode(graph);
        node1.setWidth(30);
        node1.setHeight(30);
        

        ElkNode node2 = ElkGraphUtil.createNode(graph);
        node2.setWidth(30);
        node2.setHeight(30);
        
        ElkNode node3 = ElkGraphUtil.createNode(graph);
        node3.setWidth(30);
        node3.setHeight(30);
        
        ElkNode node4 = ElkGraphUtil.createNode(graph);
        node4.setWidth(30);
        node4.setHeight(30);
        
        ElkNode nodeFinal = ElkGraphUtil.createNode(graph);
        nodeFinal.setWidth(30);
        nodeFinal.setHeight(30);

        ElkGraphUtil.createSimpleEdge(node1, node2);
        ElkGraphUtil.createSimpleEdge(node1, node3);
        ElkGraphUtil.createSimpleEdge(node1, node4);
        ElkGraphUtil.createSimpleEdge(node2, nodeFinal);
        ElkGraphUtil.createSimpleEdge(node3, nodeFinal);
        ElkGraphUtil.createSimpleEdge(node4, nodeFinal);

        
        graph.setProperty(LayeredOptions.LAYER_UNZIPPING_STRATEGY, LayerUnzippingStrategy.ALTERNATING);
        graph.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY, OrderingStrategy.PREFER_EDGES);
        node2.setProperty(LayeredOptions.LAYER_UNZIPPING_LAYER_SPLIT, 2);
        
        layeredLayout.layout(graph, new NullElkProgressMonitor());
        
        assertTrue(node1.getX() < node2.getX());
        assertTrue(node1.getX() < node4.getX());
        assertTrue(node2.getX() < node3.getX());
        assertTrue(node4.getX() < node3.getX());
        assertTrue(node3.getX() < nodeFinal.getX());
    }
    
    /**
     * Tests splitting a layer of 4 nodes into 3 layers in an A-B-C-A pattern.
     */
    @Test
    public void simpleThreeSplit() {
        ElkNode graph = ElkGraphUtil.createGraph();

        ElkNode node1 = ElkGraphUtil.createNode(graph);
        node1.setWidth(30);
        node1.setHeight(30);
        

        ElkNode node2 = ElkGraphUtil.createNode(graph);
        node2.setWidth(30);
        node2.setHeight(30);
        
        ElkNode node3 = ElkGraphUtil.createNode(graph);
        node3.setWidth(30);
        node3.setHeight(30);
        
        ElkNode node4 = ElkGraphUtil.createNode(graph);
        node4.setWidth(30);
        node4.setHeight(30);
        
        ElkNode node5 = ElkGraphUtil.createNode(graph);
        node5.setWidth(30);
        node5.setHeight(30);
        
        ElkNode nodeFinal = ElkGraphUtil.createNode(graph);
        nodeFinal.setWidth(30);
        nodeFinal.setHeight(30);

        ElkGraphUtil.createSimpleEdge(node1, node2);
        ElkGraphUtil.createSimpleEdge(node1, node3);
        ElkGraphUtil.createSimpleEdge(node1, node4);
        ElkGraphUtil.createSimpleEdge(node1, node5);
        ElkGraphUtil.createSimpleEdge(node2, nodeFinal);
        ElkGraphUtil.createSimpleEdge(node3, nodeFinal);
        ElkGraphUtil.createSimpleEdge(node4, nodeFinal);
        ElkGraphUtil.createSimpleEdge(node5, nodeFinal);
        
        graph.setProperty(LayeredOptions.LAYER_UNZIPPING_STRATEGY, LayerUnzippingStrategy.ALTERNATING);
        graph.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY, OrderingStrategy.PREFER_EDGES);
        node2.setProperty(LayeredOptions.LAYER_UNZIPPING_LAYER_SPLIT, 3);
        
        layeredLayout.layout(graph, new NullElkProgressMonitor());
        
        assertTrue(node1.getX() < node2.getX());
        assertTrue(node1.getX() < node5.getX());
        assertTrue(node2.getX() < node3.getX());
        assertTrue(node3.getX() < node4.getX());
        assertTrue(node5.getX() < node4.getX());
        assertTrue(node4.getX() < nodeFinal.getX());
    }
    
    /**
     * Tests the case that the nodes of the split layer have no further outgoing edges i.e. there is no further 
     * connected layer after them.
     */
    @Test
    public void danglingOutgoing() {
        ElkNode graph = ElkGraphUtil.createGraph();

        ElkNode node1 = ElkGraphUtil.createNode(graph);
        node1.setWidth(30);
        node1.setHeight(30);
        

        ElkNode node2 = ElkGraphUtil.createNode(graph);
        node2.setWidth(30);
        node2.setHeight(30);
        
        ElkNode node3 = ElkGraphUtil.createNode(graph);
        node3.setWidth(30);
        node3.setHeight(30);
        
        ElkNode node4 = ElkGraphUtil.createNode(graph);
        node4.setWidth(30);
        node4.setHeight(30);
        
        ElkNode node5 = ElkGraphUtil.createNode(graph);
        node5.setWidth(30);
        node5.setHeight(30);

        ElkGraphUtil.createSimpleEdge(node1, node2);
        ElkGraphUtil.createSimpleEdge(node1, node3);
        ElkGraphUtil.createSimpleEdge(node1, node4);
        ElkGraphUtil.createSimpleEdge(node1, node5);
        
        graph.setProperty(LayeredOptions.LAYER_UNZIPPING_STRATEGY, LayerUnzippingStrategy.ALTERNATING);
        graph.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY, OrderingStrategy.PREFER_EDGES);
        node2.setProperty(LayeredOptions.LAYER_UNZIPPING_LAYER_SPLIT, 2);
        
        layeredLayout.layout(graph, new NullElkProgressMonitor());
        
        assertTrue(node1.getX() < node2.getX());
        assertTrue(node1.getX() < node4.getX());
        assertTrue(node2.getX() < node3.getX());
        assertTrue(node2.getX() < node5.getX());
        assertTrue(node4.getX() < node3.getX());
        assertTrue(node4.getX() < node5.getX());
    }
    
    /**
     * Tests the case that there is no connected layer leading into the layer that is split.
     */
    @Test
    public void danglingIncoming() {
        ElkNode graph = ElkGraphUtil.createGraph();

        ElkNode node1 = ElkGraphUtil.createNode(graph);
        node1.setWidth(30);
        node1.setHeight(30);
        

        ElkNode node2 = ElkGraphUtil.createNode(graph);
        node2.setWidth(30);
        node2.setHeight(30);
        
        ElkNode node3 = ElkGraphUtil.createNode(graph);
        node3.setWidth(30);
        node3.setHeight(30);
        
        ElkNode node4 = ElkGraphUtil.createNode(graph);
        node4.setWidth(30);
        node4.setHeight(30);
        
        ElkNode node5 = ElkGraphUtil.createNode(graph);
        node5.setWidth(30);
        node5.setHeight(30);
        
        ElkNode nodeFinal = ElkGraphUtil.createNode(graph);
        nodeFinal.setWidth(30);
        nodeFinal.setHeight(30);
        
        ElkGraphUtil.createSimpleEdge(node1, nodeFinal);
        ElkGraphUtil.createSimpleEdge(node2, nodeFinal);
        ElkGraphUtil.createSimpleEdge(node3, nodeFinal);
        ElkGraphUtil.createSimpleEdge(node4, nodeFinal);
        ElkGraphUtil.createSimpleEdge(node5, nodeFinal);
        
        graph.setProperty(LayeredOptions.LAYER_UNZIPPING_STRATEGY, LayerUnzippingStrategy.ALTERNATING);
        graph.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY, OrderingStrategy.PREFER_EDGES);
        node2.setProperty(LayeredOptions.LAYER_UNZIPPING_LAYER_SPLIT, 2);
        
        layeredLayout.layout(graph, new NullElkProgressMonitor());
        
        assertTrue(node1.getX() < node2.getX());
        assertTrue(node1.getX() < node4.getX());
        assertTrue(node3.getX() < node2.getX());
        assertTrue(node3.getX() < node4.getX());
        assertTrue(node5.getX() < node2.getX());
        assertTrue(node5.getX() < node4.getX());
        
        assertTrue(node2.getX() < nodeFinal.getX());
        assertTrue(node4.getX() < nodeFinal.getX());
    }

    /**
     * Tests a two layer graph being split into a four layer graph.
     */
    @Test
    public void multipleLayersSplit() {
        ElkNode graph = ElkGraphUtil.createGraph();

        // LAYER 1
        ElkNode node1 = ElkGraphUtil.createNode(graph);
        node1.setWidth(30);
        node1.setHeight(30);
        

        ElkNode node2 = ElkGraphUtil.createNode(graph);
        node2.setWidth(30);
        node2.setHeight(30);
        
        ElkNode node3 = ElkGraphUtil.createNode(graph);
        node3.setWidth(30);
        node3.setHeight(30);
        
        ElkNode node4 = ElkGraphUtil.createNode(graph);
        node4.setWidth(30);
        node4.setHeight(30);
        
        // LAYER 2
        ElkNode node21 = ElkGraphUtil.createNode(graph);
        node21.setWidth(30);
        node21.setHeight(30);
        
        ElkNode node22 = ElkGraphUtil.createNode(graph);
        node22.setWidth(30);
        node22.setHeight(30);
        
        ElkNode node23 = ElkGraphUtil.createNode(graph);
        node23.setWidth(30);
        node23.setHeight(30);
        
        ElkGraphUtil.createSimpleEdge(node1, node21);
        
        ElkGraphUtil.createSimpleEdge(node2, node21);
        ElkGraphUtil.createSimpleEdge(node2, node22);

        ElkGraphUtil.createSimpleEdge(node3, node21);
        ElkGraphUtil.createSimpleEdge(node3, node22);
        ElkGraphUtil.createSimpleEdge(node4, node22);
        ElkGraphUtil.createSimpleEdge(node4, node23);
        
        graph.setProperty(LayeredOptions.LAYER_UNZIPPING_STRATEGY, LayerUnzippingStrategy.ALTERNATING);
        graph.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY, OrderingStrategy.PREFER_EDGES);
        node1.setProperty(LayeredOptions.LAYER_UNZIPPING_LAYER_SPLIT, 2);
        node21.setProperty(LayeredOptions.LAYER_UNZIPPING_LAYER_SPLIT, 2);
        
        layeredLayout.layout(graph, new NullElkProgressMonitor());
        
        assertTrue(node1.getX() < node2.getX());
        assertTrue(node1.getX() < node4.getX());
        assertTrue(node3.getX() < node2.getX());
        assertTrue(node3.getX() < node4.getX());

        assertTrue(node2.getX() < node21.getX());
        assertTrue(node2.getX() < node23.getX());
        assertTrue(node4.getX() < node21.getX());
        assertTrue(node4.getX() < node23.getX());
        
        assertTrue(node21.getX() < node22.getX());
        assertTrue(node23.getX() < node22.getX());
    }
    
    /**
     * Tests the case that there are multiple incoming edges connecting to a node that is in a layer being split.
     */
    @Test
    public void multipleIncomingEdges() {
        ElkNode graph = ElkGraphUtil.createGraph();

        ElkNode node1 = ElkGraphUtil.createNode(graph);
        node1.setWidth(30);
        node1.setHeight(30);
        

        ElkNode node2 = ElkGraphUtil.createNode(graph);
        node2.setWidth(30);
        node2.setHeight(30);
        
        ElkNode node3 = ElkGraphUtil.createNode(graph);
        node3.setWidth(30);
        node3.setHeight(30);
        
        ElkNode node4 = ElkGraphUtil.createNode(graph);
        node4.setWidth(30);
        node4.setHeight(30);
        
        ElkNode nodeFinal = ElkGraphUtil.createNode(graph);
        nodeFinal.setWidth(30);
        nodeFinal.setHeight(30);

        ElkGraphUtil.createSimpleEdge(node1, node2);
        ElkGraphUtil.createSimpleEdge(node1, node2);
        ElkGraphUtil.createSimpleEdge(node1, node2);
        ElkGraphUtil.createSimpleEdge(node1, node2);
        ElkGraphUtil.createSimpleEdge(node1, node2);
        ElkGraphUtil.createSimpleEdge(node1, node3);
        ElkGraphUtil.createSimpleEdge(node1, node3);
        ElkGraphUtil.createSimpleEdge(node1, node3);
        ElkGraphUtil.createSimpleEdge(node1, node3);
        ElkGraphUtil.createSimpleEdge(node1, node4);
        ElkGraphUtil.createSimpleEdge(node2, nodeFinal);
        ElkGraphUtil.createSimpleEdge(node3, nodeFinal);
        ElkGraphUtil.createSimpleEdge(node4, nodeFinal);

        
        graph.setProperty(LayeredOptions.LAYER_UNZIPPING_STRATEGY, LayerUnzippingStrategy.ALTERNATING);
        graph.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY, OrderingStrategy.PREFER_EDGES);
        node2.setProperty(LayeredOptions.LAYER_UNZIPPING_LAYER_SPLIT, 2);
        
        layeredLayout.layout(graph, new NullElkProgressMonitor());
        
        assertTrue(node1.getX() < node2.getX());
        assertTrue(node1.getX() < node4.getX());
        assertTrue(node2.getX() < node3.getX());
        assertTrue(node4.getX() < node3.getX());
        assertTrue(node3.getX() < nodeFinal.getX());
    }
    
    /**
     * Tests the case that there are multiple outgoing edges connecting to a node that is in a layer being split.
     */
    @Test
    public void multipleOutgoingEdges() {
        ElkNode graph = ElkGraphUtil.createGraph();

        ElkNode node1 = ElkGraphUtil.createNode(graph);
        node1.setWidth(30);
        node1.setHeight(30);
        

        ElkNode node2 = ElkGraphUtil.createNode(graph);
        node2.setWidth(30);
        node2.setHeight(30);
        
        ElkNode node3 = ElkGraphUtil.createNode(graph);
        node3.setWidth(30);
        node3.setHeight(30);
        
        ElkNode node4 = ElkGraphUtil.createNode(graph);
        node4.setWidth(30);
        node4.setHeight(30);
        
        ElkNode nodeFinal = ElkGraphUtil.createNode(graph);
        nodeFinal.setWidth(30);
        nodeFinal.setHeight(30);

        ElkGraphUtil.createSimpleEdge(node1, node2);
        ElkGraphUtil.createSimpleEdge(node1, node3);
        ElkGraphUtil.createSimpleEdge(node1, node4);
        ElkGraphUtil.createSimpleEdge(node2, nodeFinal);
        ElkGraphUtil.createSimpleEdge(node2, nodeFinal);
        ElkGraphUtil.createSimpleEdge(node2, nodeFinal);
        ElkGraphUtil.createSimpleEdge(node2, nodeFinal);
        ElkGraphUtil.createSimpleEdge(node3, nodeFinal);
        ElkGraphUtil.createSimpleEdge(node4, nodeFinal);
        ElkGraphUtil.createSimpleEdge(node4, nodeFinal);
        ElkGraphUtil.createSimpleEdge(node4, nodeFinal);
        ElkGraphUtil.createSimpleEdge(node4, nodeFinal);

        
        graph.setProperty(LayeredOptions.LAYER_UNZIPPING_STRATEGY, LayerUnzippingStrategy.ALTERNATING);
        graph.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY, OrderingStrategy.PREFER_EDGES);
        node2.setProperty(LayeredOptions.LAYER_UNZIPPING_LAYER_SPLIT, 2);
        
        layeredLayout.layout(graph, new NullElkProgressMonitor());
        
        assertTrue(node1.getX() < node2.getX());
        assertTrue(node1.getX() < node4.getX());
        assertTrue(node2.getX() < node3.getX());
        assertTrue(node4.getX() < node3.getX());
        assertTrue(node3.getX() < nodeFinal.getX());
    }
    
    /**
     * Tests the case that some nodes are not connected to a preceding layer and some are.
     */
    @Test
    public void mixedDanglingIncoming() {
        ElkNode graph = ElkGraphUtil.createGraph();

        ElkNode node1 = ElkGraphUtil.createNode(graph);
        node1.setWidth(30);
        node1.setHeight(30);
        

        ElkNode node2 = ElkGraphUtil.createNode(graph);
        node2.setWidth(30);
        node2.setHeight(30);
        
        ElkNode node3 = ElkGraphUtil.createNode(graph);
        node3.setWidth(30);
        node3.setHeight(30);
        
        ElkNode node4 = ElkGraphUtil.createNode(graph);
        node4.setWidth(30);
        node4.setHeight(30);
        
        ElkNode node5 = ElkGraphUtil.createNode(graph);
        node5.setWidth(30);
        node5.setHeight(30);
        
        ElkNode node6 = ElkGraphUtil.createNode(graph);
        node6.setWidth(30);
        node6.setHeight(30);
        
        ElkNode node7 = ElkGraphUtil.createNode(graph);
        node7.setWidth(30);
        node7.setHeight(30);
        
        ElkGraphUtil.createSimpleEdge(node1, node3);
        ElkGraphUtil.createSimpleEdge(node1, node4);
        
        ElkGraphUtil.createSimpleEdge(node2, node5);
        ElkGraphUtil.createSimpleEdge(node3, node5);
        ElkGraphUtil.createSimpleEdge(node3, node6);
        ElkGraphUtil.createSimpleEdge(node4, node7);
        
        graph.setProperty(LayeredOptions.LAYER_UNZIPPING_STRATEGY, LayerUnzippingStrategy.ALTERNATING);
        graph.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY, OrderingStrategy.PREFER_EDGES);
        node2.setProperty(LayeredOptions.LAYER_UNZIPPING_LAYER_SPLIT, 2);
        node5.setProperty(LayeredOptions.LAYER_UNZIPPING_LAYER_SPLIT, 2);
        
        layeredLayout.layout(graph, new NullElkProgressMonitor());
        
        assertTrue(node1.getX() < node2.getX());
        assertTrue(node1.getX() < node4.getX());
        
        assertTrue(node2.getX() < node3.getX());
        assertTrue(node4.getX() < node3.getX());
        
        assertTrue(node3.getX() < node5.getX());
        assertTrue(node3.getX() < node7.getX());
        
        assertTrue(node5.getX() < node6.getX());
        assertTrue(node7.getX() < node6.getX());
    }
    
    /**
     * Tests the case that some nodes are not connected to a following layer and some are.
     */
    @Test
    public void mixedDanglingOutgoing() {
        ElkNode graph = ElkGraphUtil.createGraph();

        ElkNode node1 = ElkGraphUtil.createNode(graph);
        node1.setWidth(30);
        node1.setHeight(30);
        

        ElkNode node2 = ElkGraphUtil.createNode(graph);
        node2.setWidth(30);
        node2.setHeight(30);
        
        ElkNode node3 = ElkGraphUtil.createNode(graph);
        node3.setWidth(30);
        node3.setHeight(30);
        
        ElkNode node4 = ElkGraphUtil.createNode(graph);
        node4.setWidth(30);
        node4.setHeight(30);
        
        ElkNode node5 = ElkGraphUtil.createNode(graph);
        node5.setWidth(30);
        node5.setHeight(30);
        
        ElkNode node6 = ElkGraphUtil.createNode(graph);
        node6.setWidth(30);
        node6.setHeight(30);
        
        ElkNode node7 = ElkGraphUtil.createNode(graph);
        node7.setWidth(30);
        node7.setHeight(30);
        
        ElkGraphUtil.createSimpleEdge(node1, node2);
        ElkGraphUtil.createSimpleEdge(node1, node3);
        ElkGraphUtil.createSimpleEdge(node1, node4);

        ElkGraphUtil.createSimpleEdge(node3, node5);
        ElkGraphUtil.createSimpleEdge(node3, node6);
        ElkGraphUtil.createSimpleEdge(node4, node7);
        
        graph.setProperty(LayeredOptions.LAYER_UNZIPPING_STRATEGY, LayerUnzippingStrategy.ALTERNATING);
        graph.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY, OrderingStrategy.PREFER_EDGES);
        node2.setProperty(LayeredOptions.LAYER_UNZIPPING_LAYER_SPLIT, 2);
        node5.setProperty(LayeredOptions.LAYER_UNZIPPING_LAYER_SPLIT, 2);
        
        layeredLayout.layout(graph, new NullElkProgressMonitor());
        
        assertTrue(node1.getX() < node2.getX());
        assertTrue(node1.getX() < node4.getX());
        
        assertTrue(node2.getX() < node3.getX());
        assertTrue(node4.getX() < node3.getX());
        
        assertTrue(node3.getX() < node5.getX());
        assertTrue(node3.getX() < node7.getX());
        
        assertTrue(node5.getX() < node6.getX());
        assertTrue(node7.getX() < node6.getX());
    }
    
    /**
     * Tests the case that one layer does a 2-split and another layer does a 3-split.
     */
    @Test
    public void mixedTwoThreeLayerSplit() {
        ElkNode graph = ElkGraphUtil.createGraph();

        // LAYER 1
        ElkNode node1 = ElkGraphUtil.createNode(graph);
        node1.setWidth(30);
        node1.setHeight(30);
        

        ElkNode node2 = ElkGraphUtil.createNode(graph);
        node2.setWidth(30);
        node2.setHeight(30);
        
        ElkNode node3 = ElkGraphUtil.createNode(graph);
        node3.setWidth(30);
        node3.setHeight(30);
        
        ElkNode node4 = ElkGraphUtil.createNode(graph);
        node4.setWidth(30);
        node4.setHeight(30);
        
        // LAYER 2
        ElkNode node21 = ElkGraphUtil.createNode(graph);
        node21.setWidth(30);
        node21.setHeight(30);
        
        ElkNode node22 = ElkGraphUtil.createNode(graph);
        node22.setWidth(30);
        node22.setHeight(30);
        
        ElkNode node23 = ElkGraphUtil.createNode(graph);
        node23.setWidth(30);
        node23.setHeight(30);
        
        ElkGraphUtil.createSimpleEdge(node1, node21);
        
        ElkGraphUtil.createSimpleEdge(node2, node21);
        ElkGraphUtil.createSimpleEdge(node2, node22);

        ElkGraphUtil.createSimpleEdge(node3, node21);
        ElkGraphUtil.createSimpleEdge(node3, node22);
        ElkGraphUtil.createSimpleEdge(node4, node22);
        ElkGraphUtil.createSimpleEdge(node4, node23);
        
        graph.setProperty(LayeredOptions.LAYER_UNZIPPING_STRATEGY, LayerUnzippingStrategy.ALTERNATING);
        graph.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY, OrderingStrategy.PREFER_EDGES);
        node1.setProperty(LayeredOptions.LAYER_UNZIPPING_LAYER_SPLIT, 3);
        node21.setProperty(LayeredOptions.LAYER_UNZIPPING_LAYER_SPLIT, 2);
        
        layeredLayout.layout(graph, new NullElkProgressMonitor());
        
        assertTrue(node1.getX() < node2.getX());
        assertTrue(node2.getX() < node3.getX());
        assertTrue(node4.getX() < node2.getX());

        assertTrue(node3.getX() < node21.getX());
        assertTrue(node3.getX() < node23.getX());
        
        assertTrue(node21.getX() < node22.getX());
        assertTrue(node23.getX() < node22.getX());
    }
    
    /**
     * Tests that the staggered placement of the nodes resets when long edge crosses the layer.
     */
    @Test
    public void resetOnLongEdges() {
        ElkNode graph = ElkGraphUtil.createGraph();

        ElkNode node1 = ElkGraphUtil.createNode(graph);
        node1.setWidth(30);
        node1.setHeight(30);
        

        ElkNode node2 = ElkGraphUtil.createNode(graph);
        node2.setWidth(30);
        node2.setHeight(30);
        
        ElkNode node3 = ElkGraphUtil.createNode(graph);
        node3.setWidth(30);
        node3.setHeight(30);
        
        ElkNode node4 = ElkGraphUtil.createNode(graph);
        node4.setWidth(30);
        node4.setHeight(30);
        
        ElkNode node5 = ElkGraphUtil.createNode(graph);
        node5.setWidth(30);
        node5.setHeight(30);
        
        ElkGraphUtil.createSimpleEdge(node1, node2);
        ElkGraphUtil.createSimpleEdge(node1, node3);
        // long edge
        ElkGraphUtil.createSimpleEdge(node1, node5);
        ElkGraphUtil.createSimpleEdge(node1, node4);
        
        ElkGraphUtil.createSimpleEdge(node4, node5);
        
        graph.setProperty(LayeredOptions.LAYER_UNZIPPING_STRATEGY, LayerUnzippingStrategy.ALTERNATING);
        graph.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY, OrderingStrategy.PREFER_EDGES);
        node2.setProperty(LayeredOptions.LAYER_UNZIPPING_LAYER_SPLIT, 2);
        node2.setProperty(LayeredOptions.LAYER_UNZIPPING_RESET_ON_LONG_EDGES, true);
        
        layeredLayout.layout(graph, new NullElkProgressMonitor());
        
        assertTrue(node1.getX() < node2.getX());
        assertTrue(node1.getX() < node4.getX());
        assertTrue(node2.getX() < node3.getX());
        assertTrue(node3.getX() < node5.getX());
        assertTrue(node4.getX() < node3.getX());
    }
    
    /**
     * Tests the case that the resetOnLongEdges option is disabled. In this case the staggering should not be reset.
     */
    @Test
    public void noResetOnLongEdges() {
        ElkNode graph = ElkGraphUtil.createGraph();

        ElkNode node1 = ElkGraphUtil.createNode(graph);
        node1.setWidth(30);
        node1.setHeight(30);
        

        ElkNode node2 = ElkGraphUtil.createNode(graph);
        node2.setWidth(30);
        node2.setHeight(30);
        
        ElkNode node3 = ElkGraphUtil.createNode(graph);
        node3.setWidth(30);
        node3.setHeight(30);
        
        ElkNode node4 = ElkGraphUtil.createNode(graph);
        node4.setWidth(30);
        node4.setHeight(30);
        
        ElkNode node5 = ElkGraphUtil.createNode(graph);
        node5.setWidth(30);
        node5.setHeight(30);
        
        ElkNode node6 = ElkGraphUtil.createNode(graph);
        node6.setWidth(30);
        node6.setHeight(30);
        
        ElkGraphUtil.createSimpleEdge(node1, node2);
        ElkGraphUtil.createSimpleEdge(node1, node3);
        // long edge
        ElkGraphUtil.createSimpleEdge(node1, node5);
        ElkGraphUtil.createSimpleEdge(node1, node4);
        
        ElkGraphUtil.createSimpleEdge(node4, node5);
        ElkGraphUtil.createSimpleEdge(node1, node6);
        
        graph.setProperty(LayeredOptions.LAYER_UNZIPPING_STRATEGY, LayerUnzippingStrategy.ALTERNATING);
        graph.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY, OrderingStrategy.PREFER_EDGES);
        node2.setProperty(LayeredOptions.LAYER_UNZIPPING_LAYER_SPLIT, 2);
        node2.setProperty(LayeredOptions.LAYER_UNZIPPING_RESET_ON_LONG_EDGES, false);
        
        layeredLayout.layout(graph, new NullElkProgressMonitor());
        
        assertTrue(node1.getX() < node2.getX());
        assertTrue(node1.getX() < node6.getX());
        assertTrue(node2.getX() < node3.getX());
        assertTrue(node2.getX() < node4.getX());
        assertTrue(node3.getX() < node5.getX());
        assertTrue(node4.getX() < node5.getX());
        assertTrue(node6.getX() < node3.getX());
        assertTrue(node6.getX() < node4.getX());
    }
}
