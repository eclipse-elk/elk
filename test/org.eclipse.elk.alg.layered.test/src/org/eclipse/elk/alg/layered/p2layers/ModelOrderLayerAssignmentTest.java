/*******************************************************************************
 * Copyright (c) 2020 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.layered.p2layers;

import static org.junit.Assert.assertEquals;

import org.eclipse.elk.alg.layered.LayeredLayoutProvider;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.layered.options.LayeringStrategy;
import org.eclipse.elk.alg.layered.options.NodePromotionStrategy;
import org.eclipse.elk.alg.layered.options.OrderingStrategy;
import org.eclipse.elk.alg.test.PlainJavaInitialization;
import org.eclipse.elk.core.math.ElkPadding;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.options.Direction;
import org.eclipse.elk.core.util.BasicProgressMonitor;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.util.ElkGraphUtil;
import org.junit.BeforeClass;
import org.junit.Test;

/**
 * Test whether the model order layer assignment together with model order node promotion works
 * for some selected examples.
 *
 */
public class ModelOrderLayerAssignmentTest {

    @BeforeClass
    public static void init() {
        PlainJavaInitialization.initializePlainJavaLayout();
    }

    @Test
    public void testForwardPromotionDisregardsEdgeLength() {
        // Desired layout
        // n1 -> n2 -> n3 -> n5
        //    =======> n4 ->        

        ElkNode parent = ElkGraphUtil.createGraph();
        ElkNode n1 = ElkGraphUtil.createNode(parent);
        ElkNode n2 = ElkGraphUtil.createNode(parent);
        ElkNode n3 = ElkGraphUtil.createNode(parent);
        ElkNode n4 = ElkGraphUtil.createNode(parent);
        ElkNode n5 = ElkGraphUtil.createNode(parent);
        n1.setDimensions(30, 30);
        n2.setDimensions(30, 30);
        n3.setDimensions(30, 30);
        n4.setDimensions(30, 30);
        n5.setDimensions(30, 30);

        parent.setProperty(CoreOptions.ALGORITHM, LayeredOptions.ALGORITHM_ID);
        parent.setProperty(CoreOptions.DIRECTION, Direction.RIGHT);
        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY, OrderingStrategy.PREFER_EDGES);
        parent.setProperty(LayeredOptions.LAYERING_STRATEGY, LayeringStrategy.LONGEST_PATH_SOURCE);
        parent.setProperty(LayeredOptions.LAYERING_NODE_PROMOTION_STRATEGY,
                NodePromotionStrategy.MODEL_ORDER_LEFT_TO_RIGHT);
        parent.setProperty(CoreOptions.PADDING, new ElkPadding(0.0));
        parent.setProperty(CoreOptions.SPACING_NODE_NODE, 10.0);
        parent.setProperty(LayeredOptions.SPACING_NODE_NODE_BETWEEN_LAYERS, 20.0);

        ElkGraphUtil.createSimpleEdge(n1, n2);
        ElkGraphUtil.createSimpleEdge(n1, n4);
        ElkGraphUtil.createSimpleEdge(n1, n4);
        ElkGraphUtil.createSimpleEdge(n2, n3);
        ElkGraphUtil.createSimpleEdge(n3, n5);
        ElkGraphUtil.createSimpleEdge(n4, n5);
        ElkGraphUtil.createLabel("n1", n1);
        ElkGraphUtil.createLabel("n2", n2);
        ElkGraphUtil.createLabel("n3", n3);
        ElkGraphUtil.createLabel("n4", n4);
        ElkGraphUtil.createLabel("n5", n5);
        
        LayeredLayoutProvider layoutProvider = new LayeredLayoutProvider();
        layoutProvider.layout(parent, new BasicProgressMonitor());

        // Assert the positions or the nodes.
        assertEquals("n1 on wrong x coordinate", 0, n1.getX(), 0.1);
        assertEquals("n1 on wrong y coordinate", 7.5, n1.getY(), 0.1);
        // Check whether position is node size of n1 + 3*edgeEdge spacing.
        assertEquals("n2 on wrong x coordinate", 30 + 30, n2.getX(), 0.1);
        assertEquals("n2 on wrong y coordinate", 0, n2.getY(), 0.1);
        // Check whether position is node size of n1 + 3*edgeEdge spacing + node size of n2 + between layer spacing.
        assertEquals("n3 on wrong x coordinate", 30 + 30 + 30 + 20, n3.getX(), 0.1);
        assertEquals("n3 on wrong y coordinate", 0, n3.getY(), 0.1);
        // Check whether n4 aligns with n3.
        assertEquals("n4 on wrong x coordinate", 110, n4.getX(), 0.1);
        assertEquals("n4 on wrong y coordinate", 30 + 10, n4.getY(), 0.1);
        assertEquals("n5 on wrong x coordinate", 30 + 30 + 30 + 20 + 30 + 20, n5.getX(), 0.1);
        assertEquals("n5 on wrong y coordinate", 5, n5.getY(), 0.1);
    }
    


    @Test
    public void testBackwardPromotionDisregardsEdgeLength() {
        // Desired layout
        // n1 -> n2 =======> n5
        //    -> n3 -> n4 ->        

        ElkNode parent = ElkGraphUtil.createGraph();
        ElkNode n1 = ElkGraphUtil.createNode(parent);
        ElkNode n2 = ElkGraphUtil.createNode(parent);
        ElkNode n3 = ElkGraphUtil.createNode(parent);
        ElkNode n4 = ElkGraphUtil.createNode(parent);
        ElkNode n5 = ElkGraphUtil.createNode(parent);
        n1.setDimensions(30, 30);
        n2.setDimensions(30, 30);
        n3.setDimensions(30, 30);
        n4.setDimensions(30, 30);
        n5.setDimensions(30, 30);
        ElkGraphUtil.createLabel("n1", n1);
        ElkGraphUtil.createLabel("n2", n2);
        ElkGraphUtil.createLabel("n3", n3);
        ElkGraphUtil.createLabel("n4", n4);
        ElkGraphUtil.createLabel("n5", n5);

        parent.setProperty(CoreOptions.ALGORITHM, LayeredOptions.ALGORITHM_ID);
        parent.setProperty(CoreOptions.DIRECTION, Direction.RIGHT);
        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY, OrderingStrategy.PREFER_EDGES);
        parent.setProperty(LayeredOptions.LAYERING_STRATEGY, LayeringStrategy.LONGEST_PATH);
        parent.setProperty(LayeredOptions.LAYERING_NODE_PROMOTION_STRATEGY,
                NodePromotionStrategy.MODEL_ORDER_RIGHT_TO_LEFT);
        parent.setProperty(CoreOptions.PADDING, new ElkPadding(0.0));
        parent.setProperty(CoreOptions.SPACING_NODE_NODE, 10.0);
        parent.setProperty(LayeredOptions.SPACING_NODE_NODE_BETWEEN_LAYERS, 20.0);

        ElkGraphUtil.createSimpleEdge(n1, n2);
        ElkGraphUtil.createSimpleEdge(n1, n3);
        ElkGraphUtil.createSimpleEdge(n2, n5);
        ElkGraphUtil.createSimpleEdge(n2, n5);
        ElkGraphUtil.createSimpleEdge(n3, n4);
        ElkGraphUtil.createSimpleEdge(n4, n5);
        
        LayeredLayoutProvider layoutProvider = new LayeredLayoutProvider();
        layoutProvider.layout(parent, new BasicProgressMonitor());

        // Assert the positions or the nodes.
        assertEquals("n1 on wrong x coordinate", 0, n1.getX(), 0.1);
        assertEquals("n1 on wrong y coordinate", 5, n1.getY(), 0.1);
        assertEquals("n2 on wrong x coordinate", 30 + 20, n2.getX(), 0.1);
        assertEquals("n2 on wrong y coordinate", 0, n2.getY(), 0.1);
        // n3 aligns wit n2.
        assertEquals("n3 on wrong x coordinate", 50, n3.getX(), 0.1);
        assertEquals("n3 on wrong y coordinate", 30 + 10, n3.getY(), 0.1);
        assertEquals("n4 on wrong x coordinate", 60 + 10 + 30, n4.getX(), 0.1);
        assertEquals("n4 on wrong y coordinate", 40, n4.getY(), 0.1);
        assertEquals("n5 on wrong x coordinate", 150, n5.getX(), 0.1);
        assertEquals("n5 on wrong y coordinate", 6, n5.getY(), 0.1);
    }

    @Test
    public void testForwardPromotionMoreNodes() {
        // Desired layout
        // n1 -> n2 -> n3 ->
        //    =======> n4 -> n6 
        //    -------> n5 ->
        //        

        ElkNode parent = ElkGraphUtil.createGraph();
        ElkNode n1 = ElkGraphUtil.createNode(parent);
        ElkNode n2 = ElkGraphUtil.createNode(parent);
        ElkNode n3 = ElkGraphUtil.createNode(parent);
        ElkNode n4 = ElkGraphUtil.createNode(parent);
        ElkNode n5 = ElkGraphUtil.createNode(parent);
        ElkNode n6 = ElkGraphUtil.createNode(parent);
        n1.setDimensions(30, 30);
        n2.setDimensions(30, 30);
        n3.setDimensions(30, 30);
        n4.setDimensions(30, 30);
        n5.setDimensions(30, 30);
        n6.setDimensions(30, 30);
        ElkGraphUtil.createLabel("n1", n1);
        ElkGraphUtil.createLabel("n2", n2);
        ElkGraphUtil.createLabel("n3", n3);
        ElkGraphUtil.createLabel("n4", n4);
        ElkGraphUtil.createLabel("n5", n5);
        ElkGraphUtil.createLabel("n6", n6);

        parent.setProperty(CoreOptions.ALGORITHM, LayeredOptions.ALGORITHM_ID);
        parent.setProperty(CoreOptions.DIRECTION, Direction.RIGHT);
        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY, OrderingStrategy.PREFER_EDGES);
        parent.setProperty(LayeredOptions.LAYERING_STRATEGY, LayeringStrategy.LONGEST_PATH_SOURCE);
        parent.setProperty(LayeredOptions.LAYERING_NODE_PROMOTION_STRATEGY,
                NodePromotionStrategy.MODEL_ORDER_LEFT_TO_RIGHT);
        parent.setProperty(CoreOptions.PADDING, new ElkPadding(0.0));
        parent.setProperty(CoreOptions.SPACING_NODE_NODE, 10.0);
        parent.setProperty(LayeredOptions.SPACING_NODE_NODE_BETWEEN_LAYERS, 20.0);

        ElkGraphUtil.createSimpleEdge(n1, n2);
        ElkGraphUtil.createSimpleEdge(n1, n4);
        ElkGraphUtil.createSimpleEdge(n1, n4);
        ElkGraphUtil.createSimpleEdge(n1, n5);
        ElkGraphUtil.createSimpleEdge(n2, n3);
        ElkGraphUtil.createSimpleEdge(n3, n6);
        ElkGraphUtil.createSimpleEdge(n4, n6);
        ElkGraphUtil.createSimpleEdge(n5, n6);
        
        LayeredLayoutProvider layoutProvider = new LayeredLayoutProvider();
        layoutProvider.layout(parent, new BasicProgressMonitor());

        // Assert the positions or the nodes.
        assertEquals("n1 on wrong x coordinate", 0, n1.getX(), 0.1);
        assertEquals("n1 on wrong y coordinate", 9, n1.getY(), 0.1);
        // Check whether position is node size of n1 + 4*edgeEdge spacing.
        assertEquals("n2 on wrong x coordinate", 30 + 40, n2.getX(), 0.1);
        assertEquals("n2 on wrong y coordinate", 0, n2.getY(), 0.1);
        // Check whether position is node size of n1 + 4*edgeEdge spacing + node size of n2 + between layer spacing.
        assertEquals("n3 on wrong x coordinate", 70 + 30 + 20, n3.getX(), 0.1);
        assertEquals("n3 on wrong y coordinate", 0, n3.getY(), 0.1);
        // Check whether n4 aligns with n3.
        assertEquals("n4 on wrong x coordinate", 120, n4.getX(), 0.1);
        assertEquals("n4 on wrong y coordinate", 30 + 10, n4.getY(), 0.1);
        // Check whether n5 aligns with n4 and n3.
        assertEquals("n5 on wrong x coordinate", 120, n5.getX(), 0.1);
        assertEquals("n5 on wrong y coordinate", 40 + 30 + 10, n5.getY(), 0.1);
        assertEquals("n6 on wrong x coordinate", 170, n6.getX(), 0.1);
        assertEquals("n6 on wrong y coordinate", 40, n6.getY(), 0.1);
    }
    
    @Test
    public void testBackwardPromotionMoreNodes() {
        // Desired layout
        // n1 -> n2 =======>
        //    -> n3 -------> n6 
        //    -> n4 -> n5 ->
        //        

        ElkNode parent = ElkGraphUtil.createGraph();
        ElkNode n1 = ElkGraphUtil.createNode(parent);
        ElkNode n2 = ElkGraphUtil.createNode(parent);
        ElkNode n3 = ElkGraphUtil.createNode(parent);
        ElkNode n4 = ElkGraphUtil.createNode(parent);
        ElkNode n5 = ElkGraphUtil.createNode(parent);
        ElkNode n6 = ElkGraphUtil.createNode(parent);
        n1.setDimensions(30, 30);
        n2.setDimensions(30, 30);
        n3.setDimensions(30, 30);
        n4.setDimensions(30, 30);
        n5.setDimensions(30, 30);
        n6.setDimensions(30, 30);
        ElkGraphUtil.createLabel("n1", n1);
        ElkGraphUtil.createLabel("n2", n2);
        ElkGraphUtil.createLabel("n3", n3);
        ElkGraphUtil.createLabel("n4", n4);
        ElkGraphUtil.createLabel("n5", n5);
        ElkGraphUtil.createLabel("n6", n6);

        parent.setProperty(CoreOptions.ALGORITHM, LayeredOptions.ALGORITHM_ID);
        parent.setProperty(CoreOptions.DIRECTION, Direction.RIGHT);
        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY, OrderingStrategy.PREFER_EDGES);
        parent.setProperty(LayeredOptions.LAYERING_STRATEGY, LayeringStrategy.LONGEST_PATH);
        parent.setProperty(LayeredOptions.LAYERING_NODE_PROMOTION_STRATEGY,
                NodePromotionStrategy.MODEL_ORDER_RIGHT_TO_LEFT);
        parent.setProperty(CoreOptions.PADDING, new ElkPadding(0.0));
        parent.setProperty(CoreOptions.SPACING_NODE_NODE, 10.0);
        parent.setProperty(LayeredOptions.SPACING_NODE_NODE_BETWEEN_LAYERS, 20.0);

        ElkGraphUtil.createSimpleEdge(n1, n2);
        ElkGraphUtil.createSimpleEdge(n1, n3);
        ElkGraphUtil.createSimpleEdge(n1, n4);
        ElkGraphUtil.createSimpleEdge(n2, n6);
        ElkGraphUtil.createSimpleEdge(n2, n6);
        ElkGraphUtil.createSimpleEdge(n3, n6);
        ElkGraphUtil.createSimpleEdge(n4, n5);
        ElkGraphUtil.createSimpleEdge(n5, n6);
        
        LayeredLayoutProvider layoutProvider = new LayeredLayoutProvider();
        layoutProvider.layout(parent, new BasicProgressMonitor());

        // Assert the positions or the nodes.
        assertEquals("n1 on wrong x coordinate", 0, n1.getX(), 0.1);
        assertEquals("n1 on wrong y coordinate", 7.5, n1.getY(), 0.1);
        assertEquals("n2 on wrong x coordinate", 60, n2.getX(), 0.1);
        assertEquals("n2 on wrong y coordinate", 0, n2.getY(), 0.1);
        assertEquals("n3 on wrong x coordinate", 60, n3.getX(), 0.1);
        assertEquals("n3 on wrong y coordinate", 40, n3.getY(), 0.1);
        assertEquals("n4 on wrong x coordinate", 60, n4.getX(), 0.1);
        assertEquals("n4 on wrong y coordinate", 80, n4.getY(), 0.1);
        assertEquals("n5 on wrong x coordinate", 110, n5.getX(), 0.1);
        assertEquals("n5 on wrong y coordinate", 80, n5.getY(), 0.1);
        assertEquals("n6 on wrong x coordinate", 170, n6.getX(), 0.1);
        assertEquals("n6 on wrong y coordinate", 9, n6.getY(), 0.1);
    }

    @Test
    public void testForwardPromotionNoEndNode() {
        // Desired layout
        // n1 -> n2 -> n3
        //    =======> n4
        //

        ElkNode parent = ElkGraphUtil.createGraph();
        ElkNode n1 = ElkGraphUtil.createNode(parent);
        ElkNode n2 = ElkGraphUtil.createNode(parent);
        ElkNode n3 = ElkGraphUtil.createNode(parent);
        ElkNode n4 = ElkGraphUtil.createNode(parent);
        n1.setDimensions(30, 30);
        n2.setDimensions(30, 30);
        n3.setDimensions(30, 30);
        n4.setDimensions(30, 30);
        ElkGraphUtil.createLabel("n1", n1);
        ElkGraphUtil.createLabel("n2", n2);
        ElkGraphUtil.createLabel("n3", n3);
        ElkGraphUtil.createLabel("n4", n4);

        parent.setProperty(CoreOptions.ALGORITHM, LayeredOptions.ALGORITHM_ID);
        parent.setProperty(CoreOptions.DIRECTION, Direction.RIGHT);
        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY, OrderingStrategy.PREFER_EDGES);
        parent.setProperty(LayeredOptions.LAYERING_STRATEGY, LayeringStrategy.LONGEST_PATH_SOURCE);
        parent.setProperty(LayeredOptions.LAYERING_NODE_PROMOTION_STRATEGY,
                NodePromotionStrategy.MODEL_ORDER_LEFT_TO_RIGHT);
        parent.setProperty(CoreOptions.PADDING, new ElkPadding(0.0));
        parent.setProperty(CoreOptions.SPACING_NODE_NODE, 10.0);
        parent.setProperty(LayeredOptions.SPACING_NODE_NODE_BETWEEN_LAYERS, 20.0);

        ElkGraphUtil.createSimpleEdge(n1, n2);
        ElkGraphUtil.createSimpleEdge(n1, n4);
        ElkGraphUtil.createSimpleEdge(n1, n4);
        ElkGraphUtil.createSimpleEdge(n2, n3);
        
        LayeredLayoutProvider layoutProvider = new LayeredLayoutProvider();
        layoutProvider.layout(parent, new BasicProgressMonitor());

        // Assert the positions or the nodes.
        assertEquals("n1 on wrong x coordinate", 0, n1.getX(), 0.1);
        assertEquals("n1 on wrong y coordinate", 7.5, n1.getY(), 0.1);
        // Check whether position is node size of n1 + 4*edgeEdge spacing.
        assertEquals("n2 on wrong x coordinate", 30 + 30, n2.getX(), 0.1);
        assertEquals("n2 on wrong y coordinate", 0, n2.getY(), 0.1);
        // Check whether position is node size of n1 + 4*edgeEdge spacing + node size of n2 + between layer spacing.
        assertEquals("n3 on wrong x coordinate", 60 + 30 + 20, n3.getX(), 0.1);
        assertEquals("n3 on wrong y coordinate", 0, n3.getY(), 0.1);
        // Check whether n4 aligns with n3.
        assertEquals("n4 on wrong x coordinate", 110, n4.getX(), 0.1);
        assertEquals("n4 on wrong y coordinate", 30 + 10, n4.getY(), 0.1);
    }

  @Test
  public void testBackwardPromotionNoEndNode() {
      // Desired layout
      // n1 =======> n4
      // n2 -> n3 ->
      //

      ElkNode parent = ElkGraphUtil.createGraph();
      ElkNode n1 = ElkGraphUtil.createNode(parent);
      ElkNode n2 = ElkGraphUtil.createNode(parent);
      ElkNode n3 = ElkGraphUtil.createNode(parent);
      ElkNode n4 = ElkGraphUtil.createNode(parent);
      n1.setDimensions(30, 30);
      n2.setDimensions(30, 30);
      n3.setDimensions(30, 30);
      n4.setDimensions(30, 30);
      ElkGraphUtil.createLabel("n1", n1);
      ElkGraphUtil.createLabel("n2", n2);
      ElkGraphUtil.createLabel("n3", n3);
      ElkGraphUtil.createLabel("n4", n4);

      parent.setProperty(CoreOptions.ALGORITHM, LayeredOptions.ALGORITHM_ID);
      parent.setProperty(CoreOptions.DIRECTION, Direction.RIGHT);
      parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY, OrderingStrategy.PREFER_EDGES);
      parent.setProperty(LayeredOptions.LAYERING_STRATEGY, LayeringStrategy.LONGEST_PATH);
      parent.setProperty(LayeredOptions.LAYERING_NODE_PROMOTION_STRATEGY,
              NodePromotionStrategy.MODEL_ORDER_RIGHT_TO_LEFT);
      parent.setProperty(CoreOptions.PADDING, new ElkPadding(0.0));
      parent.setProperty(CoreOptions.SPACING_NODE_NODE, 10.0);
      parent.setProperty(LayeredOptions.SPACING_NODE_NODE_BETWEEN_LAYERS, 20.0);

      ElkGraphUtil.createSimpleEdge(n1, n4);
      ElkGraphUtil.createSimpleEdge(n1, n4);
      ElkGraphUtil.createSimpleEdge(n2, n3);
      ElkGraphUtil.createSimpleEdge(n3, n4);
      
      LayeredLayoutProvider layoutProvider = new LayeredLayoutProvider();
      layoutProvider.layout(parent, new BasicProgressMonitor());

      // Assert the positions or the nodes.
      assertEquals("n1 on wrong x coordinate", 0, n1.getX(), 0.1);
      assertEquals("n1 on wrong y coordinate", 0, n1.getY(), 0.1);
      assertEquals("n2 on wrong x coordinate", 0, n2.getX(), 0.1);
      assertEquals("n2 on wrong y coordinate", 40, n2.getY(), 0.1);
      assertEquals("n3 on wrong x coordinate", 50, n3.getX(), 0.1);
      assertEquals("n3 on wrong y coordinate", 40, n3.getY(), 0.1);
      assertEquals("n4 on wrong x coordinate", 100, n4.getX(), 0.1);
      assertEquals("n4 on wrong y coordinate", 6, n4.getY(), 0.1);
  }

  @Test
  public void testBackwardPromotionNoEndNode2() {
      // Desired layout
      // n1 -> n2
      //    -> n3 -> n4 
      //

      ElkNode parent = ElkGraphUtil.createGraph();
      ElkNode n1 = ElkGraphUtil.createNode(parent);
      ElkNode n2 = ElkGraphUtil.createNode(parent);
      ElkNode n3 = ElkGraphUtil.createNode(parent);
      ElkNode n4 = ElkGraphUtil.createNode(parent);
      n1.setDimensions(30, 30);
      n2.setDimensions(30, 30);
      n3.setDimensions(30, 30);
      n4.setDimensions(30, 30);
      ElkGraphUtil.createLabel("n1", n1);
      ElkGraphUtil.createLabel("n2", n2);
      ElkGraphUtil.createLabel("n3", n3);
      ElkGraphUtil.createLabel("n4", n4);

      parent.setProperty(CoreOptions.ALGORITHM, LayeredOptions.ALGORITHM_ID);
      parent.setProperty(CoreOptions.DIRECTION, Direction.RIGHT);
      parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY, OrderingStrategy.PREFER_EDGES);
      parent.setProperty(LayeredOptions.LAYERING_STRATEGY, LayeringStrategy.LONGEST_PATH);
      parent.setProperty(LayeredOptions.LAYERING_NODE_PROMOTION_STRATEGY,
              NodePromotionStrategy.MODEL_ORDER_RIGHT_TO_LEFT);
      parent.setProperty(CoreOptions.PADDING, new ElkPadding(0.0));
      parent.setProperty(CoreOptions.SPACING_NODE_NODE, 10.0);
      parent.setProperty(LayeredOptions.SPACING_NODE_NODE_BETWEEN_LAYERS, 20.0);

      ElkGraphUtil.createSimpleEdge(n1, n2);
      ElkGraphUtil.createSimpleEdge(n1, n3);
      ElkGraphUtil.createSimpleEdge(n3, n4);
      
      LayeredLayoutProvider layoutProvider = new LayeredLayoutProvider();
      layoutProvider.layout(parent, new BasicProgressMonitor());

      // Assert the positions or the nodes.
      assertEquals("n1 on wrong x coordinate", 0, n1.getX(), 0.1);
      assertEquals("n1 on wrong y coordinate", 5, n1.getY(), 0.1);
      assertEquals("n2 on wrong x coordinate", 50, n2.getX(), 0.1);
      assertEquals("n2 on wrong y coordinate", 0, n2.getY(), 0.1);
      assertEquals("n3 on wrong x coordinate", 50, n3.getX(), 0.1);
      assertEquals("n3 on wrong y coordinate", 40, n3.getY(), 0.1);
      assertEquals("n4 on wrong x coordinate", 100, n4.getX(), 0.1);
      assertEquals("n4 on wrong y coordinate", 40, n4.getY(), 0.1);
  }

}
