/*******************************************************************************
 * Copyright (c) 2024 Kiel University.
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
import org.eclipse.elk.alg.layered.options.CrossingMinimizationStrategy;
import org.eclipse.elk.alg.layered.options.CycleBreakingStrategy;
import org.eclipse.elk.alg.layered.options.GreedySwitchType;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.layered.options.OrderingStrategy;
import org.eclipse.elk.alg.test.PlainJavaInitialization;
import org.eclipse.elk.core.math.ElkPadding;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.options.Direction;
import org.eclipse.elk.core.options.PortConstraints;
import org.eclipse.elk.core.options.PortSide;
import org.eclipse.elk.core.util.BasicProgressMonitor;
import org.eclipse.elk.core.util.NullElkProgressMonitor;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.ElkPort;
import org.eclipse.elk.graph.util.ElkGraphUtil;
import org.junit.BeforeClass;
import org.junit.Test;

/**
 * Test whether these concrete example did not break while using the consider model order strategy.
 *
 */
public class ConcreteSortByInputModelTest {

    @BeforeClass
    public static void init() {
        PlainJavaInitialization.initializePlainJavaLayout();
    }

    @Test
    public void testOutgoingEastIncomingWest() {
        //   p1--->p4
        // n1        n2
        //   p2<---p3

        ElkNode parent = ElkGraphUtil.createGraph();
        ElkNode n1 = ElkGraphUtil.createNode(parent);
        ElkNode n2 = ElkGraphUtil.createNode(parent);
        n1.setDimensions(20, 20);
        n2.setDimensions(20, 20);

        parent.setProperty(CoreOptions.ALGORITHM, LayeredOptions.ALGORITHM_ID);
        // Model order configuration.
        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY, OrderingStrategy.PREFER_EDGES);
        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_PORT_MODEL_ORDER, false);
        // No crossing minimization, hence only model order is used.
        parent.setProperty(LayeredOptions.CROSSING_MINIMIZATION_STRATEGY, CrossingMinimizationStrategy.NONE);
        parent.setProperty(LayeredOptions.CROSSING_MINIMIZATION_GREEDY_SWITCH_TYPE, GreedySwitchType.OFF);
        // For presentation.
        parent.setProperty(CoreOptions.DIRECTION, Direction.RIGHT);
        parent.setProperty(CoreOptions.PADDING, new ElkPadding(0.0));
        parent.setProperty(CoreOptions.SPACING_NODE_NODE, 10.0);
        parent.setProperty(LayeredOptions.SPACING_NODE_NODE_BETWEEN_LAYERS, 20.0);
        n1.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        n2.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        // Make sure that the order is correct.
        parent.setProperty(LayeredOptions.CYCLE_BREAKING_STRATEGY, CycleBreakingStrategy.MODEL_ORDER);
        
        ElkPort p1 = ElkGraphUtil.createPort(n1);
        p1.setProperty(LayeredOptions.PORT_SIDE, PortSide.EAST);
        ElkPort p2 = ElkGraphUtil.createPort(n1);
        p2.setProperty(LayeredOptions.PORT_SIDE, PortSide.EAST);
        ElkPort p3 = ElkGraphUtil.createPort(n2);
        p3.setProperty(LayeredOptions.PORT_SIDE, PortSide.WEST);
        ElkPort p4 = ElkGraphUtil.createPort(n2);
        p4.setProperty(LayeredOptions.PORT_SIDE, PortSide.WEST);

        ElkGraphUtil.createSimpleEdge(p1, p4);
        ElkGraphUtil.createSimpleEdge(p2, p3);
        
        LayeredLayoutProvider layoutProvider = new LayeredLayoutProvider();
        layoutProvider.layout(parent, new BasicProgressMonitor());
        // Assert the ordering of p1 and p2.
        assertTrue("p1 above p2", p1.getY() < p2.getY());
        // Assert the ordering of p3 and p4. 
        assertTrue("p4 above p3", p4.getY() < p3.getY());
        
        // Check for port model order, which should result in the same graph.

        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_PORT_MODEL_ORDER, true);
        
        layoutProvider.layout(parent, new BasicProgressMonitor());
        // Assert the ordering of p1 and p2.
        assertTrue("p1 above p2", p1.getY() < p2.getY());
        // Assert the ordering of p3 and p4. 
        assertTrue("p4 above p3", p4.getY() < p3.getY());
        
    }

    @Test
    public void testOutgoingEastIncomingWestMultipleNodes() {
        // n1:p1--->p4
        //            n2
        //       |--p3
        // ni:p2<-

        ElkNode parent = ElkGraphUtil.createGraph();
        ElkNode n1 = ElkGraphUtil.createNode(parent);
        ElkNode ni = ElkGraphUtil.createNode(parent);
        ElkNode n2 = ElkGraphUtil.createNode(parent);
        n1.setDimensions(20, 20);
        ni.setDimensions(20, 20);
        n2.setDimensions(20, 20);

        parent.setProperty(CoreOptions.ALGORITHM, LayeredOptions.ALGORITHM_ID);
        // Model order configuration.
        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY, OrderingStrategy.PREFER_EDGES);
        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_PORT_MODEL_ORDER, false);
        // No crossing minimization, hence only model order is used.
        parent.setProperty(LayeredOptions.CROSSING_MINIMIZATION_STRATEGY, CrossingMinimizationStrategy.NONE);
        parent.setProperty(LayeredOptions.CROSSING_MINIMIZATION_GREEDY_SWITCH_TYPE, GreedySwitchType.OFF);
        // For presentation.
        parent.setProperty(CoreOptions.DIRECTION, Direction.RIGHT);
        parent.setProperty(CoreOptions.PADDING, new ElkPadding(0.0));
        parent.setProperty(CoreOptions.SPACING_NODE_NODE, 10.0);
        parent.setProperty(LayeredOptions.SPACING_NODE_NODE_BETWEEN_LAYERS, 20.0);
        n1.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        ni.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        n2.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        // Make sure that the order is correct.
        parent.setProperty(LayeredOptions.CYCLE_BREAKING_STRATEGY, CycleBreakingStrategy.MODEL_ORDER);
        
        ElkPort p1 = ElkGraphUtil.createPort(n1);
        p1.setProperty(LayeredOptions.PORT_SIDE, PortSide.EAST);
        ElkPort p2 = ElkGraphUtil.createPort(ni);
        p2.setProperty(LayeredOptions.PORT_SIDE, PortSide.EAST);
        ElkPort p3 = ElkGraphUtil.createPort(n2);
        p3.setProperty(LayeredOptions.PORT_SIDE, PortSide.WEST);
        ElkPort p4 = ElkGraphUtil.createPort(n2);
        p4.setProperty(LayeredOptions.PORT_SIDE, PortSide.WEST);

        ElkGraphUtil.createSimpleEdge(p1, p4);
        ElkGraphUtil.createSimpleEdge(p2, p3);
        
        LayeredLayoutProvider layoutProvider = new LayeredLayoutProvider();
        layoutProvider.layout(parent, new BasicProgressMonitor());
        // Assert the ordering of n1 and ni.
        assertTrue("Node order of n1 and ni", n1.getY() < ni.getY());
        // Assert the ordering of p3 and p4. 
        assertTrue("p4 above p3", p4.getY() < p3.getY());
        
        // Check for port model order, which should result in the same graph.

        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_PORT_MODEL_ORDER, true);
        
        layoutProvider.layout(parent, new BasicProgressMonitor());
        // Assert the ordering of n1 and ni.
        assertTrue("Node order of n1 and ni", n1.getY() < ni.getY());
        // Assert the ordering of p3 and p4. 
        assertTrue("p4 above p3", p4.getY() < p3.getY());
    }

    @Test
    public void testOutgoingEastIncomingNorth() {
        //   p1---------
        // n1           |
        //   p2<----    |
        //         |    v
        //         p3  p4
        //           n2

        ElkNode parent = ElkGraphUtil.createGraph();
        ElkNode n1 = ElkGraphUtil.createNode(parent);
        ElkNode n2 = ElkGraphUtil.createNode(parent);
        n1.setDimensions(20, 20);
        n2.setDimensions(20, 20);

        parent.setProperty(CoreOptions.ALGORITHM, LayeredOptions.ALGORITHM_ID);
        // Model order configuration.
        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY, OrderingStrategy.PREFER_EDGES);
        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_PORT_MODEL_ORDER, false);
        // No crossing minimization, hence only model order is used.
        parent.setProperty(LayeredOptions.CROSSING_MINIMIZATION_STRATEGY, CrossingMinimizationStrategy.NONE);
        parent.setProperty(LayeredOptions.CROSSING_MINIMIZATION_GREEDY_SWITCH_TYPE, GreedySwitchType.OFF);
        // For presentation.
        parent.setProperty(CoreOptions.DIRECTION, Direction.RIGHT);
        parent.setProperty(CoreOptions.PADDING, new ElkPadding(0.0));
        parent.setProperty(CoreOptions.SPACING_NODE_NODE, 10.0);
        parent.setProperty(LayeredOptions.SPACING_NODE_NODE_BETWEEN_LAYERS, 20.0);
        n1.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        n2.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        // Make sure that the order is correct.
        parent.setProperty(LayeredOptions.CYCLE_BREAKING_STRATEGY, CycleBreakingStrategy.MODEL_ORDER);
        
        ElkPort p1 = ElkGraphUtil.createPort(n1);
        p1.setProperty(LayeredOptions.PORT_SIDE, PortSide.EAST);
        ElkPort p2 = ElkGraphUtil.createPort(n1);
        p2.setProperty(LayeredOptions.PORT_SIDE, PortSide.EAST);
        ElkPort p3 = ElkGraphUtil.createPort(n2);
        p3.setProperty(LayeredOptions.PORT_SIDE, PortSide.NORTH);
        ElkPort p4 = ElkGraphUtil.createPort(n2);
        p4.setProperty(LayeredOptions.PORT_SIDE, PortSide.NORTH);

        ElkGraphUtil.createSimpleEdge(p1, p4);
        ElkGraphUtil.createSimpleEdge(p2, p3);
        
        LayeredLayoutProvider layoutProvider = new LayeredLayoutProvider();
        layoutProvider.layout(parent, new BasicProgressMonitor());
        // Assert the ordering of p1 and p2.
        assertTrue("p1 above p2", p1.getY() < p2.getY());
        // Assert the ordering of p3 and p4. 
        assertTrue("p3 before p4", p3.getX() < p4.getX());
        
        // Check for port model order, which should result in the same graph.

        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_PORT_MODEL_ORDER, true);
        
        layoutProvider.layout(parent, new BasicProgressMonitor());
        // Assert the ordering of p1 and p2.
        assertTrue("p1 above p2", p1.getY() < p2.getY());
        // Assert the ordering of p3 and p4. 
        assertTrue("p3 before p4", p3.getX() < p4.getX());
    }

    @Test
    public void testOutgoingEastIncomingNorthMultipleNodes() {
        //n1:p1---------
        //              |
        //ni:p2<----    |
        //         |    v
        //         p3  p4
        //           n2

        ElkNode parent = ElkGraphUtil.createGraph();
        ElkNode n1 = ElkGraphUtil.createNode(parent);
        ElkNode ni = ElkGraphUtil.createNode(parent);
        ElkNode n2 = ElkGraphUtil.createNode(parent);
        n1.setDimensions(20, 20);
        ni.setDimensions(20, 20);
        n2.setDimensions(20, 20);

        parent.setProperty(CoreOptions.ALGORITHM, LayeredOptions.ALGORITHM_ID);
        // Model order configuration.
        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY, OrderingStrategy.PREFER_EDGES);
        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_PORT_MODEL_ORDER, false);
        // No crossing minimization, hence only model order is used.
        parent.setProperty(LayeredOptions.CROSSING_MINIMIZATION_STRATEGY, CrossingMinimizationStrategy.NONE);
        parent.setProperty(LayeredOptions.CROSSING_MINIMIZATION_GREEDY_SWITCH_TYPE, GreedySwitchType.OFF);
        // For presentation.
        parent.setProperty(CoreOptions.DIRECTION, Direction.RIGHT);
        parent.setProperty(CoreOptions.PADDING, new ElkPadding(0.0));
        parent.setProperty(CoreOptions.SPACING_NODE_NODE, 10.0);
        parent.setProperty(LayeredOptions.SPACING_NODE_NODE_BETWEEN_LAYERS, 20.0);
        n1.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        ni.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        n2.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        // Make sure that the order is correct.
        parent.setProperty(LayeredOptions.CYCLE_BREAKING_STRATEGY, CycleBreakingStrategy.MODEL_ORDER);
        
        ElkPort p1 = ElkGraphUtil.createPort(n1);
        p1.setProperty(LayeredOptions.PORT_SIDE, PortSide.EAST);
        ElkPort p2 = ElkGraphUtil.createPort(ni);
        p2.setProperty(LayeredOptions.PORT_SIDE, PortSide.EAST);
        ElkPort p3 = ElkGraphUtil.createPort(n2);
        p3.setProperty(LayeredOptions.PORT_SIDE, PortSide.NORTH);
        ElkPort p4 = ElkGraphUtil.createPort(n2);
        p4.setProperty(LayeredOptions.PORT_SIDE, PortSide.NORTH);

        ElkGraphUtil.createSimpleEdge(p1, p4);
        ElkGraphUtil.createSimpleEdge(p2, p3);
        
        LayeredLayoutProvider layoutProvider = new LayeredLayoutProvider();
        layoutProvider.layout(parent, new BasicProgressMonitor());
        // Assert the ordering of n1 and ni.
        assertTrue("Node order of n1 and ni", n1.getY() < ni.getY());
        // Assert the ordering of p3 and p4. 
        assertTrue("p3 before p4", p3.getX() < p4.getX());
        
        // Check for port model order, which should result in the same graph.

        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_PORT_MODEL_ORDER, true);
        
        layoutProvider.layout(parent, new BasicProgressMonitor());
        // Assert the ordering of n1 and ni.
        assertTrue("Node order of n1 and ni", n1.getY() < ni.getY());
        // Assert the ordering of p3 and p4. 
        assertTrue("p3 before p4", p3.getX() < p4.getX());
    }

    @Test
    public void testOutgoingEastIncomingSouth() {
        //           n2
        //         p4  p3
        //         A    |
        //   p1-----    |
        // n1           |
        //   p2<--------

        ElkNode parent = ElkGraphUtil.createGraph();
        ElkNode n1 = ElkGraphUtil.createNode(parent);
        ElkNode n2 = ElkGraphUtil.createNode(parent);
        n1.setDimensions(20, 20);
        n2.setDimensions(20, 20);

        parent.setProperty(CoreOptions.ALGORITHM, LayeredOptions.ALGORITHM_ID);
        // Model order configuration.
        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY, OrderingStrategy.PREFER_EDGES);
        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_PORT_MODEL_ORDER, false);
        // No crossing minimization, hence only model order is used.
        parent.setProperty(LayeredOptions.CROSSING_MINIMIZATION_STRATEGY, CrossingMinimizationStrategy.NONE);
        parent.setProperty(LayeredOptions.CROSSING_MINIMIZATION_GREEDY_SWITCH_TYPE, GreedySwitchType.OFF);
        // For presentation.
        parent.setProperty(CoreOptions.DIRECTION, Direction.RIGHT);
        parent.setProperty(CoreOptions.PADDING, new ElkPadding(0.0));
        parent.setProperty(CoreOptions.SPACING_NODE_NODE, 10.0);
        parent.setProperty(LayeredOptions.SPACING_NODE_NODE_BETWEEN_LAYERS, 20.0);
        n1.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        n2.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        // Make sure that the order is correct.
        parent.setProperty(LayeredOptions.CYCLE_BREAKING_STRATEGY, CycleBreakingStrategy.MODEL_ORDER);
        
        ElkPort p1 = ElkGraphUtil.createPort(n1);
        p1.setProperty(LayeredOptions.PORT_SIDE, PortSide.EAST);
        ElkPort p2 = ElkGraphUtil.createPort(n1);
        p2.setProperty(LayeredOptions.PORT_SIDE, PortSide.EAST);
        ElkPort p3 = ElkGraphUtil.createPort(n2);
        p3.setProperty(LayeredOptions.PORT_SIDE, PortSide.SOUTH);
        ElkPort p4 = ElkGraphUtil.createPort(n2);
        p4.setProperty(LayeredOptions.PORT_SIDE, PortSide.SOUTH);

        ElkGraphUtil.createSimpleEdge(p1, p4);
        ElkGraphUtil.createSimpleEdge(p2, p3);
        
        LayeredLayoutProvider layoutProvider = new LayeredLayoutProvider();
        layoutProvider.layout(parent, new BasicProgressMonitor());
        // Assert the ordering of p1 and p2.
        assertTrue("p1 above p2", p1.getY() < p2.getY());
        // Assert the ordering of p3 and p4. 
        assertTrue("p4 before p3", p4.getX() < p3.getX());
        
        // Check for port model order, which should result in the same graph.

        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_PORT_MODEL_ORDER, true);
        
        layoutProvider.layout(parent, new BasicProgressMonitor());
        // Assert the ordering of p1 and p2.
        assertTrue("p1 above p2", p1.getY() < p2.getY());
        // Assert the ordering of p3 and p4. 
        assertTrue("p4 before p3", p4.getX() < p3.getX());
    }

    @Test
    public void testOutgoingEastIncomingSouthMultipleNodes() {
        //           n2
        //         p4  p3
        //         A    |
        //n1:p1-----    |
        //              |
        //ni:p2<--------

        ElkNode parent = ElkGraphUtil.createGraph();
        ElkNode n1 = ElkGraphUtil.createNode(parent);
        ElkNode ni = ElkGraphUtil.createNode(parent);
        ElkNode n2 = ElkGraphUtil.createNode(parent);
        n1.setDimensions(20, 20);
        ni.setDimensions(20, 20);
        n2.setDimensions(20, 20);

        parent.setProperty(CoreOptions.ALGORITHM, LayeredOptions.ALGORITHM_ID);
        // Model order configuration.
        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY, OrderingStrategy.PREFER_EDGES);
        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_PORT_MODEL_ORDER, false);
        // No crossing minimization, hence only model order is used.
        parent.setProperty(LayeredOptions.CROSSING_MINIMIZATION_STRATEGY, CrossingMinimizationStrategy.NONE);
        parent.setProperty(LayeredOptions.CROSSING_MINIMIZATION_GREEDY_SWITCH_TYPE, GreedySwitchType.OFF);
        // For presentation.
        parent.setProperty(CoreOptions.DIRECTION, Direction.RIGHT);
        parent.setProperty(CoreOptions.PADDING, new ElkPadding(0.0));
        parent.setProperty(CoreOptions.SPACING_NODE_NODE, 10.0);
        parent.setProperty(LayeredOptions.SPACING_NODE_NODE_BETWEEN_LAYERS, 20.0);
        n1.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        ni.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        n2.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        // Make sure that the order is correct.
        parent.setProperty(LayeredOptions.CYCLE_BREAKING_STRATEGY, CycleBreakingStrategy.MODEL_ORDER);
        
        ElkPort p1 = ElkGraphUtil.createPort(n1);
        p1.setProperty(LayeredOptions.PORT_SIDE, PortSide.EAST);
        ElkPort p2 = ElkGraphUtil.createPort(ni);
        p2.setProperty(LayeredOptions.PORT_SIDE, PortSide.EAST);
        ElkPort p3 = ElkGraphUtil.createPort(n2);
        p3.setProperty(LayeredOptions.PORT_SIDE, PortSide.SOUTH);
        ElkPort p4 = ElkGraphUtil.createPort(n2);
        p4.setProperty(LayeredOptions.PORT_SIDE, PortSide.SOUTH);

        ElkGraphUtil.createSimpleEdge(p1, p4);
        ElkGraphUtil.createSimpleEdge(p2, p3);
        
        LayeredLayoutProvider layoutProvider = new LayeredLayoutProvider();
        layoutProvider.layout(parent, new BasicProgressMonitor());
        // Assert the ordering of n1 and ni.
        assertTrue("Node order of n1 and ni", n1.getY() < ni.getY());
        // Assert the ordering of p3 and p4. 
        assertTrue("p4 before p3", p4.getX() < p3.getX());
        
        // Check for port model order, which should result in the same graph.

        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_PORT_MODEL_ORDER, true);

        layoutProvider.layout(parent, new BasicProgressMonitor());
        // Assert the ordering of n1 and ni.
        assertTrue("Node order of n1 and ni", n1.getY() < ni.getY());
        // Assert the ordering of p3 and p4. 
        assertTrue("p4 before p3", p4.getX() < p3.getX());
    }

    @Test
    public void testOutgoingEastIncomingEast() {
        //        p3<----
        //      n2      |
        //        p4<-  |
        //           |  |
        //           |  |
        //   p1------   |
        // n1           |
        //   p2<--------

        ElkNode parent = ElkGraphUtil.createGraph();
        ElkNode n1 = ElkGraphUtil.createNode(parent);
        ElkNode n2 = ElkGraphUtil.createNode(parent);
        n1.setDimensions(20, 20);
        n2.setDimensions(20, 20);

        parent.setProperty(CoreOptions.ALGORITHM, LayeredOptions.ALGORITHM_ID);
        // Model order configuration.
        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY, OrderingStrategy.PREFER_EDGES);
        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_PORT_MODEL_ORDER, false);
        // No crossing minimization, hence only model order is used.
        parent.setProperty(LayeredOptions.CROSSING_MINIMIZATION_STRATEGY, CrossingMinimizationStrategy.NONE);
        parent.setProperty(LayeredOptions.CROSSING_MINIMIZATION_GREEDY_SWITCH_TYPE, GreedySwitchType.OFF);
        // For presentation.
        parent.setProperty(CoreOptions.DIRECTION, Direction.RIGHT);
        parent.setProperty(CoreOptions.PADDING, new ElkPadding(0.0));
        parent.setProperty(CoreOptions.SPACING_NODE_NODE, 10.0);
        parent.setProperty(LayeredOptions.SPACING_NODE_NODE_BETWEEN_LAYERS, 20.0);
        n1.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        n2.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        // Make sure that the order is correct.
        parent.setProperty(LayeredOptions.CYCLE_BREAKING_STRATEGY, CycleBreakingStrategy.MODEL_ORDER);
        
        ElkPort p1 = ElkGraphUtil.createPort(n1);
        p1.setProperty(LayeredOptions.PORT_SIDE, PortSide.EAST);
        ElkPort p2 = ElkGraphUtil.createPort(n1);
        p2.setProperty(LayeredOptions.PORT_SIDE, PortSide.EAST);
        ElkPort p3 = ElkGraphUtil.createPort(n2);
        p3.setProperty(LayeredOptions.PORT_SIDE, PortSide.EAST);
        ElkPort p4 = ElkGraphUtil.createPort(n2);
        p4.setProperty(LayeredOptions.PORT_SIDE, PortSide.EAST);

        ElkGraphUtil.createSimpleEdge(p1, p4);
        ElkGraphUtil.createSimpleEdge(p2, p3);
        
        LayeredLayoutProvider layoutProvider = new LayeredLayoutProvider();
        layoutProvider.layout(parent, new BasicProgressMonitor());
        // Assert the ordering of p1 and p2.
        assertTrue("p1 above p2", p1.getY() < p2.getY());
        // Assert the ordering of p3 and p4. 
        assertTrue("p3 above p4", p3.getY() < p4.getY());
        
        // Check for port model order, which should result in the same graph.

        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_PORT_MODEL_ORDER, true);

        layoutProvider.layout(parent, new BasicProgressMonitor());
        // Assert the ordering of p1 and p2.
        assertTrue("p1 above p2", p1.getY() < p2.getY());
        // Assert the ordering of p3 and p4. 
        assertTrue("p3 above p4", p3.getY() < p4.getY());
    }

    @Test
    public void testOutgoingEastIncomingEastMultipleNodes() {
        //        p3<----
        //      n2      |
        //        p4<-  |
        //           |  |
        //           |  |
        //n1:p1------   |
        //              |
        //ni:p2<--------

        ElkNode parent = ElkGraphUtil.createGraph();
        ElkNode n1 = ElkGraphUtil.createNode(parent);
        ElkNode ni = ElkGraphUtil.createNode(parent);
        ElkNode n2 = ElkGraphUtil.createNode(parent);
        n1.setDimensions(20, 20);
        ni.setDimensions(20, 20);
        n2.setDimensions(20, 20);

        parent.setProperty(CoreOptions.ALGORITHM, LayeredOptions.ALGORITHM_ID);
        // Model order configuration.
        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY, OrderingStrategy.PREFER_EDGES);
        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_PORT_MODEL_ORDER, false);
        // No crossing minimization, hence only model order is used.
        parent.setProperty(LayeredOptions.CROSSING_MINIMIZATION_STRATEGY, CrossingMinimizationStrategy.NONE);
        parent.setProperty(LayeredOptions.CROSSING_MINIMIZATION_GREEDY_SWITCH_TYPE, GreedySwitchType.OFF);
        // For presentation.
        parent.setProperty(CoreOptions.DIRECTION, Direction.RIGHT);
        parent.setProperty(CoreOptions.PADDING, new ElkPadding(0.0));
        parent.setProperty(CoreOptions.SPACING_NODE_NODE, 10.0);
        parent.setProperty(LayeredOptions.SPACING_NODE_NODE_BETWEEN_LAYERS, 20.0);
        n1.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        ni.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        n2.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        // Make sure that the order is correct.
        parent.setProperty(LayeredOptions.CYCLE_BREAKING_STRATEGY, CycleBreakingStrategy.MODEL_ORDER);
        
        ElkPort p1 = ElkGraphUtil.createPort(n1);
        p1.setProperty(LayeredOptions.PORT_SIDE, PortSide.EAST);
        ElkPort p2 = ElkGraphUtil.createPort(ni);
        p2.setProperty(LayeredOptions.PORT_SIDE, PortSide.EAST);
        ElkPort p3 = ElkGraphUtil.createPort(n2);
        p3.setProperty(LayeredOptions.PORT_SIDE, PortSide.EAST);
        ElkPort p4 = ElkGraphUtil.createPort(n2);
        p4.setProperty(LayeredOptions.PORT_SIDE, PortSide.EAST);

        ElkGraphUtil.createSimpleEdge(p1, p4);
        ElkGraphUtil.createSimpleEdge(p2, p3);
        
        LayeredLayoutProvider layoutProvider = new LayeredLayoutProvider();
        layoutProvider.layout(parent, new BasicProgressMonitor());
        // Assert the ordering of n1 and ni.
        assertTrue("Node order of n1 and ni", n1.getY() < ni.getY());
        // Assert the ordering of p3 and p4. 
        assertTrue("p3 above p4", p3.getY() < p4.getY());
        
        // Check for port model order, which should result in the same graph.

        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_PORT_MODEL_ORDER, true);

        layoutProvider.layout(parent, new BasicProgressMonitor());
        // Assert the ordering of n1 and ni.
        assertTrue("Node order of n1 and ni", n1.getY() < ni.getY());
        // Assert the ordering of p3 and p4. 
        assertTrue("p3 above p4", p3.getY() < p4.getY());
    }

    @Test
    public void testOutgoingNorthIncomingWest() {
        //  ------>p4
        // |         n2
        // |  |----p3
        // |  v
        //p1  p2
        //  n1

        ElkNode parent = ElkGraphUtil.createGraph();
        ElkNode n1 = ElkGraphUtil.createNode(parent);
        ElkNode n2 = ElkGraphUtil.createNode(parent);
        n1.setDimensions(20, 20);
        n2.setDimensions(20, 20);

        parent.setProperty(CoreOptions.ALGORITHM, LayeredOptions.ALGORITHM_ID);
        // Model order configuration.
        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY, OrderingStrategy.PREFER_EDGES);
        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_PORT_MODEL_ORDER, false);
        // No crossing minimization, hence only model order is used.
        parent.setProperty(LayeredOptions.CROSSING_MINIMIZATION_STRATEGY, CrossingMinimizationStrategy.NONE);
        parent.setProperty(LayeredOptions.CROSSING_MINIMIZATION_GREEDY_SWITCH_TYPE, GreedySwitchType.OFF);
        // For presentation.
        parent.setProperty(CoreOptions.DIRECTION, Direction.RIGHT);
        parent.setProperty(CoreOptions.PADDING, new ElkPadding(0.0));
        parent.setProperty(CoreOptions.SPACING_NODE_NODE, 10.0);
        parent.setProperty(LayeredOptions.SPACING_NODE_NODE_BETWEEN_LAYERS, 20.0);
        n1.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        n2.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        // Make sure that the order is correct.
        parent.setProperty(LayeredOptions.CYCLE_BREAKING_STRATEGY, CycleBreakingStrategy.MODEL_ORDER);
        
        ElkPort p1 = ElkGraphUtil.createPort(n1);
        p1.setProperty(LayeredOptions.PORT_SIDE, PortSide.NORTH);
        ElkPort p2 = ElkGraphUtil.createPort(n1);
        p2.setProperty(LayeredOptions.PORT_SIDE, PortSide.NORTH);
        ElkPort p3 = ElkGraphUtil.createPort(n2);
        p3.setProperty(LayeredOptions.PORT_SIDE, PortSide.WEST);
        ElkPort p4 = ElkGraphUtil.createPort(n2);
        p4.setProperty(LayeredOptions.PORT_SIDE, PortSide.WEST);

        ElkGraphUtil.createSimpleEdge(p1, p4);
        ElkGraphUtil.createSimpleEdge(p2, p3);
        
        LayeredLayoutProvider layoutProvider = new LayeredLayoutProvider();
        layoutProvider.layout(parent, new BasicProgressMonitor());
        // Assert the ordering of p1 and p2.
        assertTrue("p1 before p2", p1.getX() < p2.getX());
        // Assert the ordering of p3 and p4. 
        assertTrue("p4 above p3", p4.getY() < p3.getY());
        
        // Check for port model order, which should result in the same graph.

        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_PORT_MODEL_ORDER, true);

        layoutProvider.layout(parent, new BasicProgressMonitor());
        // Assert the ordering of p1 and p2.
        assertTrue("p1 before p2", p1.getX() < p2.getX());
        // Assert the ordering of p3 and p4. 
        assertTrue("p4 above p3", p4.getY() < p3.getY());
    }

    @Test
    public void testOutgoingNorthIncomingWestMultipleNodes() {
        //  |--|
        //  |  |
        // p1  |
        // n1  |
        //     |--->p4
        //            n2
        //  |-------p3
        //  v
        // p2
        // ni

        ElkNode parent = ElkGraphUtil.createGraph();
        ElkNode n1 = ElkGraphUtil.createNode(parent);
        ElkNode ni = ElkGraphUtil.createNode(parent);
        ElkNode n2 = ElkGraphUtil.createNode(parent);
        n1.setDimensions(20, 20);
        ni.setDimensions(20, 20);
        n2.setDimensions(20, 20);

        parent.setProperty(CoreOptions.ALGORITHM, LayeredOptions.ALGORITHM_ID);
        // Model order configuration.
        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY, OrderingStrategy.PREFER_EDGES);
        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_PORT_MODEL_ORDER, false);
        // No crossing minimization, hence only model order is used.
        parent.setProperty(LayeredOptions.CROSSING_MINIMIZATION_STRATEGY, CrossingMinimizationStrategy.NONE);
        parent.setProperty(LayeredOptions.CROSSING_MINIMIZATION_GREEDY_SWITCH_TYPE, GreedySwitchType.OFF);
        // For presentation.
        parent.setProperty(CoreOptions.DIRECTION, Direction.RIGHT);
        parent.setProperty(CoreOptions.PADDING, new ElkPadding(0.0));
        parent.setProperty(CoreOptions.SPACING_NODE_NODE, 10.0);
        parent.setProperty(LayeredOptions.SPACING_NODE_NODE_BETWEEN_LAYERS, 20.0);
        n1.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        ni.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        n2.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        // Make sure that the order is correct.
        parent.setProperty(LayeredOptions.CYCLE_BREAKING_STRATEGY, CycleBreakingStrategy.MODEL_ORDER);
        
        ElkPort p1 = ElkGraphUtil.createPort(n1);
        p1.setProperty(LayeredOptions.PORT_SIDE, PortSide.NORTH);
        ElkPort p2 = ElkGraphUtil.createPort(ni);
        p2.setProperty(LayeredOptions.PORT_SIDE, PortSide.NORTH);
        ElkPort p3 = ElkGraphUtil.createPort(n2);
        p3.setProperty(LayeredOptions.PORT_SIDE, PortSide.WEST);
        ElkPort p4 = ElkGraphUtil.createPort(n2);
        p4.setProperty(LayeredOptions.PORT_SIDE, PortSide.WEST);

        ElkGraphUtil.createSimpleEdge(p1, p4);
        ElkGraphUtil.createSimpleEdge(p2, p3);
        
        LayeredLayoutProvider layoutProvider = new LayeredLayoutProvider();
        layoutProvider.layout(parent, new BasicProgressMonitor());
        // Assert the ordering of n1 and ni.
        assertTrue("Node order of n1 and ni", n1.getY() < ni.getY());
        // Assert the ordering of p3 and p4. 
        assertTrue("p4 above p3", p4.getY() < p3.getY());
        
        // Check for port model order, which should result in the same graph.

        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_PORT_MODEL_ORDER, true);

        layoutProvider.layout(parent, new BasicProgressMonitor());
        // Assert the ordering of n1 and ni.
        assertTrue("Node order of n1 and ni", n1.getY() < ni.getY());
        // Assert the ordering of p3 and p4. 
        assertTrue("p4 above p3", p4.getY() < p3.getY());
    }

    @Test
    public void testOutgoingNorthIncomingNorth() {
        //  ----------
        // |          |
        // |  |----|  |
        // |  v    |  v
        //p1  p2  p3  p4
        //  n1      n2

        ElkNode parent = ElkGraphUtil.createGraph();
        ElkNode n1 = ElkGraphUtil.createNode(parent);
        ElkNode n2 = ElkGraphUtil.createNode(parent);
        n1.setDimensions(20, 20);
        n2.setDimensions(20, 20);

        parent.setProperty(CoreOptions.ALGORITHM, LayeredOptions.ALGORITHM_ID);
        // Model order configuration.
        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY, OrderingStrategy.PREFER_EDGES);
        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_PORT_MODEL_ORDER, false);
        // No crossing minimization, hence only model order is used.
        parent.setProperty(LayeredOptions.CROSSING_MINIMIZATION_STRATEGY, CrossingMinimizationStrategy.NONE);
        parent.setProperty(LayeredOptions.CROSSING_MINIMIZATION_GREEDY_SWITCH_TYPE, GreedySwitchType.OFF);
        // For presentation.
        parent.setProperty(CoreOptions.DIRECTION, Direction.RIGHT);
        parent.setProperty(CoreOptions.PADDING, new ElkPadding(0.0));
        parent.setProperty(CoreOptions.SPACING_NODE_NODE, 10.0);
        parent.setProperty(LayeredOptions.SPACING_NODE_NODE_BETWEEN_LAYERS, 20.0);
        n1.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        n2.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        // Make sure that the order is correct.
        parent.setProperty(LayeredOptions.CYCLE_BREAKING_STRATEGY, CycleBreakingStrategy.MODEL_ORDER);
        
        ElkPort p1 = ElkGraphUtil.createPort(n1);
        p1.setProperty(LayeredOptions.PORT_SIDE, PortSide.NORTH);
        ElkPort p2 = ElkGraphUtil.createPort(n1);
        p2.setProperty(LayeredOptions.PORT_SIDE, PortSide.NORTH);
        ElkPort p3 = ElkGraphUtil.createPort(n2);
        p3.setProperty(LayeredOptions.PORT_SIDE, PortSide.NORTH);
        ElkPort p4 = ElkGraphUtil.createPort(n2);
        p4.setProperty(LayeredOptions.PORT_SIDE, PortSide.NORTH);

        ElkGraphUtil.createSimpleEdge(p1, p4);
        ElkGraphUtil.createSimpleEdge(p2, p3);
        
        LayeredLayoutProvider layoutProvider = new LayeredLayoutProvider();
        layoutProvider.layout(parent, new BasicProgressMonitor());
        // Assert the ordering of p1 and p2.
        assertTrue("p1 before p2", p1.getX() < p2.getX());
        // Assert the ordering of p3 and p4. 
        assertTrue("p3 before p4", p3.getX() < p4.getX());
        
        // Check for port model order, which should result in the same graph.

        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_PORT_MODEL_ORDER, true);

        layoutProvider.layout(parent, new BasicProgressMonitor());
        // Assert the ordering of p1 and p2.
        assertTrue("p1 before p2", p1.getX() < p2.getX());
        // Assert the ordering of p3 and p4. 
        assertTrue("p3 before p4", p3.getX() < p4.getX());
    }

    @Test
    public void testOutgoingNorthIncomingNorthMultipleNodes() {
        //  |--|
        //  |  |
        // p1  |
        // n1  |
        //     |-------
        //            |
        //  |------|  |
        //  |     p3  p4
        //  v       n2
        // p2
        // ni

        ElkNode parent = ElkGraphUtil.createGraph();
        ElkNode n1 = ElkGraphUtil.createNode(parent);
        ElkNode ni = ElkGraphUtil.createNode(parent);
        ElkNode n2 = ElkGraphUtil.createNode(parent);
        n1.setDimensions(20, 20);
        ni.setDimensions(20, 20);
        n2.setDimensions(20, 20);

        parent.setProperty(CoreOptions.ALGORITHM, LayeredOptions.ALGORITHM_ID);
        // Model order configuration.
        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY, OrderingStrategy.PREFER_EDGES);
        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_PORT_MODEL_ORDER, false);
        // No crossing minimization, hence only model order is used.
        parent.setProperty(LayeredOptions.CROSSING_MINIMIZATION_STRATEGY, CrossingMinimizationStrategy.NONE);
        parent.setProperty(LayeredOptions.CROSSING_MINIMIZATION_GREEDY_SWITCH_TYPE, GreedySwitchType.OFF);
        // For presentation.
        parent.setProperty(CoreOptions.DIRECTION, Direction.RIGHT);
        parent.setProperty(CoreOptions.PADDING, new ElkPadding(0.0));
        parent.setProperty(CoreOptions.SPACING_NODE_NODE, 10.0);
        parent.setProperty(LayeredOptions.SPACING_NODE_NODE_BETWEEN_LAYERS, 20.0);
        n1.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        ni.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        n2.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        // Make sure that the order is correct.
        parent.setProperty(LayeredOptions.CYCLE_BREAKING_STRATEGY, CycleBreakingStrategy.MODEL_ORDER);
        
        ElkPort p1 = ElkGraphUtil.createPort(n1);
        p1.setProperty(LayeredOptions.PORT_SIDE, PortSide.NORTH);
        ElkPort p2 = ElkGraphUtil.createPort(ni);
        p2.setProperty(LayeredOptions.PORT_SIDE, PortSide.NORTH);
        ElkPort p3 = ElkGraphUtil.createPort(n2);
        p3.setProperty(LayeredOptions.PORT_SIDE, PortSide.NORTH);
        ElkPort p4 = ElkGraphUtil.createPort(n2);
        p4.setProperty(LayeredOptions.PORT_SIDE, PortSide.NORTH);

        ElkGraphUtil.createSimpleEdge(p1, p4);
        ElkGraphUtil.createSimpleEdge(p2, p3);
        
        LayeredLayoutProvider layoutProvider = new LayeredLayoutProvider();
        layoutProvider.layout(parent, new BasicProgressMonitor());
        // Assert the ordering of n1 and ni.
        assertTrue("Node order of n1 and ni", n1.getY() < ni.getY());
        // Assert the ordering of p3 and p4. 
        assertTrue("p3 before p4", p3.getX() < p4.getX());
        
        // Check for port model order, which should result in the same graph.

        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_PORT_MODEL_ORDER, true);

        layoutProvider.layout(parent, new BasicProgressMonitor());
        // Assert the ordering of n1 and ni.
        assertTrue("Node order of n1 and ni", n1.getY() < ni.getY());
        // Assert the ordering of p3 and p4. 
        assertTrue("p3 before p4", p3.getX() < p4.getX());
    }

    @Test
    public void testOutgoingNorthIncomingSouth() {
        //      n2
        //    p4  p3
        // ____A  |
        // |  ____|     
        // |  |
        // |  v
        //p1  p2
        //  n1

        ElkNode parent = ElkGraphUtil.createGraph();
        ElkNode n1 = ElkGraphUtil.createNode(parent);
        ElkNode n2 = ElkGraphUtil.createNode(parent);
        n1.setDimensions(20, 20);
        n2.setDimensions(20, 20);

        parent.setProperty(CoreOptions.ALGORITHM, LayeredOptions.ALGORITHM_ID);
        // Model order configuration.
        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY, OrderingStrategy.PREFER_EDGES);
        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_PORT_MODEL_ORDER, false);
        // No crossing minimization, hence only model order is used.
        parent.setProperty(LayeredOptions.CROSSING_MINIMIZATION_STRATEGY, CrossingMinimizationStrategy.NONE);
        parent.setProperty(LayeredOptions.CROSSING_MINIMIZATION_GREEDY_SWITCH_TYPE, GreedySwitchType.OFF);
        // For presentation.
        parent.setProperty(CoreOptions.DIRECTION, Direction.RIGHT);
        parent.setProperty(CoreOptions.PADDING, new ElkPadding(0.0));
        parent.setProperty(CoreOptions.SPACING_NODE_NODE, 10.0);
        parent.setProperty(LayeredOptions.SPACING_NODE_NODE_BETWEEN_LAYERS, 20.0);
        n1.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        n2.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        // Make sure that the order is correct.
        parent.setProperty(LayeredOptions.CYCLE_BREAKING_STRATEGY, CycleBreakingStrategy.MODEL_ORDER);
        
        ElkPort p1 = ElkGraphUtil.createPort(n1);
        p1.setProperty(LayeredOptions.PORT_SIDE, PortSide.NORTH);
        ElkPort p2 = ElkGraphUtil.createPort(n1);
        p2.setProperty(LayeredOptions.PORT_SIDE, PortSide.NORTH);
        ElkPort p3 = ElkGraphUtil.createPort(n2);
        p3.setProperty(LayeredOptions.PORT_SIDE, PortSide.SOUTH);
        ElkPort p4 = ElkGraphUtil.createPort(n2);
        p4.setProperty(LayeredOptions.PORT_SIDE, PortSide.SOUTH);

        ElkGraphUtil.createSimpleEdge(p1, p4);
        ElkGraphUtil.createSimpleEdge(p2, p3);
        
        LayeredLayoutProvider layoutProvider = new LayeredLayoutProvider();
        layoutProvider.layout(parent, new BasicProgressMonitor());
        // Assert the ordering of p1 and p2.
        assertTrue("p1 before p2", p1.getX() < p2.getX());
        // Assert the ordering of p3 and p4. 
        assertTrue("p4 before p3", p4.getX() < p3.getX());
        
        // Check for port model order, which should result in the same graph.

        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_PORT_MODEL_ORDER, true);

        layoutProvider.layout(parent, new BasicProgressMonitor());
        // Assert the ordering of p1 and p2.
        assertTrue("p1 before p2", p1.getX() < p2.getX());
        // Assert the ordering of p3 and p4. 
        assertTrue("p4 before p3", p4.getX() < p3.getX());
    }

    @Test
    public void testOutgoingNorthIncomingSouthMultipleNodes() {
        //          n2
        //        p4  p3
        //  |------|  |
        // p1         |
        // n1         |
        //  |---------|
        //  v
        // p2
        // ni

        ElkNode parent = ElkGraphUtil.createGraph();
        ElkNode n1 = ElkGraphUtil.createNode(parent);
        ElkNode ni = ElkGraphUtil.createNode(parent);
        ElkNode n2 = ElkGraphUtil.createNode(parent);
        n1.setDimensions(20, 20);
        ni.setDimensions(20, 20);
        n2.setDimensions(20, 20);

        parent.setProperty(CoreOptions.ALGORITHM, LayeredOptions.ALGORITHM_ID);
        // Model order configuration.
        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY, OrderingStrategy.PREFER_EDGES);
        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_PORT_MODEL_ORDER, false);
        // No crossing minimization, hence only model order is used.
        parent.setProperty(LayeredOptions.CROSSING_MINIMIZATION_STRATEGY, CrossingMinimizationStrategy.NONE);
        parent.setProperty(LayeredOptions.CROSSING_MINIMIZATION_GREEDY_SWITCH_TYPE, GreedySwitchType.OFF);
        // For presentation.
        parent.setProperty(CoreOptions.DIRECTION, Direction.RIGHT);
        parent.setProperty(CoreOptions.PADDING, new ElkPadding(0.0));
        parent.setProperty(CoreOptions.SPACING_NODE_NODE, 10.0);
        parent.setProperty(LayeredOptions.SPACING_NODE_NODE_BETWEEN_LAYERS, 20.0);
        n1.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        ni.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        n2.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        // Make sure that the order is correct.
        parent.setProperty(LayeredOptions.CYCLE_BREAKING_STRATEGY, CycleBreakingStrategy.MODEL_ORDER);
        
        ElkPort p1 = ElkGraphUtil.createPort(n1);
        p1.setProperty(LayeredOptions.PORT_SIDE, PortSide.NORTH);
        ElkPort p2 = ElkGraphUtil.createPort(ni);
        p2.setProperty(LayeredOptions.PORT_SIDE, PortSide.NORTH);
        ElkPort p3 = ElkGraphUtil.createPort(n2);
        p3.setProperty(LayeredOptions.PORT_SIDE, PortSide.SOUTH);
        ElkPort p4 = ElkGraphUtil.createPort(n2);
        p4.setProperty(LayeredOptions.PORT_SIDE, PortSide.SOUTH);

        ElkGraphUtil.createSimpleEdge(p1, p4);
        ElkGraphUtil.createSimpleEdge(p2, p3);
        
        LayeredLayoutProvider layoutProvider = new LayeredLayoutProvider();
        layoutProvider.layout(parent, new BasicProgressMonitor());
        // Assert the ordering of n1 and ni.
        assertTrue("Node order of n1 and ni", n1.getY() < ni.getY());
        // Assert the ordering of p3 and p4. 
        assertTrue("p4 before p3", p4.getX() < p3.getX());
        
        // Check for port model order, which should result in the same graph.

        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_PORT_MODEL_ORDER, true);

        layoutProvider.layout(parent, new BasicProgressMonitor());
        // Assert the ordering of n1 and ni.
        assertTrue("Node order of n1 and ni", n1.getY() < ni.getY());
        // Assert the ordering of p3 and p4. 
        assertTrue("p4 before p3", p4.getX() < p3.getX());
    }

    @Test
    public void testOutgoingNorthIncomingEast() {
        //    p3----
        //  n2      |
        //    p4<-  |
        //        | |
        //  ------  |
        // |   -----     
        // |  |
        // |  v
        //p1  p2
        //  n1

        ElkNode parent = ElkGraphUtil.createGraph();
        ElkNode n1 = ElkGraphUtil.createNode(parent);
        ElkNode n2 = ElkGraphUtil.createNode(parent);
        n1.setDimensions(20, 20);
        n2.setDimensions(20, 20);

        parent.setProperty(CoreOptions.ALGORITHM, LayeredOptions.ALGORITHM_ID);
        // Model order configuration.
        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY, OrderingStrategy.PREFER_EDGES);
        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_PORT_MODEL_ORDER, false);
        // No crossing minimization, hence only model order is used.
        parent.setProperty(LayeredOptions.CROSSING_MINIMIZATION_STRATEGY, CrossingMinimizationStrategy.NONE);
        parent.setProperty(LayeredOptions.CROSSING_MINIMIZATION_GREEDY_SWITCH_TYPE, GreedySwitchType.OFF);
        // For presentation.
        parent.setProperty(CoreOptions.DIRECTION, Direction.RIGHT);
        parent.setProperty(CoreOptions.PADDING, new ElkPadding(0.0));
        parent.setProperty(CoreOptions.SPACING_NODE_NODE, 10.0);
        parent.setProperty(LayeredOptions.SPACING_NODE_NODE_BETWEEN_LAYERS, 20.0);
        n1.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        n2.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        // Make sure that the order is correct.
        parent.setProperty(LayeredOptions.CYCLE_BREAKING_STRATEGY, CycleBreakingStrategy.MODEL_ORDER);
        
        ElkPort p1 = ElkGraphUtil.createPort(n1);
        p1.setProperty(LayeredOptions.PORT_SIDE, PortSide.NORTH);
        ElkPort p2 = ElkGraphUtil.createPort(n1);
        p2.setProperty(LayeredOptions.PORT_SIDE, PortSide.NORTH);
        ElkPort p3 = ElkGraphUtil.createPort(n2);
        p3.setProperty(LayeredOptions.PORT_SIDE, PortSide.EAST);
        ElkPort p4 = ElkGraphUtil.createPort(n2);
        p4.setProperty(LayeredOptions.PORT_SIDE, PortSide.EAST);

        ElkGraphUtil.createSimpleEdge(p1, p4);
        ElkGraphUtil.createSimpleEdge(p2, p3);
        
        LayeredLayoutProvider layoutProvider = new LayeredLayoutProvider();
        layoutProvider.layout(parent, new BasicProgressMonitor());
        // Assert the ordering of p1 and p2.
        assertTrue("p1 before p2", p1.getX() < p2.getX());
        // Assert the ordering of p3 and p4. 
        assertTrue("p3 above p4", p3.getY() < p4.getY());
        
        // Check for port model order, which should result in the same graph.

        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_PORT_MODEL_ORDER, true);

        layoutProvider.layout(parent, new BasicProgressMonitor());
        // Assert the ordering of p1 and p2.
        assertTrue("p1 before p2", p1.getX() < p2.getX());
        // Assert the ordering of p3 and p4. 
        assertTrue("p3 above p4", p3.getY() < p4.getY());
    }

    @Test
    public void testOutgoingNorthIncomingEastMultipleNodes() {
        //      p3----|
        //    n1      |
        //      p4<|  |
        //  |------|  |
        // p1         |
        // n1         |
        //  |---------|
        //  v
        // p2
        // ni

        ElkNode parent = ElkGraphUtil.createGraph();
        ElkNode n1 = ElkGraphUtil.createNode(parent);
        ElkNode ni = ElkGraphUtil.createNode(parent);
        ElkNode n2 = ElkGraphUtil.createNode(parent);
        n1.setDimensions(20, 20);
        ni.setDimensions(20, 20);
        n2.setDimensions(20, 20);

        parent.setProperty(CoreOptions.ALGORITHM, LayeredOptions.ALGORITHM_ID);
        // Model order configuration.
        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY, OrderingStrategy.PREFER_EDGES);
        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_PORT_MODEL_ORDER, false);
        // No crossing minimization, hence only model order is used.
        parent.setProperty(LayeredOptions.CROSSING_MINIMIZATION_STRATEGY, CrossingMinimizationStrategy.NONE);
        parent.setProperty(LayeredOptions.CROSSING_MINIMIZATION_GREEDY_SWITCH_TYPE, GreedySwitchType.OFF);
        // For presentation.
        parent.setProperty(CoreOptions.DIRECTION, Direction.RIGHT);
        parent.setProperty(CoreOptions.PADDING, new ElkPadding(0.0));
        parent.setProperty(CoreOptions.SPACING_NODE_NODE, 10.0);
        parent.setProperty(LayeredOptions.SPACING_NODE_NODE_BETWEEN_LAYERS, 20.0);
        n1.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        ni.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        n2.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        // Make sure that the order is correct.
        parent.setProperty(LayeredOptions.CYCLE_BREAKING_STRATEGY, CycleBreakingStrategy.MODEL_ORDER);
        
        ElkPort p1 = ElkGraphUtil.createPort(n1);
        p1.setProperty(LayeredOptions.PORT_SIDE, PortSide.NORTH);
        ElkPort p2 = ElkGraphUtil.createPort(ni);
        p2.setProperty(LayeredOptions.PORT_SIDE, PortSide.NORTH);
        ElkPort p3 = ElkGraphUtil.createPort(n2);
        p3.setProperty(LayeredOptions.PORT_SIDE, PortSide.EAST);
        ElkPort p4 = ElkGraphUtil.createPort(n2);
        p4.setProperty(LayeredOptions.PORT_SIDE, PortSide.EAST);

        ElkGraphUtil.createSimpleEdge(p1, p4);
        ElkGraphUtil.createSimpleEdge(p2, p3);
        
        LayeredLayoutProvider layoutProvider = new LayeredLayoutProvider();
        layoutProvider.layout(parent, new BasicProgressMonitor());
        // Assert the ordering of n1 and ni.
        assertTrue("Node order of n1 and ni", n1.getY() < ni.getY());
        // Assert the ordering of p3 and p4. 
        assertTrue("p3 above p4", p3.getY() < p4.getY());
        
        // Check for port model order, which should result in the same graph.

        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_PORT_MODEL_ORDER, true);

        layoutProvider.layout(parent, new BasicProgressMonitor());
        // Assert the ordering of n1 and ni.
        assertTrue("Node order of n1 and ni", n1.getY() < ni.getY());
        // Assert the ordering of p3 and p4. 
        assertTrue("p3 above p4", p3.getY() < p4.getY());
    }

    @Test
    public void testOutgoingSouthIncomingWest() {
        //  n1
        //p1  p2
        // |  A
        // |  |---p3
        // |        n2
        // |----->p4

        ElkNode parent = ElkGraphUtil.createGraph();
        ElkNode n1 = ElkGraphUtil.createNode(parent);
        ElkNode n2 = ElkGraphUtil.createNode(parent);
        n1.setDimensions(20, 20);
        n2.setDimensions(20, 20);

        parent.setProperty(CoreOptions.ALGORITHM, LayeredOptions.ALGORITHM_ID);
        // Model order configuration.
        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY, OrderingStrategy.PREFER_EDGES);
        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_PORT_MODEL_ORDER, false);
        // No crossing minimization, hence only model order is used.
        parent.setProperty(LayeredOptions.CROSSING_MINIMIZATION_STRATEGY, CrossingMinimizationStrategy.NONE);
        parent.setProperty(LayeredOptions.CROSSING_MINIMIZATION_GREEDY_SWITCH_TYPE, GreedySwitchType.OFF);
        // For presentation.
        parent.setProperty(CoreOptions.DIRECTION, Direction.RIGHT);
        parent.setProperty(CoreOptions.PADDING, new ElkPadding(0.0));
        parent.setProperty(CoreOptions.SPACING_NODE_NODE, 10.0);
        parent.setProperty(LayeredOptions.SPACING_NODE_NODE_BETWEEN_LAYERS, 20.0);
        n1.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        n2.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        // Make sure that the order is correct.
        parent.setProperty(LayeredOptions.CYCLE_BREAKING_STRATEGY, CycleBreakingStrategy.MODEL_ORDER);
        
        ElkPort p1 = ElkGraphUtil.createPort(n1);
        p1.setProperty(LayeredOptions.PORT_SIDE, PortSide.SOUTH);
        ElkPort p2 = ElkGraphUtil.createPort(n1);
        p2.setProperty(LayeredOptions.PORT_SIDE, PortSide.SOUTH);
        ElkPort p3 = ElkGraphUtil.createPort(n2);
        p3.setProperty(LayeredOptions.PORT_SIDE, PortSide.WEST);
        ElkPort p4 = ElkGraphUtil.createPort(n2);
        p4.setProperty(LayeredOptions.PORT_SIDE, PortSide.WEST);

        ElkGraphUtil.createSimpleEdge(p1, p4);
        ElkGraphUtil.createSimpleEdge(p2, p3);
        
        LayeredLayoutProvider layoutProvider = new LayeredLayoutProvider();
        layoutProvider.layout(parent, new BasicProgressMonitor());
        // Assert the ordering of p1 and p2.
        assertTrue("p1 before p2", p1.getX() < p2.getX());
        // Assert the ordering of p3 and p4. 
        assertTrue("p3 above p4", p3.getY() < p4.getY());
        
        // Check for port model order, which should result in the same graph.

        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_PORT_MODEL_ORDER, true);

        layoutProvider.layout(parent, new BasicProgressMonitor());
        // Assert the ordering of p1 and p2.
        assertTrue("p1 before p2", p1.getX() < p2.getX());
        // Assert the ordering of p3 and p4. 
        assertTrue("p3 above p4", p3.getY() < p4.getY());
    }

    @Test
    public void testOutgoingSouthIncomingWestMultipleNodes() {
        // n1
        // p1
        // |----->p4
        // ni       n2
        // p2  |--p3
        // A---|
        //
        //
        //

        ElkNode parent = ElkGraphUtil.createGraph();
        ElkNode n1 = ElkGraphUtil.createNode(parent);
        ElkNode ni = ElkGraphUtil.createNode(parent);
        ElkNode n2 = ElkGraphUtil.createNode(parent);
        n1.setDimensions(20, 20);
        ni.setDimensions(20, 20);
        n2.setDimensions(20, 20);

        parent.setProperty(CoreOptions.ALGORITHM, LayeredOptions.ALGORITHM_ID);
        // Model order configuration.
        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY, OrderingStrategy.PREFER_EDGES);
        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_PORT_MODEL_ORDER, false);
        // No crossing minimization, hence only model order is used.
        parent.setProperty(LayeredOptions.CROSSING_MINIMIZATION_STRATEGY, CrossingMinimizationStrategy.NONE);
        parent.setProperty(LayeredOptions.CROSSING_MINIMIZATION_GREEDY_SWITCH_TYPE, GreedySwitchType.OFF);
        // For presentation.
        parent.setProperty(CoreOptions.DIRECTION, Direction.RIGHT);
        parent.setProperty(CoreOptions.PADDING, new ElkPadding(0.0));
        parent.setProperty(CoreOptions.SPACING_NODE_NODE, 10.0);
        parent.setProperty(LayeredOptions.SPACING_NODE_NODE_BETWEEN_LAYERS, 20.0);
        n1.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        ni.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        n2.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        // Make sure that the order is correct.
        parent.setProperty(LayeredOptions.CYCLE_BREAKING_STRATEGY, CycleBreakingStrategy.MODEL_ORDER);
        
        ElkPort p1 = ElkGraphUtil.createPort(n1);
        p1.setProperty(LayeredOptions.PORT_SIDE, PortSide.SOUTH);
        ElkPort p2 = ElkGraphUtil.createPort(ni);
        p2.setProperty(LayeredOptions.PORT_SIDE, PortSide.SOUTH);
        ElkPort p3 = ElkGraphUtil.createPort(n2);
        p3.setProperty(LayeredOptions.PORT_SIDE, PortSide.WEST);
        ElkPort p4 = ElkGraphUtil.createPort(n2);
        p4.setProperty(LayeredOptions.PORT_SIDE, PortSide.WEST);

        ElkGraphUtil.createSimpleEdge(p1, p4);
        ElkGraphUtil.createSimpleEdge(p2, p3);
        
        LayeredLayoutProvider layoutProvider = new LayeredLayoutProvider();
        layoutProvider.layout(parent, new BasicProgressMonitor());
        // Assert the ordering of n1 and ni.
        assertTrue("Node order of n1 and ni", n1.getY() < ni.getY());
        // Assert the ordering of p3 and p4. 
        assertTrue("p4 above p3", p4.getY() < p3.getY());
        
        // Check for port model order, which should result in the same graph.

        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_PORT_MODEL_ORDER, true);

        layoutProvider.layout(parent, new BasicProgressMonitor());
        // Assert the ordering of n1 and ni.
        assertTrue("Node order of n1 and ni", n1.getY() < ni.getY());
        // Assert the ordering of p3 and p4. 
        assertTrue("p4 above p3", p4.getY() < p3.getY());
    }

    @Test
    public void testOutgoingSouthIncomingNorth() {
        //  n1
        //p1  p2
        // |  A
        // |  |---|
        // |      |
        // |---|  |
        //     v  |
        //    p4  p3
        //      n2

        ElkNode parent = ElkGraphUtil.createGraph();
        ElkNode n1 = ElkGraphUtil.createNode(parent);
        ElkNode n2 = ElkGraphUtil.createNode(parent);
        n1.setDimensions(20, 20);
        n2.setDimensions(20, 20);

        parent.setProperty(CoreOptions.ALGORITHM, LayeredOptions.ALGORITHM_ID);
        // Model order configuration.
        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY, OrderingStrategy.PREFER_EDGES);
        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_PORT_MODEL_ORDER, false);
        // No crossing minimization, hence only model order is used.
        parent.setProperty(LayeredOptions.CROSSING_MINIMIZATION_STRATEGY, CrossingMinimizationStrategy.NONE);
        parent.setProperty(LayeredOptions.CROSSING_MINIMIZATION_GREEDY_SWITCH_TYPE, GreedySwitchType.OFF);
        // For presentation.
        parent.setProperty(CoreOptions.DIRECTION, Direction.RIGHT);
        parent.setProperty(CoreOptions.PADDING, new ElkPadding(0.0));
        parent.setProperty(CoreOptions.SPACING_NODE_NODE, 10.0);
        parent.setProperty(LayeredOptions.SPACING_NODE_NODE_BETWEEN_LAYERS, 20.0);
        n1.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        n2.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        // Make sure that the order is correct.
        parent.setProperty(LayeredOptions.CYCLE_BREAKING_STRATEGY, CycleBreakingStrategy.MODEL_ORDER);
        
        ElkPort p1 = ElkGraphUtil.createPort(n1);
        p1.setProperty(LayeredOptions.PORT_SIDE, PortSide.SOUTH);
        ElkPort p2 = ElkGraphUtil.createPort(n1);
        p2.setProperty(LayeredOptions.PORT_SIDE, PortSide.SOUTH);
        ElkPort p3 = ElkGraphUtil.createPort(n2);
        p3.setProperty(LayeredOptions.PORT_SIDE, PortSide.NORTH);
        ElkPort p4 = ElkGraphUtil.createPort(n2);
        p4.setProperty(LayeredOptions.PORT_SIDE, PortSide.NORTH);

        ElkGraphUtil.createSimpleEdge(p1, p4);
        ElkGraphUtil.createSimpleEdge(p2, p3);
        
        LayeredLayoutProvider layoutProvider = new LayeredLayoutProvider();
        layoutProvider.layout(parent, new NullElkProgressMonitor());
        // Assert the ordering of p1 and p2.
        assertTrue("p1 before p2", p1.getX() < p2.getX());
        // Assert the ordering of p3 and p4. 
        assertTrue("p4 before p3", p4.getX() < p3.getX());
        
        // Check for port model order, which should result in the same graph.

        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_PORT_MODEL_ORDER, true);

        layoutProvider.layout(parent, new BasicProgressMonitor());
        // Assert the ordering of p1 and p2.
        assertTrue("p1 before p2", p1.getX() < p2.getX());
        // Assert the ordering of p3 and p4. 
        assertTrue("p4 before p3", p4.getX() < p3.getX());
    }

    @Test
    public void testOutgoingSouthIncomingNorthMultipleNodes() {
        // n1
        // p1
        // |-----------
        //             |
        // ni          |
        // p2          |
        // A--------|  |
        //         p3  p4
        //           n2
        //

        ElkNode parent = ElkGraphUtil.createGraph();
        ElkNode n1 = ElkGraphUtil.createNode(parent);
        ElkNode ni = ElkGraphUtil.createNode(parent);
        ElkNode n2 = ElkGraphUtil.createNode(parent);
        n1.setDimensions(20, 20);
        ni.setDimensions(20, 20);
        n2.setDimensions(20, 20);

        parent.setProperty(CoreOptions.ALGORITHM, LayeredOptions.ALGORITHM_ID);
        // Model order configuration.
        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY, OrderingStrategy.PREFER_EDGES);
        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_PORT_MODEL_ORDER, false);
        // No crossing minimization, hence only model order is used.
        parent.setProperty(LayeredOptions.CROSSING_MINIMIZATION_STRATEGY, CrossingMinimizationStrategy.NONE);
        parent.setProperty(LayeredOptions.CROSSING_MINIMIZATION_GREEDY_SWITCH_TYPE, GreedySwitchType.OFF);
        // For presentation.
        parent.setProperty(CoreOptions.DIRECTION, Direction.RIGHT);
        parent.setProperty(CoreOptions.PADDING, new ElkPadding(0.0));
        parent.setProperty(CoreOptions.SPACING_NODE_NODE, 10.0);
        parent.setProperty(LayeredOptions.SPACING_NODE_NODE_BETWEEN_LAYERS, 20.0);
        n1.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        ni.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        n2.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        // Make sure that the order is correct.
        parent.setProperty(LayeredOptions.CYCLE_BREAKING_STRATEGY, CycleBreakingStrategy.MODEL_ORDER);
        
        ElkPort p1 = ElkGraphUtil.createPort(n1);
        p1.setProperty(LayeredOptions.PORT_SIDE, PortSide.SOUTH);
        ElkPort p2 = ElkGraphUtil.createPort(ni);
        p2.setProperty(LayeredOptions.PORT_SIDE, PortSide.SOUTH);
        ElkPort p3 = ElkGraphUtil.createPort(n2);
        p3.setProperty(LayeredOptions.PORT_SIDE, PortSide.NORTH);
        ElkPort p4 = ElkGraphUtil.createPort(n2);
        p4.setProperty(LayeredOptions.PORT_SIDE, PortSide.NORTH);

        ElkGraphUtil.createSimpleEdge(p1, p4);
        ElkGraphUtil.createSimpleEdge(p2, p3);
        
        LayeredLayoutProvider layoutProvider = new LayeredLayoutProvider();
        layoutProvider.layout(parent, new NullElkProgressMonitor());
        // Assert the ordering of n1 and ni.
        assertTrue("Node order of n1 and ni", n1.getY() < ni.getY());
        // Assert the ordering of p3 and p4. 
        assertTrue("p3 before p4", p3.getX() < p4.getX());
        
        // Check for port model order, which should result in the same graph.

        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_PORT_MODEL_ORDER, true);

        layoutProvider.layout(parent, new BasicProgressMonitor());
        // Assert the ordering of n1 and ni.
        assertTrue("Node order of n1 and ni", n1.getY() < ni.getY());
        // Assert the ordering of p3 and p4. 
        assertTrue("p3 before p4", p3.getX() < p4.getX());
    }

    @Test
    public void testOutgoingSouthIncomingSouth() {
        //  n1     n2
        //p1  p2 p3  p4
        // |  A   |  A
        // |  |---|  |
        // |         |
        // |---------|

        ElkNode parent = ElkGraphUtil.createGraph();
        ElkNode n1 = ElkGraphUtil.createNode(parent);
        ElkNode n2 = ElkGraphUtil.createNode(parent);
        n1.setDimensions(20, 20);
        n2.setDimensions(20, 20);

        parent.setProperty(CoreOptions.ALGORITHM, LayeredOptions.ALGORITHM_ID);
        // Model order configuration.
        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY, OrderingStrategy.PREFER_EDGES);
        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_PORT_MODEL_ORDER, false);
        // No crossing minimization, hence only model order is used.
        parent.setProperty(LayeredOptions.CROSSING_MINIMIZATION_STRATEGY, CrossingMinimizationStrategy.NONE);
        parent.setProperty(LayeredOptions.CROSSING_MINIMIZATION_GREEDY_SWITCH_TYPE, GreedySwitchType.OFF);
        // For presentation.
        parent.setProperty(CoreOptions.DIRECTION, Direction.RIGHT);
        parent.setProperty(CoreOptions.PADDING, new ElkPadding(0.0));
        parent.setProperty(CoreOptions.SPACING_NODE_NODE, 10.0);
        parent.setProperty(LayeredOptions.SPACING_NODE_NODE_BETWEEN_LAYERS, 20.0);
        n1.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        n2.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        // Make sure that the order is correct.
        parent.setProperty(LayeredOptions.CYCLE_BREAKING_STRATEGY, CycleBreakingStrategy.MODEL_ORDER);
        
        ElkPort p1 = ElkGraphUtil.createPort(n1);
        p1.setProperty(LayeredOptions.PORT_SIDE, PortSide.SOUTH);
        ElkPort p2 = ElkGraphUtil.createPort(n1);
        p2.setProperty(LayeredOptions.PORT_SIDE, PortSide.SOUTH);
        ElkPort p3 = ElkGraphUtil.createPort(n2);
        p3.setProperty(LayeredOptions.PORT_SIDE, PortSide.SOUTH);
        ElkPort p4 = ElkGraphUtil.createPort(n2);
        p4.setProperty(LayeredOptions.PORT_SIDE, PortSide.SOUTH);

        ElkGraphUtil.createSimpleEdge(p1, p4);
        ElkGraphUtil.createSimpleEdge(p2, p3);
        
        LayeredLayoutProvider layoutProvider = new LayeredLayoutProvider();
        layoutProvider.layout(parent, new BasicProgressMonitor());
        // Assert the ordering of p1 and p2.
        assertTrue("p1 before p2", p1.getX() < p2.getX());
        // Assert the ordering of p3 and p4. 
        assertTrue("p3 before p4", p3.getX() < p4.getX());
        
        // Check for port model order, which should result in the same graph.

        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_PORT_MODEL_ORDER, true);

        layoutProvider.layout(parent, new BasicProgressMonitor());
        // Assert the ordering of p1 and p2.
        assertTrue("p1 before p2", p1.getX() < p2.getX());
        // Assert the ordering of p3 and p4. 
        assertTrue("p3 before p4", p3.getX() < p4.getX());
    }

    @Test
    public void testOutgoingSouthIncomingSouthMultipleNodes() {
        //n1       n2
        //p1     p4  p3
        // |      |  A
        // |------|  |
        //           |
        //ni         |
        //p2         |
        //A          |
        //-----------|
        
        ElkNode parent = ElkGraphUtil.createGraph();
        ElkNode n1 = ElkGraphUtil.createNode(parent);
        ElkNode ni = ElkGraphUtil.createNode(parent);
        ElkNode n2 = ElkGraphUtil.createNode(parent);
        n1.setDimensions(20, 20);
        ni.setDimensions(20, 20);
        n2.setDimensions(20, 20);

        parent.setProperty(CoreOptions.ALGORITHM, LayeredOptions.ALGORITHM_ID);
        // Model order configuration.
        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY, OrderingStrategy.PREFER_EDGES);
        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_PORT_MODEL_ORDER, false);
        // No crossing minimization, hence only model order is used.
        parent.setProperty(LayeredOptions.CROSSING_MINIMIZATION_STRATEGY, CrossingMinimizationStrategy.NONE);
        parent.setProperty(LayeredOptions.CROSSING_MINIMIZATION_GREEDY_SWITCH_TYPE, GreedySwitchType.OFF);
        // For presentation.
        parent.setProperty(CoreOptions.DIRECTION, Direction.RIGHT);
        parent.setProperty(CoreOptions.PADDING, new ElkPadding(0.0));
        parent.setProperty(CoreOptions.SPACING_NODE_NODE, 10.0);
        parent.setProperty(LayeredOptions.SPACING_NODE_NODE_BETWEEN_LAYERS, 20.0);
        n1.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        ni.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        n2.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        // Make sure that the order is correct.
        parent.setProperty(LayeredOptions.CYCLE_BREAKING_STRATEGY, CycleBreakingStrategy.MODEL_ORDER);
        
        ElkPort p1 = ElkGraphUtil.createPort(n1);
        p1.setProperty(LayeredOptions.PORT_SIDE, PortSide.SOUTH);
        ElkPort p2 = ElkGraphUtil.createPort(ni);
        p2.setProperty(LayeredOptions.PORT_SIDE, PortSide.SOUTH);
        ElkPort p3 = ElkGraphUtil.createPort(n2);
        p3.setProperty(LayeredOptions.PORT_SIDE, PortSide.SOUTH);
        ElkPort p4 = ElkGraphUtil.createPort(n2);
        p4.setProperty(LayeredOptions.PORT_SIDE, PortSide.SOUTH);

        ElkGraphUtil.createSimpleEdge(p1, p4);
        ElkGraphUtil.createSimpleEdge(p2, p3);
        
        LayeredLayoutProvider layoutProvider = new LayeredLayoutProvider();
        layoutProvider.layout(parent, new BasicProgressMonitor());
        // Assert the ordering of n1 and ni.
        assertTrue("Node order of n1 and ni", n1.getY() < ni.getY());
        // Assert the ordering of p3 and p4. 
        assertTrue("p4 before p3", p4.getX() < p3.getX());
        
        // Check for port model order, which should result in the same graph.

        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_PORT_MODEL_ORDER, true);

        layoutProvider.layout(parent, new BasicProgressMonitor());
        // Assert the ordering of n1 and ni.
        assertTrue("Node order of n1 and ni", n1.getY() < ni.getY());
        // Assert the ordering of p3 and p4. 
        assertTrue("p4 before p3", p4.getX() < p3.getX());
    }

    @Test
    public void testOutgoingSouthIncomingEast() {
        //
        //          p4<-|
        //  n1    n2    |
        //p1  p2    p3  |
        // |  A      |  |
        // |  |------|  |  
        // |            |
        // |------------|

        ElkNode parent = ElkGraphUtil.createGraph();
        ElkNode n1 = ElkGraphUtil.createNode(parent);
        ElkNode n2 = ElkGraphUtil.createNode(parent);
        n1.setDimensions(20, 20);
        n2.setDimensions(20, 20);

        parent.setProperty(CoreOptions.ALGORITHM, LayeredOptions.ALGORITHM_ID);
        // Model order configuration.
        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY, OrderingStrategy.PREFER_EDGES);
        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_PORT_MODEL_ORDER, false);
        // No crossing minimization, hence only model order is used.
        parent.setProperty(LayeredOptions.CROSSING_MINIMIZATION_STRATEGY, CrossingMinimizationStrategy.NONE);
        parent.setProperty(LayeredOptions.CROSSING_MINIMIZATION_GREEDY_SWITCH_TYPE, GreedySwitchType.OFF);
        // For presentation.
        parent.setProperty(CoreOptions.DIRECTION, Direction.RIGHT);
        parent.setProperty(CoreOptions.PADDING, new ElkPadding(0.0));
        parent.setProperty(CoreOptions.SPACING_NODE_NODE, 10.0);
        parent.setProperty(LayeredOptions.SPACING_NODE_NODE_BETWEEN_LAYERS, 20.0);
        n1.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        n2.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        // Make sure that the order is correct.
        parent.setProperty(LayeredOptions.CYCLE_BREAKING_STRATEGY, CycleBreakingStrategy.MODEL_ORDER);
        
        ElkPort p1 = ElkGraphUtil.createPort(n1);
        p1.setProperty(LayeredOptions.PORT_SIDE, PortSide.SOUTH);
        ElkPort p2 = ElkGraphUtil.createPort(n1);
        p2.setProperty(LayeredOptions.PORT_SIDE, PortSide.SOUTH);
        ElkPort p3 = ElkGraphUtil.createPort(n2);
        p3.setProperty(LayeredOptions.PORT_SIDE, PortSide.EAST);
        ElkPort p4 = ElkGraphUtil.createPort(n2);
        p4.setProperty(LayeredOptions.PORT_SIDE, PortSide.EAST);

        ElkGraphUtil.createSimpleEdge(p1, p4);
        ElkGraphUtil.createSimpleEdge(p2, p3);
        
        LayeredLayoutProvider layoutProvider = new LayeredLayoutProvider();
        layoutProvider.layout(parent, new BasicProgressMonitor());
        // Assert the ordering of p1 and p2.
        assertTrue("p1 before p2", p1.getX() < p2.getX());
        // Assert the ordering of p3 and p4. 
        assertTrue("p4 above p3", p4.getY() < p3.getY());
        
        // Check for port model order, which should result in the same graph.

        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_PORT_MODEL_ORDER, true);

        layoutProvider.layout(parent, new BasicProgressMonitor());
        // Assert the ordering of p1 and p2.
        assertTrue("p1 before p2", p1.getX() < p2.getX());
        // Assert the ordering of p3 and p4. 
        assertTrue("p4 above p3", p4.getY() < p3.getY());
    }

    @Test
    public void testOutgoingSouthIncomingEastMultipleNodes() {
        //       p3--|
        //n1   n2    |
        //p1     p4< |
        // |       | |
        // |-------- |
        //           |
        //ni         |
        //p2         |
        //A          |
        //-----------|

        ElkNode parent = ElkGraphUtil.createGraph();
        ElkNode n1 = ElkGraphUtil.createNode(parent);
        ElkNode ni = ElkGraphUtil.createNode(parent);
        ElkNode n2 = ElkGraphUtil.createNode(parent);
        n1.setDimensions(20, 20);
        ni.setDimensions(20, 20);
        n2.setDimensions(20, 20);

        parent.setProperty(CoreOptions.ALGORITHM, LayeredOptions.ALGORITHM_ID);
        // Model order configuration.
        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY, OrderingStrategy.PREFER_EDGES);
        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_PORT_MODEL_ORDER, false);
        // No crossing minimization, hence only model order is used.
        parent.setProperty(LayeredOptions.CROSSING_MINIMIZATION_STRATEGY, CrossingMinimizationStrategy.NONE);
        parent.setProperty(LayeredOptions.CROSSING_MINIMIZATION_GREEDY_SWITCH_TYPE, GreedySwitchType.OFF);
        // For presentation.
        parent.setProperty(CoreOptions.DIRECTION, Direction.RIGHT);
        parent.setProperty(CoreOptions.PADDING, new ElkPadding(0.0));
        parent.setProperty(CoreOptions.SPACING_NODE_NODE, 10.0);
        parent.setProperty(LayeredOptions.SPACING_NODE_NODE_BETWEEN_LAYERS, 20.0);
        n1.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        ni.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        n2.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        // Make sure that the order is correct.
        parent.setProperty(LayeredOptions.CYCLE_BREAKING_STRATEGY, CycleBreakingStrategy.MODEL_ORDER);
        
        ElkPort p1 = ElkGraphUtil.createPort(n1);
        p1.setProperty(LayeredOptions.PORT_SIDE, PortSide.SOUTH);
        ElkPort p2 = ElkGraphUtil.createPort(ni);
        p2.setProperty(LayeredOptions.PORT_SIDE, PortSide.SOUTH);
        ElkPort p3 = ElkGraphUtil.createPort(n2);
        p3.setProperty(LayeredOptions.PORT_SIDE, PortSide.EAST);
        ElkPort p4 = ElkGraphUtil.createPort(n2);
        p4.setProperty(LayeredOptions.PORT_SIDE, PortSide.EAST);

        ElkGraphUtil.createSimpleEdge(p1, p4);
        ElkGraphUtil.createSimpleEdge(p2, p3);
        
        LayeredLayoutProvider layoutProvider = new LayeredLayoutProvider();
        layoutProvider.layout(parent, new BasicProgressMonitor());
        // Assert the ordering of n1 and ni.
        assertTrue("Node order of n1 and ni", n1.getY() < ni.getY());
        // Assert the ordering of p3 and p4. 
        assertTrue("p3 above p4", p3.getY() < p4.getY());
        
        // Check for port model order, which should result in the same graph.

        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_PORT_MODEL_ORDER, true);

        layoutProvider.layout(parent, new BasicProgressMonitor());
        // Assert the ordering of n1 and ni.
        assertTrue("Node order of n1 and ni", n1.getY() < ni.getY());
        // Assert the ordering of p3 and p4. 
        assertTrue("p3 above p4", p3.getY() < p4.getY());
    }

    @Test
    public void testOutgoingWestIncomingWest() {
        // |-----p1
        // |       n1
        // |  |->p2
        // |  |
        // |  |------p3
        // |           n2
        // |-------->p4

        ElkNode parent = ElkGraphUtil.createGraph();
        ElkNode n1 = ElkGraphUtil.createNode(parent);
        ElkNode n2 = ElkGraphUtil.createNode(parent);
        n1.setDimensions(20, 20);
        n2.setDimensions(20, 20);

        parent.setProperty(CoreOptions.ALGORITHM, LayeredOptions.ALGORITHM_ID);
        // Model order configuration.
        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY, OrderingStrategy.PREFER_EDGES);
        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_PORT_MODEL_ORDER, false);
        // No crossing minimization, hence only model order is used.
        parent.setProperty(LayeredOptions.CROSSING_MINIMIZATION_STRATEGY, CrossingMinimizationStrategy.NONE);
        parent.setProperty(LayeredOptions.CROSSING_MINIMIZATION_GREEDY_SWITCH_TYPE, GreedySwitchType.OFF);
        // For presentation.
        parent.setProperty(CoreOptions.DIRECTION, Direction.RIGHT);
        parent.setProperty(CoreOptions.PADDING, new ElkPadding(0.0));
        parent.setProperty(CoreOptions.SPACING_NODE_NODE, 10.0);
        parent.setProperty(LayeredOptions.SPACING_NODE_NODE_BETWEEN_LAYERS, 20.0);
        n1.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        n2.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        // Make sure that the order is correct.
        parent.setProperty(LayeredOptions.CYCLE_BREAKING_STRATEGY, CycleBreakingStrategy.MODEL_ORDER);
        
        ElkPort p1 = ElkGraphUtil.createPort(n1);
        p1.setProperty(LayeredOptions.PORT_SIDE, PortSide.WEST);
        ElkPort p2 = ElkGraphUtil.createPort(n1);
        p2.setProperty(LayeredOptions.PORT_SIDE, PortSide.WEST);
        ElkPort p3 = ElkGraphUtil.createPort(n2);
        p3.setProperty(LayeredOptions.PORT_SIDE, PortSide.WEST);
        ElkPort p4 = ElkGraphUtil.createPort(n2);
        p4.setProperty(LayeredOptions.PORT_SIDE, PortSide.WEST);

        ElkGraphUtil.createSimpleEdge(p1, p4);
        ElkGraphUtil.createSimpleEdge(p2, p3);
        
        LayeredLayoutProvider layoutProvider = new LayeredLayoutProvider();
        layoutProvider.layout(parent, new BasicProgressMonitor());
        // Assert the ordering of p1 and p2.
        assertTrue("p1 above p2", p1.getY() < p2.getY());
        // Assert the ordering of p3 and p4. 
        assertTrue("p3 above p4", p3.getY() < p4.getY());
        
        // Check for port model order, which should result in the same graph.

        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_PORT_MODEL_ORDER, true);

        layoutProvider.layout(parent, new BasicProgressMonitor());
        // Assert the ordering of p1 and p2.
        assertTrue("p1 above p2", p1.getY() < p2.getY());
        // Assert the ordering of p3 and p4. 
        assertTrue("p3 above p4", p3.getY() < p4.getY());
    }

    @Test
    public void testOutgoingWestIncomingWestMultipleNodes() {
        // |-p1:n1
        // |       
        // |-------->p4
        //             n2
        // |>p2:ni |-p3
        // |       |
        // |-------|

        ElkNode parent = ElkGraphUtil.createGraph();
        ElkNode n1 = ElkGraphUtil.createNode(parent);
        ElkNode ni = ElkGraphUtil.createNode(parent);
        ElkNode n2 = ElkGraphUtil.createNode(parent);
        n1.setDimensions(20, 20);
        ni.setDimensions(20, 20);
        n2.setDimensions(20, 20);

        parent.setProperty(CoreOptions.ALGORITHM, LayeredOptions.ALGORITHM_ID);
        // Model order configuration.
        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY, OrderingStrategy.PREFER_EDGES);
        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_PORT_MODEL_ORDER, false);
        // No crossing minimization, hence only model order is used.
        parent.setProperty(LayeredOptions.CROSSING_MINIMIZATION_STRATEGY, CrossingMinimizationStrategy.NONE);
        parent.setProperty(LayeredOptions.CROSSING_MINIMIZATION_GREEDY_SWITCH_TYPE, GreedySwitchType.OFF);
        // For presentation.
        parent.setProperty(CoreOptions.DIRECTION, Direction.RIGHT);
        parent.setProperty(CoreOptions.PADDING, new ElkPadding(0.0));
        parent.setProperty(CoreOptions.SPACING_NODE_NODE, 10.0);
        parent.setProperty(LayeredOptions.SPACING_NODE_NODE_BETWEEN_LAYERS, 20.0);
        n1.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        ni.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        n2.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        // Make sure that the order is correct.
        parent.setProperty(LayeredOptions.CYCLE_BREAKING_STRATEGY, CycleBreakingStrategy.MODEL_ORDER);
        
        ElkPort p1 = ElkGraphUtil.createPort(n1);
        p1.setProperty(LayeredOptions.PORT_SIDE, PortSide.WEST);
        ElkPort p2 = ElkGraphUtil.createPort(ni);
        p2.setProperty(LayeredOptions.PORT_SIDE, PortSide.WEST);
        ElkPort p3 = ElkGraphUtil.createPort(n2);
        p3.setProperty(LayeredOptions.PORT_SIDE, PortSide.WEST);
        ElkPort p4 = ElkGraphUtil.createPort(n2);
        p4.setProperty(LayeredOptions.PORT_SIDE, PortSide.WEST);

        ElkGraphUtil.createSimpleEdge(p1, p4);
        ElkGraphUtil.createSimpleEdge(p2, p3);
        
        LayeredLayoutProvider layoutProvider = new LayeredLayoutProvider();
        layoutProvider.layout(parent, new BasicProgressMonitor());
        // Assert the ordering of n1 and ni.
        assertTrue("Node order of n1 and ni", n1.getY() < ni.getY());
        // Assert the ordering of p3 and p4. 
        assertTrue("p4 above p3", p4.getY() < p3.getY());
        
        // Check for port model order, which should result in the same graph.

        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_PORT_MODEL_ORDER, true);

        layoutProvider.layout(parent, new BasicProgressMonitor());
        // Assert the ordering of n1 and ni.
        assertTrue("Node order of n1 and ni", n1.getY() < ni.getY());
        // Assert the ordering of p3 and p4. 
        assertTrue("p4 above p3", p4.getY() < p3.getY());
    }

    @Test
    public void testOutgoingWestIncomingNorth() {
        // |-----p1
        // |       n1
        // |  |->p2
        // |  |
        // |  |----------|
        // |             |
        // |---------v   |
        //           p4  p3
        //             n2

        ElkNode parent = ElkGraphUtil.createGraph();
        ElkNode n1 = ElkGraphUtil.createNode(parent);
        ElkNode n2 = ElkGraphUtil.createNode(parent);
        n1.setDimensions(20, 20);
        n2.setDimensions(20, 20);

        parent.setProperty(CoreOptions.ALGORITHM, LayeredOptions.ALGORITHM_ID);
        // Model order configuration.
        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY, OrderingStrategy.PREFER_EDGES);
        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_PORT_MODEL_ORDER, false);
        // No crossing minimization, hence only model order is used.
        parent.setProperty(LayeredOptions.CROSSING_MINIMIZATION_STRATEGY, CrossingMinimizationStrategy.NONE);
        parent.setProperty(LayeredOptions.CROSSING_MINIMIZATION_GREEDY_SWITCH_TYPE, GreedySwitchType.OFF);
        // For presentation.
        parent.setProperty(CoreOptions.DIRECTION, Direction.RIGHT);
        parent.setProperty(CoreOptions.PADDING, new ElkPadding(0.0));
        parent.setProperty(CoreOptions.SPACING_NODE_NODE, 10.0);
        parent.setProperty(LayeredOptions.SPACING_NODE_NODE_BETWEEN_LAYERS, 20.0);
        n1.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        n2.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        // Make sure that the order is correct.
        parent.setProperty(LayeredOptions.CYCLE_BREAKING_STRATEGY, CycleBreakingStrategy.MODEL_ORDER);
        
        ElkPort p1 = ElkGraphUtil.createPort(n1);
        p1.setProperty(LayeredOptions.PORT_SIDE, PortSide.WEST);
        ElkPort p2 = ElkGraphUtil.createPort(n1);
        p2.setProperty(LayeredOptions.PORT_SIDE, PortSide.WEST);
        ElkPort p3 = ElkGraphUtil.createPort(n2);
        p3.setProperty(LayeredOptions.PORT_SIDE, PortSide.NORTH);
        ElkPort p4 = ElkGraphUtil.createPort(n2);
        p4.setProperty(LayeredOptions.PORT_SIDE, PortSide.NORTH);

        ElkGraphUtil.createSimpleEdge(p1, p4);
        ElkGraphUtil.createSimpleEdge(p2, p3);
        
        LayeredLayoutProvider layoutProvider = new LayeredLayoutProvider();
        layoutProvider.layout(parent, new BasicProgressMonitor());
        // Assert the ordering of p1 and p2.
        assertTrue("p1 above p2", p1.getY() < p2.getY());
        // Assert the ordering of p3 and p4. 
        assertTrue("p4 before p3", p4.getX() < p3.getX());
        
        // Check for port model order, which should result in the same graph.

        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_PORT_MODEL_ORDER, true);

        layoutProvider.layout(parent, new BasicProgressMonitor());
        // Assert the ordering of p1 and p2.
        assertTrue("p1 above p2", p1.getY() < p2.getY());
        // Assert the ordering of p3 and p4. 
        assertTrue("p4 before p3", p4.getX() < p3.getX());
    }

    @Test
    public void testOutgoingWestIncomingNorthMultipleNodes() {
        // |-p1:n1
        // |       
        // |-----------|
        //             |
        // |>p2:ni     |
        // |           |
        // |--------|  v
        //         p3  p4
        //           n2

        ElkNode parent = ElkGraphUtil.createGraph();
        ElkNode n1 = ElkGraphUtil.createNode(parent);
        ElkNode ni = ElkGraphUtil.createNode(parent);
        ElkNode n2 = ElkGraphUtil.createNode(parent);
        n1.setDimensions(20, 20);
        ni.setDimensions(20, 20);
        n2.setDimensions(20, 20);

        parent.setProperty(CoreOptions.ALGORITHM, LayeredOptions.ALGORITHM_ID);
        // Model order configuration.
        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY, OrderingStrategy.PREFER_EDGES);
        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_PORT_MODEL_ORDER, false);
        // No crossing minimization, hence only model order is used.
        parent.setProperty(LayeredOptions.CROSSING_MINIMIZATION_STRATEGY, CrossingMinimizationStrategy.NONE);
        parent.setProperty(LayeredOptions.CROSSING_MINIMIZATION_GREEDY_SWITCH_TYPE, GreedySwitchType.OFF);
        // For presentation.
        parent.setProperty(CoreOptions.DIRECTION, Direction.RIGHT);
        parent.setProperty(CoreOptions.PADDING, new ElkPadding(0.0));
        parent.setProperty(CoreOptions.SPACING_NODE_NODE, 10.0);
        parent.setProperty(LayeredOptions.SPACING_NODE_NODE_BETWEEN_LAYERS, 20.0);
        n1.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        ni.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        n2.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        // Make sure that the order is correct.
        parent.setProperty(LayeredOptions.CYCLE_BREAKING_STRATEGY, CycleBreakingStrategy.MODEL_ORDER);
        
        ElkPort p1 = ElkGraphUtil.createPort(n1);
        p1.setProperty(LayeredOptions.PORT_SIDE, PortSide.WEST);
        ElkPort p2 = ElkGraphUtil.createPort(ni);
        p2.setProperty(LayeredOptions.PORT_SIDE, PortSide.WEST);
        ElkPort p3 = ElkGraphUtil.createPort(n2);
        p3.setProperty(LayeredOptions.PORT_SIDE, PortSide.NORTH);
        ElkPort p4 = ElkGraphUtil.createPort(n2);
        p4.setProperty(LayeredOptions.PORT_SIDE, PortSide.NORTH);

        ElkGraphUtil.createSimpleEdge(p1, p4);
        ElkGraphUtil.createSimpleEdge(p2, p3);
        
        LayeredLayoutProvider layoutProvider = new LayeredLayoutProvider();
        layoutProvider.layout(parent, new BasicProgressMonitor());
        // Assert the ordering of n1 and ni.
        assertTrue("Node order of n1 and ni", n1.getY() < ni.getY());
        // Assert the ordering of p3 and p4. 
        assertTrue("p3 before p4", p3.getX() < p4.getX());
        
        // Check for port model order, which should result in the same graph.

        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_PORT_MODEL_ORDER, true);

        layoutProvider.layout(parent, new BasicProgressMonitor());
        // Assert the ordering of n1 and ni.
        assertTrue("Node order of n1 and ni", n1.getY() < ni.getY());
        // Assert the ordering of p3 and p4. 
        assertTrue("p3 before p4", p3.getX() < p4.getX());
    }

    @Test
    public void testOutgoingWestIncomingSouth() {
        // |-----p1
        // |       n1
        // |  |->p2    n2
        // |  |      p3  p4
        // |  |-------|  A
        // |             |
        // |-------------|

        ElkNode parent = ElkGraphUtil.createGraph();
        ElkNode n1 = ElkGraphUtil.createNode(parent);
        ElkNode n2 = ElkGraphUtil.createNode(parent);
        n1.setDimensions(20, 20);
        n2.setDimensions(20, 20);

        parent.setProperty(CoreOptions.ALGORITHM, LayeredOptions.ALGORITHM_ID);
        // Model order configuration.
        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY, OrderingStrategy.PREFER_EDGES);
        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_PORT_MODEL_ORDER, false);
        // No crossing minimization, hence only model order is used.
        parent.setProperty(LayeredOptions.CROSSING_MINIMIZATION_STRATEGY, CrossingMinimizationStrategy.NONE);
        parent.setProperty(LayeredOptions.CROSSING_MINIMIZATION_GREEDY_SWITCH_TYPE, GreedySwitchType.OFF);
        // For presentation.
        parent.setProperty(CoreOptions.DIRECTION, Direction.RIGHT);
        parent.setProperty(CoreOptions.PADDING, new ElkPadding(0.0));
        parent.setProperty(CoreOptions.SPACING_NODE_NODE, 10.0);
        parent.setProperty(LayeredOptions.SPACING_NODE_NODE_BETWEEN_LAYERS, 20.0);
        n1.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        n2.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        // Make sure that the order is correct.
        parent.setProperty(LayeredOptions.CYCLE_BREAKING_STRATEGY, CycleBreakingStrategy.MODEL_ORDER);
        
        ElkPort p1 = ElkGraphUtil.createPort(n1);
        p1.setProperty(LayeredOptions.PORT_SIDE, PortSide.WEST);
        ElkPort p2 = ElkGraphUtil.createPort(n1);
        p2.setProperty(LayeredOptions.PORT_SIDE, PortSide.WEST);
        ElkPort p3 = ElkGraphUtil.createPort(n2);
        p3.setProperty(LayeredOptions.PORT_SIDE, PortSide.SOUTH);
        ElkPort p4 = ElkGraphUtil.createPort(n2);
        p4.setProperty(LayeredOptions.PORT_SIDE, PortSide.SOUTH);

        ElkGraphUtil.createSimpleEdge(p1, p4);
        ElkGraphUtil.createSimpleEdge(p2, p3);
        
        LayeredLayoutProvider layoutProvider = new LayeredLayoutProvider();
        layoutProvider.layout(parent, new BasicProgressMonitor());
        // Assert the ordering of p1 and p2.
        assertTrue("p1 above p2", p1.getY() < p2.getY());
        // Assert the ordering of p3 and p4. 
        assertTrue("p3 before p4", p3.getX() < p4.getX());
        
        // Check for port model order, which should result in the same graph.

        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_PORT_MODEL_ORDER, true);

        layoutProvider.layout(parent, new BasicProgressMonitor());
        // Assert the ordering of p1 and p2.
        assertTrue("p1 above p2", p1.getY() < p2.getY());
        // Assert the ordering of p3 and p4. 
        assertTrue("p3 before p4", p3.getX() < p4.getX());
    }

    @Test
    public void testOutgoingWestIncomingSouthMultipleNodes() {
        // |-----p1:n1    n2
        // |            p4  p3
        // |            A   |
        // |------------|   |
        //                  |
        // |-----p2:ni      |
        // |                |
        // |----------------|

        ElkNode parent = ElkGraphUtil.createGraph();
        ElkNode n1 = ElkGraphUtil.createNode(parent);
        ElkNode ni = ElkGraphUtil.createNode(parent);
        ElkNode n2 = ElkGraphUtil.createNode(parent);
        n1.setDimensions(20, 20);
        ni.setDimensions(20, 20);
        n2.setDimensions(20, 20);

        parent.setProperty(CoreOptions.ALGORITHM, LayeredOptions.ALGORITHM_ID);
        // Model order configuration.
        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY, OrderingStrategy.PREFER_EDGES);
        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_PORT_MODEL_ORDER, false);
        // No crossing minimization, hence only model order is used.
        parent.setProperty(LayeredOptions.CROSSING_MINIMIZATION_STRATEGY, CrossingMinimizationStrategy.NONE);
        parent.setProperty(LayeredOptions.CROSSING_MINIMIZATION_GREEDY_SWITCH_TYPE, GreedySwitchType.OFF);
        // For presentation.
        parent.setProperty(CoreOptions.DIRECTION, Direction.RIGHT);
        parent.setProperty(CoreOptions.PADDING, new ElkPadding(0.0));
        parent.setProperty(CoreOptions.SPACING_NODE_NODE, 10.0);
        parent.setProperty(LayeredOptions.SPACING_NODE_NODE_BETWEEN_LAYERS, 20.0);
        n1.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        ni.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        n2.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        // Make sure that the order is correct.
        parent.setProperty(LayeredOptions.CYCLE_BREAKING_STRATEGY, CycleBreakingStrategy.MODEL_ORDER);
        
        ElkPort p1 = ElkGraphUtil.createPort(n1);
        p1.setProperty(LayeredOptions.PORT_SIDE, PortSide.WEST);
        ElkPort p2 = ElkGraphUtil.createPort(ni);
        p2.setProperty(LayeredOptions.PORT_SIDE, PortSide.WEST);
        ElkPort p3 = ElkGraphUtil.createPort(n2);
        p3.setProperty(LayeredOptions.PORT_SIDE, PortSide.SOUTH);
        ElkPort p4 = ElkGraphUtil.createPort(n2);
        p4.setProperty(LayeredOptions.PORT_SIDE, PortSide.SOUTH);

        ElkGraphUtil.createSimpleEdge(p1, p4);
        ElkGraphUtil.createSimpleEdge(p2, p3);
        
        LayeredLayoutProvider layoutProvider = new LayeredLayoutProvider();
        layoutProvider.layout(parent, new BasicProgressMonitor());
        // Assert the ordering of n1 and ni.
        assertTrue("Node order of n1 and ni", n1.getY() < ni.getY());
        // Assert the ordering of p3 and p4. 
        assertTrue("p4 before p3", p4.getX() < p3.getX());
        
        // Check for port model order, which should result in the same graph.

        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_PORT_MODEL_ORDER, true);

        layoutProvider.layout(parent, new BasicProgressMonitor());
        // Assert the ordering of n1 and ni.
        assertTrue("Node order of n1 and ni", n1.getY() < ni.getY());
        // Assert the ordering of p3 and p4. 
        assertTrue("p4 before p3", p4.getX() < p3.getX());
    }

    @Test
    public void testOutgoingWestIncomingEast() {
        // |-----p1     p4<----|
        // |       n1 n2       |
        // |  |->p2     p3--|  |
        // |  |             |  |
        // |  |-------------|  |
        // |                   |
        // |-------------------|

        ElkNode parent = ElkGraphUtil.createGraph();
        ElkNode n1 = ElkGraphUtil.createNode(parent);
        ElkNode n2 = ElkGraphUtil.createNode(parent);
        n1.setDimensions(20, 20);
        n2.setDimensions(20, 20);

        parent.setProperty(CoreOptions.ALGORITHM, LayeredOptions.ALGORITHM_ID);
        // Model order configuration.
        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY, OrderingStrategy.PREFER_EDGES);
        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_PORT_MODEL_ORDER, false);
        // No crossing minimization, hence only model order is used.
        parent.setProperty(LayeredOptions.CROSSING_MINIMIZATION_STRATEGY, CrossingMinimizationStrategy.NONE);
        parent.setProperty(LayeredOptions.CROSSING_MINIMIZATION_GREEDY_SWITCH_TYPE, GreedySwitchType.OFF);
        // For presentation.
        parent.setProperty(CoreOptions.DIRECTION, Direction.RIGHT);
        parent.setProperty(CoreOptions.PADDING, new ElkPadding(0.0));
        parent.setProperty(CoreOptions.SPACING_NODE_NODE, 10.0);
        parent.setProperty(LayeredOptions.SPACING_NODE_NODE_BETWEEN_LAYERS, 20.0);
        n1.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        n2.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        // Make sure that the order is correct.
        parent.setProperty(LayeredOptions.CYCLE_BREAKING_STRATEGY, CycleBreakingStrategy.MODEL_ORDER);
        
        ElkPort p1 = ElkGraphUtil.createPort(n1);
        p1.setProperty(LayeredOptions.PORT_SIDE, PortSide.WEST);
        ElkPort p2 = ElkGraphUtil.createPort(n1);
        p2.setProperty(LayeredOptions.PORT_SIDE, PortSide.WEST);
        ElkPort p3 = ElkGraphUtil.createPort(n2);
        p3.setProperty(LayeredOptions.PORT_SIDE, PortSide.EAST);
        ElkPort p4 = ElkGraphUtil.createPort(n2);
        p4.setProperty(LayeredOptions.PORT_SIDE, PortSide.EAST);

        ElkGraphUtil.createSimpleEdge(p1, p4);
        ElkGraphUtil.createSimpleEdge(p2, p3);
        
        LayeredLayoutProvider layoutProvider = new LayeredLayoutProvider();
        layoutProvider.layout(parent, new BasicProgressMonitor());
        // Assert the ordering of p1 and p2.
        assertTrue("p1 above p2", p1.getY() < p2.getY());
        // Assert the ordering of p3 and p4. 
        assertTrue("p4 above p3", p4.getY() < p3.getY());
        
        // Check for port model order, which should result in the same graph.

        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_PORT_MODEL_ORDER, true);

        layoutProvider.layout(parent, new BasicProgressMonitor());
        // Assert the ordering of p1 and p2.
        assertTrue("p1 above p2", p1.getY() < p2.getY());
        // Assert the ordering of p3 and p4. 
        assertTrue("p4 above p3", p4.getY() < p3.getY());
    }

    @Test
    public void testOutgoingWestIncomingEastMultipleNodes() {
        //                 p3-----|
        //               n2       |
        //    |--p1:n1     p4<-|  |
        //    |                |  |
        //    |----------------|  |
        //                        |
        //    |->p2:ni            |
        //    |                   |
        //    |-------------------|

        ElkNode parent = ElkGraphUtil.createGraph();
        ElkNode n1 = ElkGraphUtil.createNode(parent);
        ElkNode ni = ElkGraphUtil.createNode(parent);
        ElkNode n2 = ElkGraphUtil.createNode(parent);
        n1.setDimensions(20, 20);
        ni.setDimensions(20, 20);
        n2.setDimensions(20, 20);

        parent.setProperty(CoreOptions.ALGORITHM, LayeredOptions.ALGORITHM_ID);
        // Model order configuration.
        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY, OrderingStrategy.PREFER_EDGES);
        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_PORT_MODEL_ORDER, false);
        // No crossing minimization, hence only model order is used.
        parent.setProperty(LayeredOptions.CROSSING_MINIMIZATION_STRATEGY, CrossingMinimizationStrategy.NONE);
        parent.setProperty(LayeredOptions.CROSSING_MINIMIZATION_GREEDY_SWITCH_TYPE, GreedySwitchType.OFF);
        // For presentation.
        parent.setProperty(CoreOptions.DIRECTION, Direction.RIGHT);
        parent.setProperty(CoreOptions.PADDING, new ElkPadding(0.0));
        parent.setProperty(CoreOptions.SPACING_NODE_NODE, 10.0);
        parent.setProperty(LayeredOptions.SPACING_NODE_NODE_BETWEEN_LAYERS, 20.0);
        n1.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        ni.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        n2.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        // Make sure that the order is correct.
        parent.setProperty(LayeredOptions.CYCLE_BREAKING_STRATEGY, CycleBreakingStrategy.MODEL_ORDER);
        
        ElkPort p1 = ElkGraphUtil.createPort(n1);
        p1.setProperty(LayeredOptions.PORT_SIDE, PortSide.WEST);
        ElkPort p2 = ElkGraphUtil.createPort(ni);
        p2.setProperty(LayeredOptions.PORT_SIDE, PortSide.WEST);
        ElkPort p3 = ElkGraphUtil.createPort(n2);
        p3.setProperty(LayeredOptions.PORT_SIDE, PortSide.EAST);
        ElkPort p4 = ElkGraphUtil.createPort(n2);
        p4.setProperty(LayeredOptions.PORT_SIDE, PortSide.EAST);

        ElkGraphUtil.createSimpleEdge(p1, p4);
        ElkGraphUtil.createSimpleEdge(p2, p3);
        
        LayeredLayoutProvider layoutProvider = new LayeredLayoutProvider();
        layoutProvider.layout(parent, new BasicProgressMonitor());
        // Assert the ordering of n1 and ni.
        assertTrue("Node order of n1 and ni", n1.getY() < ni.getY());
        // Assert the ordering of p3 and p4. 
        assertTrue("p3 above p4", p3.getY() < p4.getY());
        
        // Check for port model order, which should result in the same graph.

        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_PORT_MODEL_ORDER, true);

        layoutProvider.layout(parent, new BasicProgressMonitor());
        // Assert the ordering of n1 and ni.
        assertTrue("Node order of n1 and ni", n1.getY() < ni.getY());
        // Assert the ordering of p3 and p4. 
        assertTrue("p3 above p4", p3.getY() < p4.getY());
    }
}