/*******************************************************************************
 * Copyright (c) 2023 Kiel University and others.
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
package org.eclipse.elk.alg.layered.issues;

import static org.junit.Assert.assertEquals;

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
import org.eclipse.elk.core.util.BasicProgressMonitor;
import org.eclipse.elk.graph.ElkEdge;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.util.ElkGraphUtil;
import org.junit.BeforeClass;
import org.junit.Test;

/**
 * Test whether the model order layer assignment together with model order node promotion works
 * for some selected examples.
 *
 */
public class Issue871Test {

    @BeforeClass
    public static void init() {
        PlainJavaInitialization.initializePlainJavaLayout();
    }

    @Test
    public void testFeedbackEdgeBelow() {
        // Desired layout
        // n1 -> n2 -> n4 -
        //    -> n3       |
        //   |____________|
        ElkNode parent = ElkGraphUtil.createGraph();
        ElkNode n1 = ElkGraphUtil.createNode(parent);
        ElkNode n2 = ElkGraphUtil.createNode(parent);
        ElkNode n3 = ElkGraphUtil.createNode(parent);
        ElkNode n4 = ElkGraphUtil.createNode(parent);
        n1.setDimensions(30, 30);
        n2.setDimensions(30, 30);
        n3.setDimensions(30, 30);
        n4.setDimensions(30, 30);

        parent.setProperty(CoreOptions.ALGORITHM, LayeredOptions.ALGORITHM_ID);
        parent.setProperty(CoreOptions.DIRECTION, Direction.RIGHT);
        parent.setProperty(LayeredOptions.CYCLE_BREAKING_STRATEGY, CycleBreakingStrategy.MODEL_ORDER);
        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY, OrderingStrategy.PREFER_EDGES);
        parent.setProperty(LayeredOptions.CROSSING_MINIMIZATION_STRATEGY, CrossingMinimizationStrategy.NONE);
        parent.setProperty(LayeredOptions.CROSSING_MINIMIZATION_GREEDY_SWITCH_TYPE, GreedySwitchType.OFF);
        parent.setProperty(CoreOptions.PADDING, new ElkPadding(0.0));
        parent.setProperty(CoreOptions.SPACING_NODE_NODE, 10.0);
        parent.setProperty(LayeredOptions.SPACING_NODE_NODE_BETWEEN_LAYERS, 20.0);
        parent.setProperty(LayeredOptions.FEEDBACK_EDGES, true);

        ElkGraphUtil.createSimpleEdge(n1, n2);
        ElkGraphUtil.createSimpleEdge(n1, n3);
        ElkGraphUtil.createSimpleEdge(n2, n4);
        ElkGraphUtil.createSimpleEdge(n4, n3); // Feedback edge
        ElkGraphUtil.createLabel("n1", n1);
        ElkGraphUtil.createLabel("n2", n2);
        ElkGraphUtil.createLabel("n3", n3);
        ElkGraphUtil.createLabel("n4", n4);
        
        LayeredLayoutProvider layoutProvider = new LayeredLayoutProvider();
        layoutProvider.layout(parent, new BasicProgressMonitor());

        // Check whether n4 aligns with n2.
        assertEquals("n4 and n2 should align", n2.getY(), n4.getY(), 0.1);
    }

    @Test
    public void testFeedbackEdgeBasic() {
        // Desired layout
        // n1 -> n2 -> n3-
        //   |____________|
        ElkNode parent = ElkGraphUtil.createGraph();
        ElkNode n1 = ElkGraphUtil.createNode(parent);
        ElkNode n2 = ElkGraphUtil.createNode(parent);
        ElkNode n3 = ElkGraphUtil.createNode(parent);
        n1.setDimensions(30, 30);
        n2.setDimensions(30, 30);
        n3.setDimensions(30, 30);

        parent.setProperty(CoreOptions.ALGORITHM, LayeredOptions.ALGORITHM_ID);
        parent.setProperty(CoreOptions.DIRECTION, Direction.RIGHT);
        parent.setProperty(LayeredOptions.CYCLE_BREAKING_STRATEGY, CycleBreakingStrategy.MODEL_ORDER);
        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY, OrderingStrategy.PREFER_EDGES);
        parent.setProperty(LayeredOptions.CROSSING_MINIMIZATION_STRATEGY, CrossingMinimizationStrategy.NONE);
        parent.setProperty(LayeredOptions.CROSSING_MINIMIZATION_GREEDY_SWITCH_TYPE, GreedySwitchType.OFF);
        parent.setProperty(CoreOptions.PADDING, new ElkPadding(0.0));
        parent.setProperty(CoreOptions.SPACING_NODE_NODE, 10.0);
        parent.setProperty(LayeredOptions.SPACING_NODE_NODE_BETWEEN_LAYERS, 20.0);
        parent.setProperty(LayeredOptions.FEEDBACK_EDGES, true);

        ElkGraphUtil.createSimpleEdge(n1, n2);
        ElkGraphUtil.createSimpleEdge(n2, n3);
        ElkGraphUtil.createSimpleEdge(n3, n2); // Feedback edge
        ElkGraphUtil.createLabel("n1", n1);
        ElkGraphUtil.createLabel("n2", n2);
        ElkGraphUtil.createLabel("n3", n3);
        
        LayeredLayoutProvider layoutProvider = new LayeredLayoutProvider();
        layoutProvider.layout(parent, new BasicProgressMonitor());

        // Check whether n3 aligns with n2.
        assertEquals("n3 and n2 should align", n2.getY(), n3.getY(), 0.1);
    }

    @Test
    public void testNoFeedbackEdgesStillWorking() {
        // Desired layout
        // n1 -1-> n2 ------>
        //           <- n3 -> n4
        //   --2------------>
        ElkNode parent = ElkGraphUtil.createGraph();
        ElkNode n1 = ElkGraphUtil.createNode(parent);
        ElkNode n2 = ElkGraphUtil.createNode(parent);
        ElkNode n3 = ElkGraphUtil.createNode(parent);
        ElkNode n4 = ElkGraphUtil.createNode(parent);
        n1.setDimensions(30, 30);
        n2.setDimensions(30, 30);
        n3.setDimensions(30, 30);
        n4.setDimensions(30, 30);

        parent.setProperty(CoreOptions.ALGORITHM, LayeredOptions.ALGORITHM_ID);
        parent.setProperty(CoreOptions.DIRECTION, Direction.RIGHT);
        parent.setProperty(LayeredOptions.CYCLE_BREAKING_STRATEGY, CycleBreakingStrategy.MODEL_ORDER);
        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY, OrderingStrategy.PREFER_EDGES);
        parent.setProperty(LayeredOptions.CROSSING_MINIMIZATION_STRATEGY, CrossingMinimizationStrategy.NONE);
        parent.setProperty(LayeredOptions.CROSSING_MINIMIZATION_GREEDY_SWITCH_TYPE, GreedySwitchType.OFF);
        parent.setProperty(CoreOptions.PADDING, new ElkPadding(0.0));
        parent.setProperty(CoreOptions.SPACING_NODE_NODE, 10.0);
        parent.setProperty(LayeredOptions.SPACING_NODE_NODE_BETWEEN_LAYERS, 20.0);

        ElkEdge e1 = ElkGraphUtil.createSimpleEdge(n1, n2);
        ElkGraphUtil.createLabel("1", e1);
        ElkEdge e2 = ElkGraphUtil.createSimpleEdge(n1, n4);
        ElkGraphUtil.createLabel("2", e2);
        ElkGraphUtil.createSimpleEdge(n2, n4);
        ElkGraphUtil.createSimpleEdge(n3, n2);
        ElkGraphUtil.createSimpleEdge(n3, n4);
        ElkGraphUtil.createLabel("n1", n1);
        ElkGraphUtil.createLabel("n2", n2);
        ElkGraphUtil.createLabel("n3", n3);
        ElkGraphUtil.createLabel("n4", n4);
        
        LayeredLayoutProvider layoutProvider = new LayeredLayoutProvider();
        layoutProvider.layout(parent, new BasicProgressMonitor());

        // Check whether nodes are on the expected positions.
        assertEquals("n1 X coordinate is wrong", 0, n1.getX(), 0.1);
        assertEquals("n1 Y coordinate is wrong", 31, n1.getY(), 0.1);
        assertEquals("n2 X coordinate is wrong", 70.0, n2.getX(), 0.1);
        assertEquals("n2 Y coordinate is wrong", 6, n2.getY(), 0.1);
        assertEquals("n3 X coordinate is wrong", 120, n3.getX(), 0.1);
        assertEquals("n3 Y coordinate is wrong", 11, n3.getY(), 0.1);
        assertEquals("n4 X coordinate is wrong", 170, n4.getX(), 0.1);
        assertEquals("n4 Y coordinate is wrong", 11, n4.getY(), 0.1);
    }
}
