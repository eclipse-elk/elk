/*******************************************************************************
 * Copyright (c) 2026 Kiel University and others.
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
package org.eclipse.elk.alg.layered;

import static org.junit.Assert.assertThrows;
import static org.junit.Assert.assertTrue;

import org.eclipse.elk.alg.layered.options.CycleBreakingStrategy;
import org.eclipse.elk.alg.layered.options.GroupOrderStrategy;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.layered.options.OrderingStrategy;
import org.eclipse.elk.alg.test.PlainJavaInitialization;
import org.eclipse.elk.core.RecursiveGraphLayoutEngine;
import org.eclipse.elk.core.UnsupportedConfigurationException;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.options.HierarchyHandling;
import org.eclipse.elk.core.util.BasicProgressMonitor;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.properties.IProperty;
import org.eclipse.elk.graph.util.ElkGraphUtil;
import org.junit.BeforeClass;
import org.junit.Test;

/**
 * Tests that ELK Layered refuses model order strategies configured on a nested graph that is laid out as part of its
 * parent's hierarchy ({@link HierarchyHandling#INCLUDE_CHILDREN}). The model order information these strategies read
 * is only computed for the top-level graph of such a hierarchy, and model order strategies are not supported with
 * {@link HierarchyHandling#INCLUDE_CHILDREN}. Before the check, these configurations failed with a
 * {@link NullPointerException} from deep inside a layout phase.
 */
public class NestedModelOrderConfigurationTest {

    private static final String NESTED_ID = "p";

    @BeforeClass
    public static void init() {
        PlainJavaInitialization.initializePlainJavaLayout();
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Graphs

    /** A root graph with the given hierarchy handling. */
    private static ElkNode root(final HierarchyHandling hierarchyHandling) {
        ElkNode root = ElkGraphUtil.createGraph();
        root.setProperty(CoreOptions.ALGORITHM, LayeredOptions.ALGORITHM_ID);
        root.setProperty(CoreOptions.HIERARCHY_HANDLING, hierarchyHandling);
        return root;
    }

    /** A node {@code p} in the given graph with {@code nodes} children {@code n0}, {@code n1}, .... */
    private static ElkNode nested(final ElkNode graph, final int nodes) {
        ElkNode nested = node(graph, NESTED_ID);
        for (int i = 0; i < nodes; i++) {
            node(nested, "n" + i);
        }
        return nested;
    }

    private static ElkNode node(final ElkNode parent, final String id) {
        ElkNode node = ElkGraphUtil.createNode(parent);
        node.setIdentifier(id);
        node.setDimensions(30, 30);
        return node;
    }

    /** Connects the first two children of the given node with {@code n0 -> n1} and, if requested, {@code n1 -> n0}. */
    private static void connect(final ElkNode parent, final boolean cycle) {
        ElkNode n0 = parent.getChildren().get(0);
        ElkNode n1 = parent.getChildren().get(1);
        ElkGraphUtil.createSimpleEdge(n0, n1);
        if (cycle) {
            ElkGraphUtil.createSimpleEdge(n1, n0);
        }
    }

    private static void layout(final ElkNode graph) {
        new RecursiveGraphLayoutEngine().layout(graph, new BasicProgressMonitor());
    }

    private static void assertRefused(final ElkNode graph, final IProperty<?> option, final Object value) {
        UnsupportedConfigurationException e =
                assertThrows(UnsupportedConfigurationException.class, () -> layout(graph));
        String message = e.getMessage();
        assertTrue(message, message.startsWith("Node '" + NESTED_ID + "' sets " + option.getId() + " " + value));
        assertTrue(message, message.contains(HierarchyHandling.INCLUDE_CHILDREN.toString()));
        assertTrue(message, message.contains(HierarchyHandling.SEPARATE_CHILDREN.toString()));
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Refused configurations

    @Test
    public void testModelOrderCycleBreakingRefused() {
        ElkNode root = root(HierarchyHandling.INCLUDE_CHILDREN);
        ElkNode nested = nested(root, 1);
        nested.setProperty(LayeredOptions.CYCLE_BREAKING_STRATEGY, CycleBreakingStrategy.MODEL_ORDER);
        assertRefused(root, LayeredOptions.CYCLE_BREAKING_STRATEGY, CycleBreakingStrategy.MODEL_ORDER);
    }

    @Test
    public void testOtherModelOrderCycleBreakersRefused() {
        for (CycleBreakingStrategy strategy : new CycleBreakingStrategy[] { CycleBreakingStrategy.GREEDY_MODEL_ORDER,
                CycleBreakingStrategy.SCC_CONNECTIVITY, CycleBreakingStrategy.SCC_NODE_TYPE }) {
            ElkNode root = root(HierarchyHandling.INCLUDE_CHILDREN);
            ElkNode nested = nested(root, 2);
            connect(nested, true);
            nested.setProperty(LayeredOptions.CYCLE_BREAKING_STRATEGY, strategy);
            assertRefused(root, LayeredOptions.CYCLE_BREAKING_STRATEGY, strategy);
        }
    }

    @Test
    public void testDepthFirstNodeOrderWithEnforcedGroupOrderRefused() {
        ElkNode root = root(HierarchyHandling.INCLUDE_CHILDREN);
        ElkNode nested = nested(root, 2);
        connect(nested, false);
        nested.setProperty(LayeredOptions.CYCLE_BREAKING_STRATEGY, CycleBreakingStrategy.DFS_NODE_ORDER);
        nested.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_GROUP_MODEL_ORDER_CB_GROUP_ORDER_STRATEGY,
                GroupOrderStrategy.ENFORCED);
        assertRefused(root, LayeredOptions.CYCLE_BREAKING_STRATEGY, CycleBreakingStrategy.DFS_NODE_ORDER);
    }

    @Test
    public void testBreadthFirstNodeOrderWithEnforcedGroupOrderRefused() {
        ElkNode root = root(HierarchyHandling.INCLUDE_CHILDREN);
        ElkNode nested = nested(root, 2);
        connect(nested, true);
        nested.setProperty(LayeredOptions.CYCLE_BREAKING_STRATEGY, CycleBreakingStrategy.BFS_NODE_ORDER);
        nested.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_GROUP_MODEL_ORDER_CB_GROUP_ORDER_STRATEGY,
                GroupOrderStrategy.ENFORCED);
        assertRefused(root, LayeredOptions.CYCLE_BREAKING_STRATEGY, CycleBreakingStrategy.BFS_NODE_ORDER);
    }

    @Test
    public void testConsiderModelOrderRefused() {
        for (OrderingStrategy strategy : new OrderingStrategy[] { OrderingStrategy.NODES_AND_EDGES,
                OrderingStrategy.PREFER_EDGES, OrderingStrategy.PREFER_NODES }) {
            ElkNode root = root(HierarchyHandling.INCLUDE_CHILDREN);
            ElkNode nested = nested(root, 2);
            nested.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY, strategy);
            assertRefused(root, LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY, strategy);
        }
    }

    /**
     * A single nested node never makes the model order readers compare two nodes, so this configuration used to lay
     * out. It is refused now as well, since the configuration is the same one that fails as soon as there are two.
     */
    @Test
    public void testConsiderModelOrderWithSingleNodeRefused() {
        ElkNode root = root(HierarchyHandling.INCLUDE_CHILDREN);
        ElkNode nested = nested(root, 1);
        nested.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY, OrderingStrategy.PREFER_NODES);
        assertRefused(root, LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY, OrderingStrategy.PREFER_NODES);
    }

    @Test
    public void testForceNodeModelOrderRefused() {
        ElkNode root = root(HierarchyHandling.INCLUDE_CHILDREN);
        ElkNode nested = nested(root, 2);
        nested.setProperty(LayeredOptions.CROSSING_MINIMIZATION_FORCE_NODE_MODEL_ORDER, true);
        assertRefused(root, LayeredOptions.CROSSING_MINIMIZATION_FORCE_NODE_MODEL_ORDER, true);
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Supported configurations

    /** The top-level graph of an {@link HierarchyHandling#INCLUDE_CHILDREN} hierarchy has its model order data. */
    @Test
    public void testModelOrderCycleBreakingOnRootLaysOut() {
        ElkNode root = root(HierarchyHandling.INCLUDE_CHILDREN);
        root.setProperty(LayeredOptions.CYCLE_BREAKING_STRATEGY, CycleBreakingStrategy.MODEL_ORDER);
        nested(root, 2);
        node(root, "a");
        node(root, "b");
        ElkNode a = root.getChildren().get(1);
        ElkNode b = root.getChildren().get(2);
        ElkGraphUtil.createSimpleEdge(a, b);
        ElkGraphUtil.createSimpleEdge(b, a);

        layout(root);

        // The model order cycle breaker reverses b -> a, the edge against the model order
        assertTrue(a.getX() < b.getX());
    }

    /** A nested graph laid out on its own is imported flat and has its model order data. */
    @Test
    public void testModelOrderCycleBreakingWithSeparateChildrenLaysOut() {
        ElkNode root = root(HierarchyHandling.SEPARATE_CHILDREN);
        ElkNode nested = nested(root, 2);
        connect(nested, true);
        nested.setProperty(LayeredOptions.CYCLE_BREAKING_STRATEGY, CycleBreakingStrategy.MODEL_ORDER);

        layout(root);

        // The model order cycle breaker reverses n1 -> n0, the edge against the model order
        assertTrue(nested.getChildren().get(0).getX() < nested.getChildren().get(1).getX());
    }

    /** Without the enforced group order, the node order cycle breakers read only the node model order. */
    @Test
    public void testDepthFirstNodeOrderWithoutEnforcedGroupOrderLaysOut() {
        ElkNode root = root(HierarchyHandling.INCLUDE_CHILDREN);
        ElkNode nested = nested(root, 2);
        connect(nested, false);
        nested.setProperty(LayeredOptions.CYCLE_BREAKING_STRATEGY, CycleBreakingStrategy.DFS_NODE_ORDER);

        layout(root);

        assertTrue(nested.getChildren().get(0).getX() < nested.getChildren().get(1).getX());
    }

}
