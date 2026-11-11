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
package org.eclipse.elk.alg.layered.p2layers;

import static org.junit.Assert.assertThrows;
import static org.junit.Assert.assertTrue;

import org.eclipse.elk.alg.layered.options.CycleBreakingStrategy;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.layered.options.LayeringStrategy;
import org.eclipse.elk.alg.test.PlainJavaInitialization;
import org.eclipse.elk.core.RecursiveGraphLayoutEngine;
import org.eclipse.elk.core.UnsupportedGraphException;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.util.BasicProgressMonitor;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.util.ElkGraphUtil;
import org.junit.BeforeClass;
import org.junit.Test;

/**
 * Tests that the BF and DF model order layer assigners refuse a layering whose edges do not all point to a later
 * layer. Such edges arise when an edge runs against the model order, which a cycle breaker other than the model
 * order one leaves in place. Before the check, they reached crossing minimization and could fail there with an
 * {@link IndexOutOfBoundsException}.
 */
public class ModelOrderLayeringEdgeDirectionTest {

    @BeforeClass
    public static void init() {
        PlainJavaInitialization.initializePlainJavaLayout();
    }

    /** Two nodes and one edge from the second to the first node in model order: {@code n0 <- n1}. */
    private static ElkNode backwardEdgeGraph(final LayeringStrategy layering, final CycleBreakingStrategy cycles) {
        ElkNode parent = graph(layering, cycles);
        ElkNode n0 = node(parent);
        ElkNode n1 = node(parent);
        ElkGraphUtil.createSimpleEdge(n1, n0);
        return parent;
    }

    /** {@code n0 -> n1 <- n2}: BF puts n2 into the layer of n1, DF starts a new strip with n2 in the first layer. */
    private static ElkNode joiningEdgeGraph(final LayeringStrategy layering) {
        ElkNode parent = graph(layering, CycleBreakingStrategy.GREEDY);
        ElkNode n0 = node(parent);
        ElkNode n1 = node(parent);
        ElkNode n2 = node(parent);
        ElkGraphUtil.createSimpleEdge(n0, n1);
        ElkGraphUtil.createSimpleEdge(n2, n1);
        return parent;
    }

    private static ElkNode graph(final LayeringStrategy layering, final CycleBreakingStrategy cycles) {
        ElkNode parent = ElkGraphUtil.createGraph();
        parent.setProperty(CoreOptions.ALGORITHM, LayeredOptions.ALGORITHM_ID);
        parent.setProperty(LayeredOptions.LAYERING_STRATEGY, layering);
        parent.setProperty(LayeredOptions.CYCLE_BREAKING_STRATEGY, cycles);
        return parent;
    }

    private static ElkNode node(final ElkNode parent) {
        ElkNode node = ElkGraphUtil.createNode(parent);
        node.setDimensions(30, 30);
        return node;
    }

    private static void layout(final ElkNode graph) {
        new RecursiveGraphLayoutEngine().layout(graph, new BasicProgressMonitor());
    }

    private static void assertRefused(final ElkNode graph, final String layerer) {
        UnsupportedGraphException e = assertThrows(UnsupportedGraphException.class, () -> layout(graph));
        assertTrue(e.getMessage(), e.getMessage().startsWith("The " + layerer + " model order layer assigner"));
        assertTrue(e.getMessage(), e.getMessage().contains(LayeredOptions.CYCLE_BREAKING_STRATEGY.getId()));
    }

    @Test
    public void testBreadthFirstRefusesBackwardEdge() {
        assertRefused(backwardEdgeGraph(LayeringStrategy.BF_MODEL_ORDER, CycleBreakingStrategy.GREEDY), "BF");
    }

    @Test
    public void testDepthFirstRefusesBackwardEdge() {
        assertRefused(backwardEdgeGraph(LayeringStrategy.DF_MODEL_ORDER, CycleBreakingStrategy.GREEDY), "DF");
    }

    @Test
    public void testBreadthFirstRefusesJoiningEdge() {
        assertRefused(joiningEdgeGraph(LayeringStrategy.BF_MODEL_ORDER), "BF");
    }

    @Test
    public void testDepthFirstKeepsValidLayeringOfJoiningEdge() {
        // n2 has no incoming edge and starts a new strip in the first layer, so its edge still points forward
        ElkNode graph = joiningEdgeGraph(LayeringStrategy.DF_MODEL_ORDER);
        layout(graph);
        ElkNode n1 = graph.getChildren().get(1);
        ElkNode n2 = graph.getChildren().get(2);
        assertTrue(n2.getX() < n1.getX());
    }

    @Test
    public void testModelOrderCycleBreakingReversesBackwardEdge() {
        for (LayeringStrategy layering : new LayeringStrategy[] {
                LayeringStrategy.BF_MODEL_ORDER, LayeringStrategy.DF_MODEL_ORDER }) {
            ElkNode graph = backwardEdgeGraph(layering, CycleBreakingStrategy.MODEL_ORDER);
            layout(graph);
            ElkNode n0 = graph.getChildren().get(0);
            ElkNode n1 = graph.getChildren().get(1);
            assertTrue(layering.toString(), n0.getX() < n1.getX());
        }
    }
}
