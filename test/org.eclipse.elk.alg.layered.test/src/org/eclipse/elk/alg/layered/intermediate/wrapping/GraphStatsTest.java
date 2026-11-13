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
package org.eclipse.elk.alg.layered.intermediate.wrapping;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;

import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.options.CuttingStrategy;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.layered.options.WrappingStrategy;
import org.eclipse.elk.alg.test.PlainJavaInitialization;
import org.eclipse.elk.core.RecursiveGraphLayoutEngine;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.util.BasicProgressMonitor;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.util.ElkGraphUtil;
import org.junit.BeforeClass;
import org.junit.Test;

import com.google.common.collect.Lists;

/**
 * Tests which cut indices {@link GraphStats} allows, with and without
 * {@link LayeredOptions#WRAPPING_VALIDIFY_FORBIDDEN_INDICES}, and that a wrapped layout keeps the cuts that are not
 * forbidden.
 */
public class GraphStatsTest {

    @BeforeClass
    public static void init() {
        PlainJavaInitialization.initializePlainJavaLayout();
    }

    /** A path {@code n0 -> n1 -> n2 -> n3 -> n4} with one node per layer: every cut but the first is allowed. */
    private static LGraph pathGraph() {
        LGraph graph = new LGraph();
        LNode previous = node(graph, layer(graph));
        for (int i = 1; i < 5; i++) {
            LNode next = node(graph, layer(graph));
            edge(previous, next);
            previous = next;
        }
        return graph;
    }

    /**
     * {@code a -> b, a -> c, b -> d, c -> d, d -> e, e -> f} in the layers {@code [a] [b, c] [d] [e] [f]}: the edges
     * into layers 1 and 2 have different sources or targets, so only cuts before layers 3 and 4 are allowed.
     */
    private static LGraph diamondGraph() {
        LGraph graph = new LGraph();
        LNode a = node(graph, layer(graph));
        Layer layer1 = layer(graph);
        LNode b = node(graph, layer1);
        LNode c = node(graph, layer1);
        LNode d = node(graph, layer(graph));
        LNode e = node(graph, layer(graph));
        LNode f = node(graph, layer(graph));
        edge(a, b);
        edge(a, c);
        edge(b, d);
        edge(c, d);
        edge(d, e);
        edge(e, f);
        return graph;
    }

    /**
     * Lays out a path of six nodes with {@code SINGLE_EDGE} wrapping and a manual cut before the fourth node, with the
     * given forbidden indices if they are not {@code null}.
     */
    private static ElkNode layoutWrappedPath(final List<Integer> forbiddenIndices) {
        ElkNode graph = ElkGraphUtil.createGraph();
        graph.setProperty(CoreOptions.ALGORITHM, LayeredOptions.ALGORITHM_ID);
        graph.setProperty(LayeredOptions.WRAPPING_STRATEGY, WrappingStrategy.SINGLE_EDGE);
        graph.setProperty(LayeredOptions.WRAPPING_CUTTING_STRATEGY, CuttingStrategy.MANUAL);
        graph.setProperty(LayeredOptions.WRAPPING_CUTTING_CUTS, Lists.newArrayList(3));
        if (forbiddenIndices != null) {
            graph.setProperty(LayeredOptions.WRAPPING_VALIDIFY_FORBIDDEN_INDICES, forbiddenIndices);
        }

        ElkNode previous = null;
        for (int i = 0; i < 6; i++) {
            ElkNode node = ElkGraphUtil.createNode(graph);
            node.setDimensions(30, 30);
            if (previous != null) {
                ElkGraphUtil.createSimpleEdge(previous, node);
            }
            previous = node;
        }

        new RecursiveGraphLayoutEngine().layout(graph, new BasicProgressMonitor());
        return graph;
    }

    /** Whether the fourth node of a path laid out by {@link #layoutWrappedPath(List)} starts a new chunk. */
    private static boolean isWrappedBeforeFourthNode(final ElkNode graph) {
        ElkNode third = graph.getChildren().get(2);
        ElkNode fourth = graph.getChildren().get(3);
        return fourth.getX() < third.getX() && fourth.getY() > third.getY();
    }

    private static Layer layer(final LGraph graph) {
        Layer layer = new Layer(graph);
        graph.getLayers().add(layer);
        return layer;
    }

    private static LNode node(final LGraph graph, final Layer layer) {
        LNode node = new LNode(graph);
        node.setLayer(layer);
        return node;
    }

    private static void edge(final LNode source, final LNode target) {
        LPort sourcePort = new LPort();
        sourcePort.setNode(source);
        LPort targetPort = new LPort();
        targetPort.setNode(target);
        LEdge edge = new LEdge();
        edge.setSource(sourcePort);
        edge.setTarget(targetPort);
    }

    @Test
    public void testCutsAllowedWithoutForbiddenIndices() {
        assertArrayEquals(new boolean[] { false, true, true, true, true },
                new GraphStats(pathGraph()).getCutsAllowed());
        assertArrayEquals(new boolean[] { false, false, false, true, true },
                new GraphStats(diamondGraph()).getCutsAllowed());
    }

    @Test
    public void testForbiddenIndicesForbidOnlyTheListedCuts() {
        LGraph graph = pathGraph();
        graph.setProperty(LayeredOptions.WRAPPING_VALIDIFY_FORBIDDEN_INDICES, Lists.newArrayList(2, 4));
        assertArrayEquals(new boolean[] { false, true, false, true, false }, new GraphStats(graph).getCutsAllowed());
    }

    @Test
    public void testForbiddenIndicesKeepCutsForbiddenByTheGraph() {
        LGraph graph = diamondGraph();
        graph.setProperty(LayeredOptions.WRAPPING_VALIDIFY_FORBIDDEN_INDICES, Lists.newArrayList(3));
        assertArrayEquals(new boolean[] { false, false, false, false, true }, new GraphStats(graph).getCutsAllowed());
    }

    @Test
    public void testForbiddenIndicesOutsideTheLayeringAreIgnored() {
        LGraph graph = pathGraph();
        graph.setProperty(LayeredOptions.WRAPPING_VALIDIFY_FORBIDDEN_INDICES, Lists.newArrayList(-1, 0, 5));
        assertArrayEquals(new boolean[] { false, true, true, true, true }, new GraphStats(graph).getCutsAllowed());
    }

    @Test
    public void testManualCutWithoutForbiddenIndices() {
        assertTrue(isWrappedBeforeFourthNode(layoutWrappedPath(null)));
    }

    @Test
    public void testManualCutAllowedNextToForbiddenIndex() {
        assertTrue(isWrappedBeforeFourthNode(layoutWrappedPath(Lists.newArrayList(1))));
    }

    @Test
    public void testManualCutAtForbiddenIndexIsMoved() {
        // The greedy validify strategy moves the cut to the next allowed index
        ElkNode graph = layoutWrappedPath(Lists.newArrayList(3));
        assertFalse(isWrappedBeforeFourthNode(graph));
        ElkNode fourth = graph.getChildren().get(3);
        ElkNode fifth = graph.getChildren().get(4);
        assertTrue(fifth.getX() < fourth.getX() && fifth.getY() > fourth.getY());
    }
}
