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

import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertTrue;

import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.test.PlainJavaInitialization;
import org.eclipse.elk.core.RecursiveGraphLayoutEngine;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.options.Direction;
import org.eclipse.elk.core.util.BasicProgressMonitor;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.util.ElkGraphUtil;
import org.junit.BeforeClass;
import org.junit.Test;

/**
 * Tests the in-layer constraints of the {@link InteractiveLayeredGraphVisitor} on the graph
 * {@code n0 -> n1, n0 -> n2, n1 -> n3, n2 -> n3, n3 -> n4, n3 -> n5}, laid out, visited and laid out again as an
 * interactive editor does.
 */
public class InteractiveLayeredGraphVisitorTest {

    private static final Direction[] DIRECTIONS = { Direction.RIGHT, Direction.LEFT, Direction.DOWN, Direction.UP };

    @BeforeClass
    public static void init() {
        PlainJavaInitialization.initializePlainJavaLayout();
    }

    private static ElkNode[] graph(final Direction direction) {
        ElkNode parent = ElkGraphUtil.createGraph();
        parent.setProperty(CoreOptions.ALGORITHM, LayeredOptions.ALGORITHM_ID);
        parent.setProperty(CoreOptions.DIRECTION, direction);
        parent.setProperty(LayeredOptions.GENERATE_POSITION_AND_LAYER_IDS, true);
        ElkNode[] nodes = new ElkNode[6];
        for (int i = 0; i < nodes.length; i++) {
            nodes[i] = ElkGraphUtil.createNode(parent);
            nodes[i].setIdentifier("n" + i);
            nodes[i].setDimensions(30, 30);
        }
        int[][] edges = { { 0, 1 }, { 0, 2 }, { 1, 3 }, { 2, 3 }, { 3, 4 }, { 3, 5 } };
        for (int[] edge : edges) {
            ElkGraphUtil.createSimpleEdge(nodes[edge[0]], nodes[edge[1]]);
        }
        return nodes;
    }

    private static void layoutVisitLayout(final ElkNode[] nodes) {
        ElkNode parent = nodes[0].getParent();
        new RecursiveGraphLayoutEngine().layout(parent, new BasicProgressMonitor());
        new InteractiveLayeredGraphVisitor().visit(parent);
        new RecursiveGraphLayoutEngine().layout(parent, new BasicProgressMonitor());
    }

    private static boolean horizontal(final Direction direction) {
        return direction == Direction.RIGHT || direction == Direction.LEFT;
    }

    /** The coordinate along the layout direction, which tells the layers apart. */
    private static double layerCoordinate(final ElkNode node, final Direction direction) {
        return horizontal(direction) ? node.getX() : node.getY();
    }

    /** The coordinate within a layer. */
    private static double inLayerCoordinate(final ElkNode node, final Direction direction) {
        return horizontal(direction) ? node.getY() : node.getX();
    }

    @Test
    public void testSuccessorConstraintToConnectedNodeIsIgnored() {
        // n3 -> n4 is an edge, so n3 and n4 can never share a layer. The visitor used to throw an
        // IndexOutOfBoundsException because n3 was left out of the order of its layer.
        for (Direction direction : DIRECTIONS) {
            ElkNode[] n = graph(direction);
            n[3].setProperty(LayeredOptions.CROSSING_MINIMIZATION_IN_LAYER_SUCC_OF, "n4");
            layoutVisitLayout(n);
            assertNotEquals(direction.toString(), layerCoordinate(n[3], direction), layerCoordinate(n[4], direction),
                    0.0);
        }
    }

    @Test
    public void testPredecessorAndSuccessorConstraintsToConnectedNodesAreIgnored() {
        // Both targets are connected to n3 by edges. The visitor used to throw a NullPointerException.
        for (Direction direction : DIRECTIONS) {
            ElkNode[] n = graph(direction);
            n[3].setProperty(LayeredOptions.CROSSING_MINIMIZATION_IN_LAYER_PRED_OF, "n2");
            n[3].setProperty(LayeredOptions.CROSSING_MINIMIZATION_IN_LAYER_SUCC_OF, "n4");
            layoutVisitLayout(n);
            assertNotEquals(direction.toString(), layerCoordinate(n[3], direction), layerCoordinate(n[2], direction),
                    0.0);
            assertNotEquals(direction.toString(), layerCoordinate(n[3], direction), layerCoordinate(n[4], direction),
                    0.0);
        }
    }

    @Test
    public void testSuccessorConstraintInLayerIsKept() {
        // Without the constraint n4 comes first in its layer; with it, n4 follows n5.
        for (Direction direction : DIRECTIONS) {
            ElkNode[] n = graph(direction);
            n[4].setProperty(LayeredOptions.CROSSING_MINIMIZATION_IN_LAYER_SUCC_OF, "n5");
            layoutVisitLayout(n);
            assertTrue(direction.toString(), inLayerCoordinate(n[5], direction) < inLayerCoordinate(n[4], direction));
        }
    }
}
