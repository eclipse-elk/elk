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
package org.eclipse.elk.alg.rectpacking;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.List;

import org.eclipse.elk.alg.rectpacking.options.RectPackingOptions;
import org.eclipse.elk.alg.test.PlainJavaInitialization;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.util.BasicProgressMonitor;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.util.ElkGraphUtil;
import org.junit.BeforeClass;
import org.junit.Test;

/**
 * Tests the order the interactive node reorderer gives to nodes that share a desired position. Such nodes keep their
 * model order, also in lists long enough that a comparator that is not a total order makes the sort throw or return an
 * inconsistent order.
 */
public class InteractiveNodeReordererTest {

    /** Desired positions of 32 nodes for which a sort with an inconsistent tie comparator threw. */
    private static final int[] THROWING_POSITIONS = { 1, 0, 0, 0, 1, 0, 2, 0, 0, 1, 2, 2, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0,
            2, 0, 1, 0, 1, 1, 1, 0, 2, 1 };

    /** Desired positions of 32 nodes for which a sort with an inconsistent tie comparator split the ties of 2. */
    private static final int[] SPLIT_TIES_POSITIONS = { 2, 2, 0, 2, 2, 0, 2, 2, 2, 2, 0, 2, 2, 2, 2, 0, 0, 2, 2, 0, 0, 2,
            2, 0, 2, 2, 2, 0, 0, 0, 2, 2 };

    @BeforeClass
    public static void init() {
        PlainJavaInitialization.initializePlainJavaLayout();
    }

    private ElkNode createParent() {
        ElkNode parent = ElkGraphUtil.createGraph();
        parent.setProperty(CoreOptions.ALGORITHM, RectPackingOptions.ALGORITHM_ID);
        parent.setProperty(RectPackingOptions.INTERACTIVE, true);
        return parent;
    }

    private ElkNode createNode(final ElkNode parent) {
        ElkNode node = ElkGraphUtil.createNode(parent);
        node.setDimensions(20, 20);
        return node;
    }

    private ElkNode createNode(final ElkNode parent, final int desiredPosition) {
        ElkNode node = createNode(parent);
        node.setProperty(RectPackingOptions.DESIRED_POSITION, desiredPosition);
        return node;
    }

    private List<ElkNode> createNodes(final ElkNode parent, final int[] desiredPositions) {
        List<ElkNode> nodes = new ArrayList<>();
        for (int desiredPosition : desiredPositions) {
            nodes.add(createNode(parent, desiredPosition));
        }
        return nodes;
    }

    private void layout(final ElkNode parent) {
        new RectPackingLayoutProvider().layout(parent, new BasicProgressMonitor());
    }

    /**
     * Asserts that of any two nodes with the same desired position, the one earlier in the given model order gets the
     * lower current position.
     */
    private void assertTiesKeepModelOrder(final List<ElkNode> nodesInModelOrder) {
        for (int i = 0; i < nodesInModelOrder.size(); i++) {
            ElkNode earlier = nodesInModelOrder.get(i);
            for (int j = i + 1; j < nodesInModelOrder.size(); j++) {
                ElkNode later = nodesInModelOrder.get(j);
                if (earlier.getProperty(RectPackingOptions.DESIRED_POSITION)
                        .equals(later.getProperty(RectPackingOptions.DESIRED_POSITION))) {
                    int earlierPosition = earlier.getProperty(RectPackingOptions.CURRENT_POSITION);
                    int laterPosition = later.getProperty(RectPackingOptions.CURRENT_POSITION);
                    assertTrue("Node " + i + " (current position " + earlierPosition + ") should precede node " + j
                            + " (current position " + laterPosition + "), which has the same desired position",
                            earlierPosition < laterPosition);
                }
            }
        }
    }

    @Test
    public void testTiesInShortListKeepModelOrder() {
        ElkNode parent = createParent();
        ElkNode x = createNode(parent);
        ElkNode a = createNode(parent, 1);
        ElkNode y = createNode(parent);
        ElkNode b = createNode(parent, 1);
        ElkNode z = createNode(parent);
        ElkNode c = createNode(parent, 1);

        layout(parent);

        assertEquals(0, (int) x.getProperty(RectPackingOptions.CURRENT_POSITION));
        assertEquals(1, (int) a.getProperty(RectPackingOptions.CURRENT_POSITION));
        assertEquals(2, (int) b.getProperty(RectPackingOptions.CURRENT_POSITION));
        assertEquals(3, (int) c.getProperty(RectPackingOptions.CURRENT_POSITION));
        assertEquals(4, (int) y.getProperty(RectPackingOptions.CURRENT_POSITION));
        assertEquals(5, (int) z.getProperty(RectPackingOptions.CURRENT_POSITION));
    }

    @Test
    public void testManyTiesDoNotBreakTheSort() {
        ElkNode parent = createParent();
        List<ElkNode> nodes = createNodes(parent, THROWING_POSITIONS);

        layout(parent);

        assertTiesKeepModelOrder(nodes);
    }

    @Test
    public void testTiesInLongListKeepModelOrder() {
        ElkNode parent = createParent();
        List<ElkNode> nodes = createNodes(parent, SPLIT_TIES_POSITIONS);

        layout(parent);

        assertTiesKeepModelOrder(nodes);
    }
}
