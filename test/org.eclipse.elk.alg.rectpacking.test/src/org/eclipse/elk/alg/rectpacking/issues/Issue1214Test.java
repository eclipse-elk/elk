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
package org.eclipse.elk.alg.rectpacking.issues;

import static org.junit.Assert.assertEquals;

import org.eclipse.elk.alg.rectpacking.RectPackingLayoutProvider;
import org.eclipse.elk.alg.rectpacking.options.RectPackingOptions;
import org.eclipse.elk.alg.test.PlainJavaInitialization;
import org.eclipse.elk.core.math.ElkPadding;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.util.BasicProgressMonitor;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.util.ElkGraphUtil;
import org.junit.BeforeClass;
import org.junit.Test;

/**
 * Test whether the rectpacking algorithm can place nodes without a size if there is no spacing between the nodes. All
 * options to place such a node have the same drawing, which used to cause a {@link NullPointerException}.
 */
public class Issue1214Test {

    @BeforeClass
    public static void init() {
        PlainJavaInitialization.initializePlainJavaLayout();
    }

    private ElkNode createParent() {
        ElkNode parent = ElkGraphUtil.createGraph();
        parent.setProperty(CoreOptions.ALGORITHM, RectPackingOptions.ALGORITHM_ID);
        parent.setProperty(CoreOptions.SPACING_NODE_NODE, 0.0);
        parent.setProperty(CoreOptions.PADDING, new ElkPadding(10.0));
        return parent;
    }

    private void layout(final ElkNode parent) {
        new RectPackingLayoutProvider().layout(parent, new BasicProgressMonitor());
    }

    @Test
    public void testTwoNodesWithoutSize() {
        ElkNode parent = createParent();
        ElkGraphUtil.createNode(parent).setDimensions(0, 0);
        ElkGraphUtil.createNode(parent).setDimensions(0, 0);

        layout(parent);

        assertEquals(20.0, parent.getWidth(), 0);
        assertEquals(20.0, parent.getHeight(), 0);
    }

    @Test
    public void testThreeNodesWithoutSize() {
        ElkNode parent = createParent();
        ElkGraphUtil.createNode(parent).setDimensions(0, 0);
        ElkGraphUtil.createNode(parent).setDimensions(0, 0);
        ElkGraphUtil.createNode(parent).setDimensions(0, 0);

        layout(parent);

        assertEquals(20.0, parent.getWidth(), 0);
        assertEquals(20.0, parent.getHeight(), 0);
    }

    @Test
    public void testNodeWithoutSizeAfterNodeWithSize() {
        ElkNode parent = createParent();
        ElkNode regular = ElkGraphUtil.createNode(parent);
        regular.setDimensions(10, 10);
        ElkGraphUtil.createNode(parent).setDimensions(0, 0);

        layout(parent);

        assertEquals(10.0, regular.getX(), 0);
        assertEquals(10.0, regular.getY(), 0);
        assertEquals(30.0, parent.getWidth(), 0);
    }

    @Test
    public void testNodeWithoutSizeBeforeNodeWithSize() {
        ElkNode parent = createParent();
        ElkGraphUtil.createNode(parent).setDimensions(0, 0);
        ElkNode regular = ElkGraphUtil.createNode(parent);
        regular.setDimensions(10, 10);

        layout(parent);

        assertEquals(30.0, parent.getWidth(), 0);
    }
}
