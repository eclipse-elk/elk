/*******************************************************************************
 * Copyright (c) 2020 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.rectpacking.issues;

import static org.junit.Assert.assertEquals;

import org.eclipse.elk.alg.rectpacking.RectPackingLayoutProvider;
import org.eclipse.elk.alg.rectpacking.options.RectPackingOptions;
import org.eclipse.elk.alg.rectpacking.p2packing.PackingStrategy;
import org.eclipse.elk.alg.test.PlainJavaInitialization;
import org.eclipse.elk.core.math.ElkPadding;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.util.BasicProgressMonitor;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.util.ElkGraphUtil;
import org.junit.BeforeClass;
import org.junit.Test;

/**
 * Test whether the rectpacking algorithm reduces the approximated width of the drawing if it does not need the width.
 *
 */
public class Issue583Test {

    @BeforeClass
    public static void init() {
        PlainJavaInitialization.initializePlainJavaLayout();
    }

    @Test
    public void testDecreaseApproximatedWidth() {
        // Create a graph that has a higher approximated width than actual width

        ElkNode parent = ElkGraphUtil.createGraph();
        ElkNode n1 = ElkGraphUtil.createNode(parent);
        ElkNode n2 = ElkGraphUtil.createNode(parent);
        ElkNode n3 = ElkGraphUtil.createNode(parent);
        ElkNode n4 = ElkGraphUtil.createNode(parent);
        parent.setProperty(CoreOptions.ALGORITHM, RectPackingOptions.ALGORITHM_ID);
        parent.setProperty(CoreOptions.SPACING_NODE_NODE, 0.0);
        parent.setProperty(CoreOptions.PADDING, new ElkPadding(0.0));
        n1.setDimensions(30, 30);
        n2.setDimensions(10, 10);
        n3.setDimensions(30, 30);
        // n4 is added to disable last row optimization
        n4.setDimensions(40, 10);
        
        // ___
        //|   | O
        //|___| ___
        //     |   |
        //     |___|
        // ______
        //|______|
        // Expected drawing with width approximation only
        parent.setProperty(RectPackingOptions.PACKING_STRATEGY, PackingStrategy.NONE);
        RectPackingLayoutProvider layoutProvider = new RectPackingLayoutProvider();
        layoutProvider.layout(parent, new BasicProgressMonitor());
        assertEquals("", 60.0, parent.getWidth(), 1);
        
        // ___
        //|   | O  | Approximated maximum width
        //|___|    |
        // ___     |
        //|   |    |
        //|___|    |
        // ______  |
        //|______| |
        // Expected drawing with full algorithm
        parent.setProperty(RectPackingOptions.PACKING_STRATEGY, PackingStrategy.COMPACTION);
        layoutProvider = new RectPackingLayoutProvider();
        layoutProvider.layout(parent, new BasicProgressMonitor());
        assertEquals("", 40.0, parent.getWidth(), 1);
    }

}
