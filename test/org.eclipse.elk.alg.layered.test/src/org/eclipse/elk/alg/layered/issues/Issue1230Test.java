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
package org.eclipse.elk.alg.layered.issues;

import static org.junit.Assert.assertTrue;

import java.util.Arrays;
import java.util.Collection;

import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.core.RecursiveGraphLayoutEngine;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.options.HierarchyHandling;
import org.eclipse.elk.core.options.PortConstraints;
import org.eclipse.elk.core.options.PortSide;
import org.eclipse.elk.core.util.NullElkProgressMonitor;
import org.eclipse.elk.graph.ElkEdge;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.ElkPort;
import org.eclipse.elk.graph.util.ElkGraphUtil;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameters;

/**
 * Characterizes port order with inside self loops, including the current bug in issue 1230.
 */
@RunWith(Parameterized.class)
public class Issue1230Test {

    @Parameters(name = "{0}, inside loop = {1}")
    public static Collection<Object[]> configurations() {
        return Arrays.asList(new Object[][] {
            {HierarchyHandling.INCLUDE_CHILDREN, true},
            {HierarchyHandling.INCLUDE_CHILDREN, false},
            {HierarchyHandling.SEPARATE_CHILDREN, true},
            {HierarchyHandling.SEPARATE_CHILDREN, false}
        });
    }

    private final HierarchyHandling hierarchyHandling;
    private final boolean insideLoop;

    public Issue1230Test(final HierarchyHandling hierarchyHandling, final boolean insideLoop) {
        this.hierarchyHandling = hierarchyHandling;
        this.insideLoop = insideLoop;
    }

    @Test
    public void testFixedPortOrder() {
        testFixedPortOrder(true);
    }

    @Test
    public void testFixedPortOrderWithoutIndices() {
        testFixedPortOrder(false);
    }

    private void testFixedPortOrder(final boolean explicitIndices) {
        ElkNode graph = ElkGraphUtil.createGraph();
        graph.setProperty(CoreOptions.ALGORITHM, LayeredOptions.ALGORITHM_ID);
        graph.setProperty(LayeredOptions.HIERARCHY_HANDLING, hierarchyHandling);
        ElkNode closure = ElkGraphUtil.createNode(graph);
        closure.setDimensions(100, 100);
        closure.setProperty(CoreOptions.ALGORITHM, LayeredOptions.ALGORITHM_ID);
        closure.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_ORDER);
        closure.setProperty(LayeredOptions.INSIDE_SELF_LOOPS_ACTIVATE, true);

        // Indices run clockwise: WEST ports therefore appear in reverse index order from top to bottom.
        ElkPort w3 = port(closure, PortSide.WEST, 0);
        ElkPort w2 = port(closure, PortSide.WEST, 1);
        ElkPort w1 = port(closure, PortSide.WEST, 2);
        ElkPort e1 = port(closure, PortSide.EAST, 3);
        ElkPort e2 = port(closure, PortSide.EAST, 4);
        ElkPort e3 = port(closure, PortSide.EAST, 5);
        if (!explicitIndices) {
            // Without indices, the initial positions define the order on each side.
            w1.setLocation(-4, 20);
            w2.setLocation(-4, 40);
            w3.setLocation(-4, 60);
            e1.setLocation(100, 20);
            e2.setLocation(100, 40);
            e3.setLocation(100, 60);
            for (ElkPort port : closure.getPorts()) {
                port.setProperty(LayeredOptions.PORT_INDEX, null);
            }
        }
        if (insideLoop) {
            ElkEdge edge = ElkGraphUtil.createSimpleEdge(w1, e3);
            edge.setProperty(LayeredOptions.INSIDE_SELF_LOOPS_YO, true);
        }

        new RecursiveGraphLayoutEngine().layout(graph, new NullElkProgressMonitor());

        if (hierarchyHandling == HierarchyHandling.INCLUDE_CHILDREN && insideLoop) {
            // Issue 1230: record the current incorrect order until the production fix is applied.
            // The fix must replace these expectations with w1, w2, w3, as in the other cases.
            if (explicitIndices) {
                assertTrue("Current WEST order with indices is w2, w1, w3",
                        w2.getY() < w1.getY() && w1.getY() < w3.getY());
            } else {
                assertTrue("Current WEST order without indices is w1, w3, w2",
                        w1.getY() < w3.getY() && w3.getY() < w2.getY());
            }
        } else {
            assertTrue("WEST ports must retain their order", w1.getY() < w2.getY() && w2.getY() < w3.getY());
        }
        assertTrue("EAST ports must retain their order", e1.getY() < e2.getY() && e2.getY() < e3.getY());
        for (ElkPort port : Arrays.asList(w1, w2, w3)) {
            assertTrue("WEST ports must retain their side", port.getProperty(LayeredOptions.PORT_SIDE) == PortSide.WEST);
        }
        for (ElkPort port : Arrays.asList(e1, e2, e3)) {
            assertTrue("EAST ports must retain their side", port.getProperty(LayeredOptions.PORT_SIDE) == PortSide.EAST);
        }
    }

    private ElkPort port(final ElkNode node, final PortSide side, final int index) {
        ElkPort port = ElkGraphUtil.createPort(node);
        port.setDimensions(4, 4);
        port.setProperty(LayeredOptions.PORT_SIDE, side);
        port.setProperty(LayeredOptions.PORT_INDEX, index);
        return port;
    }
}
