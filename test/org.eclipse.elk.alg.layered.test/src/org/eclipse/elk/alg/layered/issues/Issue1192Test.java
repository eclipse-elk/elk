/*******************************************************************************
 * Copyright (c) 2026 Kiel University and others.
 *
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.alg.layered.issues;

import static org.junit.Assert.assertFalse;

import java.util.Arrays;
import java.util.Collection;

import org.eclipse.elk.alg.test.PlainJavaInitialization;
import org.eclipse.elk.core.RecursiveGraphLayoutEngine;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.options.Direction;
import org.eclipse.elk.core.options.HierarchyHandling;
import org.eclipse.elk.core.options.PortConstraints;
import org.eclipse.elk.core.options.PortSide;
import org.eclipse.elk.core.util.NullElkProgressMonitor;
import org.eclipse.elk.graph.ElkConnectableShape;
import org.eclipse.elk.graph.ElkEdge;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.ElkPort;
import org.eclipse.elk.graph.util.ElkGraphUtil;
import org.junit.BeforeClass;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameter;
import org.junit.runners.Parameterized.Parameters;

/**
 * Test for issue 1192: layered + INCLUDE_CHILDREN threw
 * "Expected N hierarchical ports, but found only 0" for a ported compound with cross-hierarchy
 * edges, because external port dummies of NORTH/SOUTH ports (which carry no layer constraint)
 * may be placed in the boundary layers that the hierarchical layer sweep assumed to consist of
 * sweep-facing-side dummies only.
 */
@RunWith(Parameterized.class)
public class Issue1192Test {

    @BeforeClass
    public static void init() {
        PlainJavaInitialization.initializePlainJavaLayout();
    }

    @Parameters(name = "{0}")
    public static Collection<Object[]> data() {
        return Arrays.asList(new Object[][] {
            { PortConstraints.FIXED_ORDER }, { PortConstraints.FIXED_SIDE } });
    }

    @Parameter
    public PortConstraints portConstraints;

    @Test
    public void test() {
        ElkNode root = ElkGraphUtil.createGraph();
        root.setProperty(CoreOptions.ALGORITHM, "org.eclipse.elk.layered");
        root.setProperty(CoreOptions.DIRECTION, Direction.RIGHT);
        root.setProperty(CoreOptions.HIERARCHY_HANDLING, HierarchyHandling.INCLUDE_CHILDREN);

        // Compound node `core` with BOTH ports AND child compounds.
        ElkNode core = node(root, 120, 60);
        core.setProperty(CoreOptions.PORT_CONSTRAINTS, portConstraints);
        ElkPort pNorth1 = port(core, PortSide.NORTH);
        ElkPort pNorth2 = port(core, PortSide.NORTH);
        port(core, PortSide.EAST); // referenced by NO edge — yet load-bearing for the crash
        ElkPort pEast2 = port(core, PortSide.EAST);

        ElkNode boxA = node(core, 120, 60);
        ElkNode leafA = node(boxA, 80, 40);
        ElkNode boxB = node(core, 120, 60);
        ElkNode leafB = node(boxB, 80, 40);
        ElkNode plain = node(core, 80, 40);

        // Edges fully inside `core`.
        edge(core, boxA, boxB);
        edge(core, leafA, plain);

        // Edges crossing core's boundary, contained at the root graph.
        edge(root, pNorth1, boxA);
        edge(root, pNorth2, boxB);
        edge(root, leafB, pEast2);

        new RecursiveGraphLayoutEngine().layout(root, new NullElkProgressMonitor());

        // All edges must actually be routed.
        for (ElkNode container : new ElkNode[] { root, core }) {
            for (ElkEdge e : container.getContainedEdges()) {
                assertFalse("edge " + e + " has no routing", e.getSections().isEmpty());
            }
        }
    }

    private static ElkNode node(ElkNode parent, double w, double h) {
        ElkNode n = ElkGraphUtil.createNode(parent);
        n.setDimensions(w, h);
        return n;
    }

    private static void edge(ElkNode container, ElkConnectableShape src, ElkConnectableShape tgt) {
        ElkEdge e = ElkGraphUtil.createEdge(container);
        e.getSources().add(src);
        e.getTargets().add(tgt);
    }

    private static ElkPort port(ElkNode parent, PortSide side) {
        ElkPort p = ElkGraphUtil.createPort(parent);
        p.setProperty(CoreOptions.PORT_SIDE, side);
        return p;
    }
}
