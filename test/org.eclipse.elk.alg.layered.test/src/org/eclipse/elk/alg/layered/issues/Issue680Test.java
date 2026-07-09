/*******************************************************************************
 * Copyright (c) 2020 Kiel University and others.
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

import static org.hamcrest.Matchers.closeTo;
import static org.junit.Assert.assertThat;

import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.test.PlainJavaInitialization;
import org.eclipse.elk.core.RecursiveGraphLayoutEngine;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.options.Direction;
import org.eclipse.elk.core.options.EdgeRouting;
import org.eclipse.elk.core.util.NullElkProgressMonitor;
import org.eclipse.elk.core.util.Triple;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.ElkPort;
import org.eclipse.elk.graph.util.ElkGraphUtil;
import org.junit.BeforeClass;
import org.junit.Test;

/**
 * Test for issue 680.
 */
public class Issue680Test {

    @BeforeClass
    public static void init() {
        PlainJavaInitialization.initializePlainJavaLayout();
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Tests

    public Triple<ElkNode, ElkNode, ElkNode> testGraph() {
        ElkNode graph = ElkGraphUtil.createGraph();
        graph.setProperty(CoreOptions.ALGORITHM, LayeredOptions.ALGORITHM_ID);
        graph.setProperty(LayeredOptions.EDGE_ROUTING, EdgeRouting.ORTHOGONAL);
        graph.setProperty(LayeredOptions.DIRECTION, Direction.DOWN);

        ElkNode parent = ElkGraphUtil.createNode(graph);
        parent.setProperty(CoreOptions.ALGORITHM, LayeredOptions.ALGORITHM_ID);
        parent.setProperty(LayeredOptions.EDGE_ROUTING, EdgeRouting.ORTHOGONAL);
        parent.setProperty(LayeredOptions.DIRECTION, Direction.DOWN);
        ElkPort p1 = ElkGraphUtil.createPort(parent);
        p1.setWidth(15);
        p1.setHeight(165);
        p1.setProperty(LayeredOptions.PORT_BORDER_OFFSET, -20.0);
        ElkPort p2 = ElkGraphUtil.createPort(parent);
        p2.setWidth(15);
        p2.setHeight(166);
        p2.setProperty(LayeredOptions.PORT_BORDER_OFFSET, -22.0);

        ElkNode child = ElkGraphUtil.createNode(parent);
        child.setWidth(40.265625);
        child.setHeight(75.5);
        ElkPort childP1 = ElkGraphUtil.createPort(child);
        childP1.setWidth(15);
        childP1.setHeight(33);
        childP1.setProperty(LayeredOptions.PORT_BORDER_OFFSET, -8.0);
        ElkPort childP2 = ElkGraphUtil.createPort(child);
        childP2.setWidth(15);
        childP2.setHeight(34);
        childP2.setProperty(LayeredOptions.PORT_BORDER_OFFSET, -8.0);

        ElkGraphUtil.createSimpleEdge(p1, childP1);
        ElkGraphUtil.createSimpleEdge(childP2, p2);

        return new Triple<>(graph, parent, child);
    }

    private static final double EPSILON = 1e-5;

    @Test
    public void test() {
        Triple<ElkNode, ElkNode, ElkNode> nodes = testGraph();

        new RecursiveGraphLayoutEngine().layout(nodes.getFirst(), new NullElkProgressMonitor());

        assertThat(nodes.getSecond().getY(), closeTo(157.0, EPSILON));
        assertThat(nodes.getThird().getY(), closeTo(57.0, EPSILON));
    }

}
