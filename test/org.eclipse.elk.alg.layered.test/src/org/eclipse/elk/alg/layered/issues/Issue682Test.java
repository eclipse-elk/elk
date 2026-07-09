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

import java.util.Arrays;
import java.util.Collection;
import java.util.EnumSet;

import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.test.PlainJavaInitialization;
import org.eclipse.elk.core.RecursiveGraphLayoutEngine;
import org.eclipse.elk.core.math.ElkPadding;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.options.Direction;
import org.eclipse.elk.core.options.EdgeRouting;
import org.eclipse.elk.core.options.NodeLabelPlacement;
import org.eclipse.elk.core.options.SizeConstraint;
import org.eclipse.elk.core.util.NullElkProgressMonitor;
import org.eclipse.elk.core.util.Triple;
import org.eclipse.elk.graph.ElkLabel;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.util.ElkGraphUtil;
import org.junit.BeforeClass;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameter;
import org.junit.runners.Parameterized.Parameters;

/**
 * Test for issue 682.
 */
@RunWith(Parameterized.class)
public class Issue682Test {

    @BeforeClass
    public static void init() {
        PlainJavaInitialization.initializePlainJavaLayout();
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Tests

    @Parameters(name = "{0}")
    public static Collection<Object[]> data() {
        return Arrays
                .asList(new Object[][] { { Direction.RIGHT }, { Direction.DOWN }, { Direction.UP }, { Direction.LEFT },

                });
    }

    @Parameter(0)
    public Direction layoutDirection;

    public Triple<ElkNode, ElkNode, ElkLabel> testGraph() {
        ElkNode graph = ElkGraphUtil.createGraph();
        graph.setProperty(CoreOptions.ALGORITHM, LayeredOptions.ALGORITHM_ID);
        graph.setProperty(LayeredOptions.EDGE_ROUTING, EdgeRouting.ORTHOGONAL);
        graph.setProperty(LayeredOptions.DIRECTION, layoutDirection);

        graph.setProperty(LayeredOptions.NODE_LABELS_PADDING, new ElkPadding(21, 32, 43, 54));

        ElkNode parent = ElkGraphUtil.createNode(graph);
        parent.setProperty(LayeredOptions.NODE_SIZE_CONSTRAINTS, EnumSet.of(SizeConstraint.NODE_LABELS));
        parent.setProperty(LayeredOptions.NODE_LABELS_PLACEMENT, NodeLabelPlacement.insideTopCenter());
        ElkLabel label = ElkGraphUtil.createLabel("foobar", parent);
        label.setWidth(23); // Arbitrary
        label.setHeight(22);

        return new Triple<>(graph, parent, label);
    }

    private static final double EPSILON = 1e-5;

    @Test
    public void test() {
        Triple<ElkNode, ElkNode, ElkLabel> uut = testGraph();

        new RecursiveGraphLayoutEngine().layout(uut.getFirst(), new NullElkProgressMonitor());

        // Label position
        ElkLabel label = uut.getThird();
        assertThat(label.getX(), closeTo(54.0, EPSILON));
        assertThat(label.getY(), closeTo(21.0, EPSILON));
        // Node dimension
        assertThat(uut.getSecond().getWidth(), closeTo(label.getX() + label.getWidth() + 32, EPSILON));

        // TODO the following line fails, it looks like the label's height is applied twice.
        // assertThat(uut.getSecond().getHeight(), closeTo(label.getY() + label.getHeight() + 43, EPSILON));
    }

}
