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

import static org.junit.Assert.assertTrue;

import java.util.List;

import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.test.framework.LayoutTestRunner;
import org.eclipse.elk.alg.test.framework.annotations.Algorithm;
import org.eclipse.elk.alg.test.framework.annotations.DefaultConfiguration;
import org.eclipse.elk.alg.test.framework.annotations.GraphResourceProvider;
import org.eclipse.elk.alg.test.framework.io.AbstractResourcePath;
import org.eclipse.elk.alg.test.framework.io.ModelResourcePath;
import org.eclipse.elk.core.options.PortSide;
import org.eclipse.elk.graph.ElkLabel;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.ElkPort;
import org.eclipse.elk.graph.ElkShape;
import org.junit.Test;
import org.junit.runner.RunWith;

import com.google.common.collect.Lists;
import com.google.common.math.DoubleMath;

/**
 * Further test for issue 405. The difference compared to {@link Issue405Test} is that the models being tested here are
 * not completely symmetric when it comes to their port label placement. This is due to the way port label placement
 * works at the moment.
 * 
 * As a consequence this test has to be more permissive: the internal {@link #checkPositionEquals(ElkShape, ElkShape)}
 * method decides based on the {@link PortSide} which coordinate (x or y) to check for equality.
 */
@RunWith(LayoutTestRunner.class)
@Algorithm(LayeredOptions.ALGORITHM_ID)
@DefaultConfiguration
public class Issue405_2Test {

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Sources

    @GraphResourceProvider
    public List<AbstractResourcePath> testGraphs() {
        // Note that there are further test cases for #405. However, they cannot be tested in this test case
        // due to assumptions made here regarding the symmetry of computed port and port label positions.
        // They are checked in 'Issue405Test2'.
        return Lists.newArrayList(
                new ModelResourcePath("tickets/layered/405_differentPortLabelPositionsNSWEOutside.elkt"),
                new ModelResourcePath("tickets/layered/405_differentPortLabelPositionsTwoWestOutsideNextToPort.elkt"));
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Tests

    private static final double DOUBLE_EPSILON = 10e-4;
 
    private <T extends ElkShape> boolean checkPositionEquals(T first, T second) {
        final PortSide ps =
                first instanceof ElkLabel ? ((ElkLabel) first).getParent().getProperty(LayeredOptions.PORT_SIDE)
                        : first.getProperty(LayeredOptions.PORT_SIDE);
        if (PortSide.SIDES_EAST_WEST.contains(ps)) {
            return DoubleMath.fuzzyEquals(first.getX(), second.getX(), DOUBLE_EPSILON);
        } else {
            return DoubleMath.fuzzyEquals(first.getY(), second.getY(), DOUBLE_EPSILON);
        }
    }

    @Test
    public void testPortAndLabelPositionsEqual(final ElkNode graph) {

        // There are four (non-hierarchical) nodes within the diagram for which different layout directions are used.
        // Each node has a number of ports, each of which has a single label.
        // While we do not know in which "order" the ports are placed, we expect the collective of port positions 
        // (and port label positions) to match each other for all of the various layout directions. 
        final List<String> topNodes = Lists.newArrayList("Right", "Down", "Left", "Up");

        // Now check for each subsequent pair of these "layout-direction nodes" that their ports and labels match.
        topNodes.stream()
            .map(name -> Issue405Test.getRelevantPortsAndPortLabels(graph, name))
            .reduce((previous, current) -> {
                for (ElkPort port : previous.getFirst()) {
                    boolean foundExactlyOneMatch = current.getFirst().stream()
                            .map(other -> checkPositionEquals(port, other))
                            .reduce(false, (a, b) -> a ^ b);
                    assertTrue("Unexpected port position " + port, foundExactlyOneMatch);
                }
    
                for (ElkLabel label : previous.getSecond()) {
                    boolean foundExactlyOneMatch = current.getSecond().stream()
                            .map(other -> checkPositionEquals(label, other))
                            .reduce(false, (a, b) -> a ^ b);
                    assertTrue("Unexpected port label position " + label, foundExactlyOneMatch);
    
                }
                return current;
            });
    }

}
