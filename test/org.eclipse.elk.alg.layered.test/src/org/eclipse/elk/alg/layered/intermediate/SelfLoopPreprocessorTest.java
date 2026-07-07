/*******************************************************************************
 * Copyright (c) 2019 Kiel University and others.
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
package org.eclipse.elk.alg.layered.intermediate;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.intermediate.loops.SelfLoopTestGraphCreator;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.core.options.PortConstraints;
import org.eclipse.elk.core.util.NullElkProgressMonitor;
import org.junit.Before;
import org.junit.Test;

/**
 * Tests whether the {@link SelfLoopPreprocessor} correctly hides ports and self loop edges.
 */
public class SelfLoopPreprocessorTest {

    private SelfLoopPreProcessor slPreprocessor;

    @Before
    public void createProcessor() {
        slPreprocessor = new SelfLoopPreProcessor();
    }

    @Test
    public void testFreePorts() {
        testPortsHidden(PortConstraints.FREE);
    }

    @Test
    public void testFixedSide() {
        testPortsHidden(PortConstraints.FIXED_SIDE);
    }

    @Test
    public void testFixedOrder() {
        testPortsUntouched(PortConstraints.FIXED_ORDER);
    }

    @Test
    public void testFixedRatio() {
        testPortsUntouched(PortConstraints.FIXED_RATIO);
    }

    @Test
    public void testFixedPos() {
        testPortsUntouched(PortConstraints.FIXED_POS);
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Actual Tests

    private void testPortsHidden(final PortConstraints portConstraints) {
        TestGraph testGraph = new TestGraph(portConstraints);
        slPreprocessor.process(testGraph.testGraph, new NullElkProgressMonitor());

        // Only the regular and unconnected ports should be left
        assertEquals(testGraph.regularPorts.size() + testGraph.unconnectedPorts.size(),
                testGraph.selfLoopNode.getPorts().size());

        assertTrue(testGraph.selfLoopNode.getPorts().containsAll(testGraph.regularPorts));
        assertTrue(testGraph.selfLoopNode.getPorts().containsAll(testGraph.unconnectedPorts));
        
        testHiddenEdges(testGraph);
    }

    private void testPortsUntouched(PortConstraints portConstraints) {
        TestGraph testGraph = new TestGraph(portConstraints);
        slPreprocessor.process(testGraph.testGraph, new NullElkProgressMonitor());

        // All ports should be left
        assertEquals(
                testGraph.regularPorts.size() + testGraph.onlySelfLoopPorts.size() + testGraph.unconnectedPorts.size(),
                testGraph.selfLoopNode.getPorts().size());

        testHiddenEdges(testGraph);
    }

    private void testHiddenEdges(TestGraph testGraph) {
        // Not a single self loop edge should be left
        for (LEdge lEdge : testGraph.selfLoopNode.getConnectedEdges()) {
            assertFalse(lEdge.isSelfLoop());
        }

        // A regular edge should be left
        assertTrue(testGraph.selfLoopNode.getConnectedEdges().iterator().hasNext());
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Test Graph Creation

    private static class TestGraph {
        private LGraph testGraph;
        private LNode selfLoopNode;
        private Set<LPort> regularPorts = new HashSet<>();
        private Set<LPort> onlySelfLoopPorts = new HashSet<>();
        private Set<LPort> unconnectedPorts = new HashSet<>();

        private TestGraph(final PortConstraints portConstraints) {
            testGraph = SelfLoopTestGraphCreator.basicGraphWithoutSelfLoops();

            selfLoopNode = testGraph.getLayerlessNodes().get(1);
            selfLoopNode.setProperty(LayeredOptions.PORT_CONSTRAINTS, portConstraints);
            selfLoopNode.getSize().set(50, 50);

            addPorts();
            addEdges();
        }

        private void addPorts() {
            SelfLoopTestGraphCreator.ports(selfLoopNode, 2, 2, 2, 1);

            // The ports need proper coordinates (the first will be placed at 10, the second at 40)
            double north = 10, east = 10, south = 10, west = 10;
            for (LPort lPort : selfLoopNode.getPorts()) {
                switch (lPort.getSide()) {
                case NORTH:
                    lPort.getPosition().set(north, 0);
                    north += 30;
                    break;

                case EAST:
                    lPort.getPosition().set(selfLoopNode.getSize().x, east);
                    east += 30;
                    break;

                case SOUTH:
                    lPort.getPosition().set(south, selfLoopNode.getSize().y);
                    south += 30;
                    break;

                case WEST:
                    lPort.getPosition().set(0, west);
                    west += 30;
                    break;

                default:
                    fail();
                }
            }
        }

        private void addEdges() {
            List<LPort> ports = selfLoopNode.getPorts();

            SelfLoopTestGraphCreator.edge(ports.get(0), ports.get(1));
            SelfLoopTestGraphCreator.edge(ports.get(1), ports.get(2));
            SelfLoopTestGraphCreator.edge(ports.get(3), ports.get(4));
            SelfLoopTestGraphCreator.edge(ports.get(4), ports.get(5));

            regularPorts.add(ports.get(0));

            onlySelfLoopPorts.add(ports.get(1));
            onlySelfLoopPorts.add(ports.get(2));
            onlySelfLoopPorts.add(ports.get(3));
            onlySelfLoopPorts.add(ports.get(4));
            onlySelfLoopPorts.add(ports.get(5));

            unconnectedPorts.add(ports.get(6));
            unconnectedPorts.add(ports.get(7));
        }
    }

}
