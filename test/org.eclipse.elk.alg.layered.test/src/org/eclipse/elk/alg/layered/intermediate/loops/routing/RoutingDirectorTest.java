/*******************************************************************************
 * Copyright (c) 2026 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.alg.layered.intermediate.loops.routing;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertSame;

import java.lang.reflect.Field;
import java.util.EnumSet;

import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.intermediate.loops.SelfHyperLoop;
import org.eclipse.elk.alg.layered.intermediate.loops.SelfLoopHolder;
import org.eclipse.elk.alg.layered.intermediate.loops.SelfLoopPort;
import org.eclipse.elk.alg.layered.intermediate.loops.SelfLoopTestGraphCreator;
import org.eclipse.elk.core.options.PortSide;
import org.junit.Test;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;

/**
 * Tests whether {@link RoutingDirector} determines opposing two-side self-loop routes independent of multimap key
 * iteration order.
 */
public class RoutingDirectorTest {

    @Test
    public void testOpposingTwoSideLoopUsesDeterministicSideOrder() throws Exception {
        LGraph graph = new LGraph();
        LNode node = SelfLoopTestGraphCreator.node(graph);

        // Create one unconnected port on each side in the default clockwise order. This makes the penalties for the two
        // possible east-west loop routes equal, so the chosen route depends only on side ordering.
        SelfLoopTestGraphCreator.ports(node, 1, 1, 1, 1);

        LPort eastPort = node.getPorts().get(1);
        LPort westPort = node.getPorts().get(3);
        SelfLoopTestGraphCreator.edge(eastPort, westPort);

        SelfLoopHolder slHolder = SelfLoopHolder.install(node);
        SelfHyperLoop slLoop = slHolder.getSLHyperLoops().get(0);
        slLoop.computePortsPerSide();

        // Simulate an arbitrary key-set iteration order opposite to the deterministic clockwise order.
        ListMultimap<PortSide, SelfLoopPort> reversedSides = ArrayListMultimap.create();
        reversedSides.putAll(PortSide.WEST, slLoop.getSLPortsBySide(PortSide.WEST));
        reversedSides.putAll(PortSide.EAST, slLoop.getSLPortsBySide(PortSide.EAST));
        setPortsBySide(slLoop, reversedSides);

        new RoutingDirector().determineLoopRoutes(slHolder);

        assertSame(eastPort, slLoop.getLeftmostPort().getLPort());
        assertSame(westPort, slLoop.getRightmostPort().getLPort());
        assertEquals(EnumSet.of(PortSide.EAST, PortSide.SOUTH, PortSide.WEST), slLoop.getOccupiedPortSides());
    }

    private static void setPortsBySide(final SelfHyperLoop slLoop, final ListMultimap<PortSide, SelfLoopPort> portsBySide)
            throws Exception {
        Field field = SelfHyperLoop.class.getDeclaredField("slPortsBySide");
        field.setAccessible(true);
        field.set(slLoop, portsBySide);
    }
}
