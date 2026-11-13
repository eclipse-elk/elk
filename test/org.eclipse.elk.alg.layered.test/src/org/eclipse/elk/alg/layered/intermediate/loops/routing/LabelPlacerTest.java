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
package org.eclipse.elk.alg.layered.intermediate.loops.routing;

import static org.junit.Assert.assertEquals;

import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LLabel;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.intermediate.loops.SelfHyperLoop;
import org.eclipse.elk.alg.layered.intermediate.loops.SelfHyperLoopLabels;
import org.eclipse.elk.alg.layered.intermediate.loops.SelfHyperLoopLabels.Alignment;
import org.eclipse.elk.alg.layered.intermediate.loops.SelfLoopHolder;
import org.eclipse.elk.alg.layered.intermediate.loops.SelfLoopTestGraphCreator;
import org.eclipse.elk.alg.layered.intermediate.loops.SelfLoopType;
import org.eclipse.elk.core.options.PortSide;
import org.eclipse.elk.core.util.BasicProgressMonitor;
import org.junit.Test;

/**
 * Tests where {@link LabelPlacer} puts the labels of four-sided self loops. Such a loop leaves a gap between its
 * rightmost and its leftmost port. Its label goes on the southern side if either of those ports is on the northern
 * side, and on the northern side otherwise.
 */
public class LabelPlacerTest {

    /**
     * Creates a node with one port per given side, in that order, and a labelled self hyper loop that connects the
     * ports at the given indices. Then routes the loop and places its label.
     */
    private static SelfHyperLoop placeLabel(final PortSide[] sides, final int... connectedPorts) {
        LGraph graph = new LGraph();
        LNode node = SelfLoopTestGraphCreator.node(graph);
        for (PortSide side : sides) {
            LPort port = new LPort();
            port.setNode(node);
            port.setSide(side);
        }

        for (int i = 1; i < connectedPorts.length; i++) {
            SelfLoopTestGraphCreator.edge(
                    node.getPorts().get(connectedPorts[i - 1]), node.getPorts().get(connectedPorts[i]));
        }
        LEdge labelledEdge = node.getPorts().get(connectedPorts[0]).getOutgoingEdges().get(0);
        labelledEdge.getLabels().add(new LLabel("label"));

        SelfLoopHolder slHolder = SelfLoopHolder.install(node);
        SelfHyperLoop slLoop = slHolder.getSLHyperLoops().get(0);
        slLoop.computePortsPerSide();
        assertEquals(SelfLoopType.FOUR_SIDES, slLoop.getSelfLoopType());

        new RoutingDirector().determineLoopRoutes(slHolder);
        new LabelPlacer().placeLabels(slHolder, null, new BasicProgressMonitor());
        return slLoop;
    }

    private static void assertLoop(final SelfHyperLoop slLoop, final PortSide leftmostSide,
            final PortSide rightmostSide, final PortSide labelSide) {

        assertEquals(leftmostSide, slLoop.getLeftmostPort().getLPort().getSide());
        assertEquals(rightmostSide, slLoop.getRightmostPort().getLPort().getSide());

        SelfHyperLoopLabels slLabels = slLoop.getSLLabels();
        assertEquals(labelSide, slLabels.getSide());
        assertEquals(Alignment.CENTER, slLabels.getAlignment());
    }

    @Test
    public void testFourSidesLoopEndingOnNorthernSide() {
        // The unconnected northern port makes the loop leave its gap between the northern and the eastern port, so
        // the loop starts on the eastern side and ends on the northern side, which it does not cover completely
        SelfHyperLoop slLoop = placeLabel(
                new PortSide[] { PortSide.NORTH, PortSide.NORTH, PortSide.EAST, PortSide.SOUTH, PortSide.WEST },
                0, 2, 3, 4);
        assertLoop(slLoop, PortSide.EAST, PortSide.NORTH, PortSide.SOUTH);
    }

    @Test
    public void testFourSidesLoopStartingOnNorthernSide() {
        // The unconnected western port puts the gap between the western and the northern port
        SelfHyperLoop slLoop = placeLabel(
                new PortSide[] { PortSide.NORTH, PortSide.EAST, PortSide.SOUTH, PortSide.WEST, PortSide.WEST },
                0, 1, 2, 3);
        assertLoop(slLoop, PortSide.NORTH, PortSide.WEST, PortSide.SOUTH);
    }

    @Test
    public void testFourSidesLoopCoveringNorthernSide() {
        // The unconnected eastern port puts the gap between the eastern and the southern port
        SelfHyperLoop slLoop = placeLabel(
                new PortSide[] { PortSide.NORTH, PortSide.EAST, PortSide.EAST, PortSide.SOUTH, PortSide.WEST },
                0, 1, 3, 4);
        assertLoop(slLoop, PortSide.SOUTH, PortSide.EAST, PortSide.NORTH);
    }
}
