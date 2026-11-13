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
package org.eclipse.elk.alg.layered.intermediate.loops.ordering;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.intermediate.SelfLoopPreProcessor;
import org.eclipse.elk.alg.layered.intermediate.loops.SelfHyperLoop;
import org.eclipse.elk.alg.layered.intermediate.loops.SelfLoopHolder;
import org.eclipse.elk.alg.layered.intermediate.loops.SelfLoopTestGraphCreator;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.layered.options.SelfLoopOrderingStrategy;
import org.eclipse.elk.core.options.PortConstraints;
import org.eclipse.elk.core.util.BasicProgressMonitor;
import org.junit.Test;

/**
 * Tests how {@link PortRestorer} splits the ports of a stacked one-sided self loop. The ports are sorted by net flow,
 * and the second group starts at the first port with a positive net flow, else at the first port with a non-negative
 * net flow, else at the center port, but never at the first or the last port.
 */
public class PortRestorerTest {

    /**
     * Creates a node with the given number of northern ports and connects them by self loops between the given pairs
     * of port indices.
     */
    private static LNode node(final int portCount, final int[][] edges) {
        LGraph graph = new LGraph();
        LNode node = SelfLoopTestGraphCreator.node(graph);
        node.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        node.setProperty(LayeredOptions.EDGE_ROUTING_SELF_LOOP_ORDERING, SelfLoopOrderingStrategy.STACKED);
        SelfLoopTestGraphCreator.ports(node, portCount, 0, 0, 0);
        for (int[] edge : edges) {
            SelfLoopTestGraphCreator.edge(node.getPorts().get(edge[0]), node.getPorts().get(edge[1]));
        }
        return node;
    }

    /** Hides the node's self loop ports as the pre-processor does and restores them. */
    private static void restore(final LNode node) {
        new SelfLoopPreProcessor().process(node.getGraph(), new BasicProgressMonitor());
        SelfLoopHolder slHolder = node.getProperty(InternalProperties.SELF_LOOP_HOLDER);
        assertTrue(slHolder.arePortsHidden());
        assertTrue(node.getPorts().isEmpty());

        slHolder.getSLHyperLoops().forEach(SelfHyperLoop::computePortsPerSide);
        new PortRestorer().restorePorts(slHolder, new BasicProgressMonitor());
    }

    @Test
    public void testSplitAtFirstNonNegativeNetFlow() {
        // a -> b -> c -> d: net flows -1, 0, 0, 1. The only port with a positive net flow is the last one, so the
        // second group starts at b, the first port with a non-negative net flow. The center port would be c.
        LNode node = node(4, new int[][] { { 0, 1 }, { 1, 2 }, { 2, 3 } });
        List<LPort> p = new ArrayList<>(node.getPorts());
        restore(node);
        assertEquals(Arrays.asList(p.get(0), p.get(3), p.get(2), p.get(1)), node.getPorts());
    }

    @Test
    public void testSplitAtFirstPositiveNetFlow() {
        // a -> b, a -> c, d -> c: net flows -2, 1, 2, -1, sorted a, d, b, c. The second group starts at b.
        LNode node = node(4, new int[][] { { 0, 1 }, { 0, 2 }, { 3, 2 } });
        List<LPort> p = new ArrayList<>(node.getPorts());
        restore(node);
        assertEquals(Arrays.asList(p.get(3), p.get(0), p.get(2), p.get(1)), node.getPorts());
    }
}
