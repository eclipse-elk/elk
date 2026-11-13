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
package org.eclipse.elk.alg.layered.p5edges.splines;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertSame;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.Collections;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Set;

import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.test.framework.LayoutTestRunner;
import org.eclipse.elk.alg.test.framework.annotations.Algorithm;
import org.eclipse.elk.alg.test.framework.annotations.DefaultConfiguration;
import org.eclipse.elk.alg.test.framework.annotations.GraphProvider;
import org.eclipse.elk.alg.test.framework.annotations.TestAfterProcessor;
import org.eclipse.elk.core.options.EdgeRouting;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.ElkPort;
import org.eclipse.elk.graph.util.ElkGraphUtil;
import org.junit.runner.RunWith;

/**
 * Tests that the edges of a spline hyperedge segment are kept in the order in which the {@link SplineEdgeRouter}
 * finds them at the hyperedge's port. The first of them represents the segment, for example as the edge of its
 * vertical segment during compaction. With hash-ordered sets, the representative edge depended on identity hashes.
 */
@RunWith(LayoutTestRunner.class)
@Algorithm(LayeredOptions.ALGORITHM_ID)
@DefaultConfiguration()
public class SplineSegmentEdgeOrderTest {

    /** Number of edges leaving the source port; about half of them run upwards and half downwards. */
    private static final int EDGE_COUNT = 10;

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Sources

    /**
     * Creates a node whose only port connects to {@link #EDGE_COUNT} nodes of the next layer, routed as splines.
     */
    @GraphProvider
    public ElkNode produceGraph() {
        ElkNode graph = ElkGraphUtil.createGraph();
        graph.setProperty(LayeredOptions.EDGE_ROUTING, EdgeRouting.SPLINES);

        ElkNode source = node(graph);
        ElkPort sourcePort = ElkGraphUtil.createPort(source);
        for (int i = 0; i < EDGE_COUNT; i++) {
            ElkGraphUtil.createSimpleEdge(sourcePort, ElkGraphUtil.createPort(node(graph)));
        }

        return graph;
    }

    private static ElkNode node(final ElkNode graph) {
        ElkNode node = ElkGraphUtil.createNode(graph);
        node.setDimensions(30, 30);
        return node;
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Tests

    /**
     * Checks for each hyperedge segment that its edges iterate in the order of the port's connected edges, and so
     * that the representative edge is the first of them.
     */
    @TestAfterProcessor(SplineEdgeRouter.class)
    public void testHyperedgeSegmentEdgeOrder(final Object graph) {
        int hyperedgeSegments = 0;
        for (SplineSegment segment : hyperedgeSegments((LGraph) graph)) {
            List<LEdge> segmentEdges = new ArrayList<>(segment.edges);
            LPort port = segmentEdges.get(0).getSource();

            List<LEdge> expectedEdges = new ArrayList<>();
            for (LEdge edge : port.getConnectedEdges()) {
                if (segment.edges.contains(edge)) {
                    expectedEdges.add(edge);
                }
            }

            assertEquals(segment.edges.size(), expectedEdges.size());
            assertSame(expectedEdges.get(0), segment.edges.iterator().next());
            assertEquals(expectedEdges, segmentEdges);
            hyperedgeSegments++;
        }

        // The upward and the downward edges each form a hyperedge
        assertEquals(2, hyperedgeSegments);
    }

    /** Collects the distinct segments of the graph's splines that represent more than one edge. */
    private static List<SplineSegment> hyperedgeSegments(final LGraph lGraph) {
        Set<SplineSegment> seen = Collections.newSetFromMap(new IdentityHashMap<>());
        List<SplineSegment> segments = new ArrayList<>();
        for (Layer layer : lGraph) {
            for (LNode node : layer) {
                for (LEdge edge : node.getOutgoingEdges()) {
                    List<SplineSegment> spline = edge.getProperty(InternalProperties.SPLINE_ROUTE_START);
                    assertTrue(spline != null);
                    for (SplineSegment segment : spline) {
                        if (segment.edges.size() > 1 && seen.add(segment)) {
                            segments.add(segment);
                        }
                    }
                }
            }
        }
        return segments;
    }
}
