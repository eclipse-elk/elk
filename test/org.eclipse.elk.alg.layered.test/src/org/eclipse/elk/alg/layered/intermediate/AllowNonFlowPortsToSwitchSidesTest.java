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
package org.eclipse.elk.alg.layered.intermediate;

import static org.junit.Assert.assertEquals;

import java.awt.geom.Line2D;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;

import org.eclipse.elk.alg.layered.LayeredLayoutProvider;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.test.PlainJavaInitialization;
import org.eclipse.elk.core.math.KVectorChain;
import org.eclipse.elk.core.options.PortConstraints;
import org.eclipse.elk.core.options.PortSide;
import org.eclipse.elk.core.util.ElkUtil;
import org.eclipse.elk.core.util.NullElkProgressMonitor;
import org.eclipse.elk.graph.ElkEdge;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.ElkPort;
import org.eclipse.elk.graph.util.ElkGraphUtil;
import org.junit.Before;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameter;
import org.junit.runners.Parameterized.Parameters;

import com.google.common.collect.Lists;

/**
 * Tests that ports properly switch sides when they are allowed and an edge crossing can be prevented.
 */
@RunWith(Parameterized.class)
public class AllowNonFlowPortsToSwitchSidesTest {

    @Parameters(name = "{0} {1} {2} -> {3}")
    public static Collection<Object[]> data() {
        return Arrays.asList(new Object[][] {
                // Arguments: 
                // port constraints, node1MaySwitchPortSides, node2MaySwitchPortSides, expected: edgesIntersect
                { PortConstraints.FIXED_SIDE, false, false, true },
                { PortConstraints.FIXED_POS, true, true, true },
                
                // As soon as one port is allowed to switch sides, the edge crossing can be prevented
                { PortConstraints.FIXED_SIDE, true, false, false },
                { PortConstraints.FIXED_SIDE, false, true, false }, 
                { PortConstraints.FIXED_SIDE, true, true, false },
                
                // TODO The behavior for FIXED_ORDER is not clearly defined, currently ports will switch if allowed
                { PortConstraints.FIXED_ORDER, true, true, false },
        });
    }

    LayeredLayoutProvider layeredLayout;
    
    @Parameter(0)
    public PortConstraints portConstraints;
    @Parameter(1)
    public boolean node1MaySwitchPortSides;
    @Parameter(2)
    public boolean node2MaySwitchPortSides;
    @Parameter(3)
    public boolean edgesIntersect;
    
    @Before
    public void setUp() {
        PlainJavaInitialization.initializePlainJavaLayout();
        layeredLayout = new LayeredLayoutProvider();
    }
    
    @Test
    public void testOnlySwitchPortSidesIfPermitted() {
        final ElkNode graph = createSimpleGraph(portConstraints, node1MaySwitchPortSides, node2MaySwitchPortSides);
        final ElkEdge e1 = graph.getContainedEdges().get(0);
        final ElkEdge e2 = graph.getContainedEdges().get(1);
        layeredLayout.layout(graph, new NullElkProgressMonitor());
        assertEquals(edgesIntersect, intersect(e1, e2));
    }
    
    /**
     * Create a simple test graph.
     * 
     * Note that port positions are set for the test case FIXED_POS.
     */
    private static ElkNode createSimpleGraph(final PortConstraints portConstraints,
            final boolean node1MaySwitchPortSides, final boolean node2MaySwitchPortSides) {
        ElkNode graph = ElkGraphUtil.createGraph();

        ElkNode node1 = ElkGraphUtil.createNode(graph);
        node1.setWidth(30);
        node1.setHeight(30);
        node1.setProperty(LayeredOptions.PORT_CONSTRAINTS, portConstraints);
        ElkPort portNorth1 = ElkGraphUtil.createPort(node1);
        portNorth1.setX(15.0);
        portNorth1.setY(-1.0);
        portNorth1.setProperty(LayeredOptions.PORT_SIDE, PortSide.NORTH);
        portNorth1.setProperty(LayeredOptions.ALLOW_NON_FLOW_PORTS_TO_SWITCH_SIDES, node1MaySwitchPortSides);
        ElkPort portSouth1 = ElkGraphUtil.createPort(node1);
        portSouth1.setX(15.0);
        portSouth1.setY(30.0);
        portSouth1.setProperty(LayeredOptions.PORT_SIDE, PortSide.SOUTH);
        portSouth1.setProperty(LayeredOptions.ALLOW_NON_FLOW_PORTS_TO_SWITCH_SIDES, node1MaySwitchPortSides);

        ElkNode node2 = ElkGraphUtil.createNode(graph);
        node2.setWidth(30);
        node2.setHeight(30);
        node2.setProperty(LayeredOptions.PORT_CONSTRAINTS, portConstraints);
        ElkPort portNorth2 = ElkGraphUtil.createPort(node2);
        portNorth2.setX(15.0);
        portNorth2.setY(-1.0);
        portNorth2.setProperty(LayeredOptions.PORT_SIDE, PortSide.NORTH);
        portNorth2.setProperty(LayeredOptions.ALLOW_NON_FLOW_PORTS_TO_SWITCH_SIDES, node2MaySwitchPortSides);
        ElkPort portSouth2 = ElkGraphUtil.createPort(node2);
        portSouth2.setX(15.0);
        portSouth2.setY(30.0);
        portSouth2.setProperty(LayeredOptions.PORT_SIDE, PortSide.SOUTH);
        portSouth2.setProperty(LayeredOptions.ALLOW_NON_FLOW_PORTS_TO_SWITCH_SIDES, node2MaySwitchPortSides);

        // Create two edges that, with fixed port sides, will have to cross.
        ElkGraphUtil.createSimpleEdge(portNorth1, portSouth2);
        ElkGraphUtil.createSimpleEdge(portSouth1, portNorth2);

        return graph;
    }

    private static List<Line2D> edgeToSegments(final ElkEdge e) {
        KVectorChain chain = ElkUtil.createVectorChain(ElkGraphUtil.firstEdgeSection(e, false, false));
        List<Line2D> segments = Lists.newArrayList();
        chain.stream().reduce((prev, cur) -> {
            segments.add(new Line2D.Double(prev.x, prev.y, cur.x, cur.y));
            return cur;
        });
        return segments;
    }
    
    /**
     * @return true if the two edges' paths cross at some point. Note that the implementation assumes straight edge
     *         segments, i.e. does not necessarily work for splines.
     */
    public static boolean intersect(final ElkEdge e1, final ElkEdge e2) {
        List<Line2D> segments1 = edgeToSegments(e1);
        List<Line2D> segments2 = edgeToSegments(e2);
        for (Line2D s1 : segments1) {
            for (Line2D s2 : segments2) {
                if (s1.intersectsLine(s2)) {
                    return true;
                }
            }
        }
        return false;
    }
}
