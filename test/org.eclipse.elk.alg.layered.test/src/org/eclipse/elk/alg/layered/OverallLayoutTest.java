/*******************************************************************************
 * Copyright (c) 2012, 2020 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.layered;

import static org.junit.Assert.assertTrue;

import java.util.Iterator;
import java.util.LinkedList;

import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.test.PlainJavaInitialization;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.options.SizeConstraint;
import org.eclipse.elk.core.util.BasicProgressMonitor;
import org.eclipse.elk.core.util.IElkProgressMonitor;
import org.eclipse.elk.graph.ElkEdge;
import org.eclipse.elk.graph.ElkEdgeSection;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.util.ElkGraphUtil;
import org.junit.Before;
import org.junit.Test;

/**
 * Test the overall layout process including graph import and layout transfer.
 */
public class OverallLayoutTest {
    private IElkProgressMonitor simpleMonitor;
    private ElkNode simpleGraph;

    /**
     * Set up the test class.
     */
    @Before
    public void setUp() {
        PlainJavaInitialization.initializePlainJavaLayout();
        LayeredLayoutProvider layoutProvider = new LayeredLayoutProvider();
        simpleGraph = createSimpleGraph();
        layoutProvider.layout(simpleGraph, new BasicProgressMonitor());
    }

    /**
     * Tests if the coordinates of the graph's nodes are positive. Doesn't work for hierarchical
     * graphs.
     */
    @Test
    public void testNodeCoordinates() {
        for (ElkNode node : simpleGraph.getChildren()) {
            assertTrue(node.getX() > 0);
            assertTrue(node.getY() > 0);
        }
    }

    /**
     * Tests if the source and target coordinates of all edges are > 0. Doesn't work for
     * hierarchical graphs.
     */
    @Test
    public void testEdgeCoordinates() {
        for (ElkEdge edge : simpleGraph.getContainedEdges()) {
            for (ElkEdgeSection section : edge.getSections()) {
                assertTrue(section.getStartX() > 0);
                assertTrue(section.getStartY() > 0);
                assertTrue(section.getEndX() > 0);
                assertTrue(section.getEndY() > 0);
            }
        }
    }

    /**
     * Tests if the graph size is positive.
     */
    @Test
    public void testGraphSize() {
        assertTrue(simpleGraph.getWidth() > 0);
        assertTrue(simpleGraph.getHeight() > 0);
    }

    /**
     * Tests if the edge routing produced perfectly orthogonal edges. This test only makes sense if
     * the orthogonal edge router is used.
     */
    @Test
    public void testEdgeOrthogonality() {
        LinkedList<KVector> bendPoints = new LinkedList<>();
        
        for (ElkEdge edge : simpleGraph.getContainedEdges()) {
            for (ElkEdgeSection section : edge.getSections()) {
                // Assemble list of bend points
                bendPoints.clear();
                bendPoints.add(new KVector(section.getStartX(), section.getStartY()));
                section.getBendPoints().stream()
                    .map(bp -> new KVector(bp.getX(), bp.getY()))
                    .forEach(vec -> bendPoints.add(vec));
                bendPoints.add(new KVector(section.getEndX(), section.getEndY()));
                
                if (bendPoints.size() >= 2) {
                    // Check whether the bend points run orthogonally
                    Iterator<KVector> bendPointIterator = bendPoints.iterator();
                    
                    KVector prevBendPoint = bendPointIterator.next();
                    while (bendPointIterator.hasNext()) {
                        KVector currBendPoint = bendPointIterator.next();
                        
                        assertTrue(prevBendPoint.x == currBendPoint.x
                                || prevBendPoint.y == currBendPoint.y);
                        
                        prevBendPoint = currBendPoint;
                    }
                }
            }
        }
    }

    /**
     * Create a simple test graph. The graph has at least two nodes and an edge, the nodes have
     * predefined sizes with fixed size constraint, but neither nodes nor edges have predefined
     * coordinates, i.e. they are all set to (0,0). The graph size is initially 0 as well.
     */
    private static ElkNode createSimpleGraph() {
        ElkNode graph = ElkGraphUtil.createGraph();

        ElkNode node1 = ElkGraphUtil.createNode(graph);
        node1.setWidth(30);
        node1.setHeight(30);
        node1.setProperty(LayeredOptions.NODE_SIZE_CONSTRAINTS, SizeConstraint.fixed());

        ElkNode node2 = ElkGraphUtil.createNode(graph);
        node2.setWidth(30);
        node2.setHeight(30);
        node2.setProperty(LayeredOptions.NODE_SIZE_CONSTRAINTS, SizeConstraint.fixed());

        ElkGraphUtil.createSimpleEdge(node1, node2);

        return graph;
    }

}
