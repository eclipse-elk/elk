/*******************************************************************************
 * Copyright (c) 2016 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.graph.util;

import static org.junit.Assert.*;

import org.eclipse.elk.graph.ElkConnectableShape;
import org.eclipse.elk.graph.ElkEdge;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.ElkPort;
import org.junit.Test;

/**
 * Tests for the {@link ElkGraphUtil} class.
 */
public class ElkGraphUtilTest {
    
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // connectableShapeToX
    
    @Test
    public void testConnectableShapeToNode() {
        ElkNode node = ElkGraphUtil.createNode(null);
        assertEquals(node, ElkGraphUtil.connectableShapeToNode(node));
        
        ElkPort port = ElkGraphUtil.createPort(node);
        assertEquals(node, ElkGraphUtil.connectableShapeToNode(port));
    }
    
    @Test(expected = NullPointerException.class)
    public void testConnectableShapeToNodeWithNullPointer() {
        ElkGraphUtil.connectableShapeToNode(null);
    }
    
    @Test
    public void testConnectableShapeToPort() {
        ElkNode node = ElkGraphUtil.createNode(null);
        assertNull(ElkGraphUtil.connectableShapeToPort(node));
        
        ElkPort port = ElkGraphUtil.createPort(node);
        assertEquals(port, ElkGraphUtil.connectableShapeToPort(port));
    }
    
    @Test(expected = NullPointerException.class)
    public void testConnectableShapeToPortWithNullPointer() {
        ElkGraphUtil.connectableShapeToPort(null);
    }
    
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // findBestEdgeContainment

    @Test
    public void testFindBestEdgeContainment() {
        // Create a basic hierarchical graph that covers all possibilities
        ElkNode graph1 = ElkGraphUtil.createGraph();
        
        ElkNode node1 = ElkGraphUtil.createNode(graph1);
        ElkNode node2 = ElkGraphUtil.createNode(graph1);
        ElkNode node3 = ElkGraphUtil.createNode(graph1);
        ElkNode node3_1 = ElkGraphUtil.createNode(node3);
        ElkNode node3_1_1 = ElkGraphUtil.createNode(node3_1);
        ElkNode node4 = ElkGraphUtil.createNode(graph1);
        ElkNode node4_1 = ElkGraphUtil.createNode(node4);
        
        // Create a second graph to test special cases later
        ElkNode graph2 = ElkGraphUtil.createGraph();
        
        ElkNode nodeA = ElkGraphUtil.createNode(graph2);
        
        // Edge which connects two top-level nodes
        ElkEdge sameLevelEdge = createEdgeWithoutContainment(node1, node2);
        assertEquals(graph1, ElkGraphUtil.findBestEdgeContainment(sameLevelEdge));
        
        // Self loop
        ElkEdge selfLoopEdge = createEdgeWithoutContainment(node1, node1);
        assertEquals(graph1, ElkGraphUtil.findBestEdgeContainment(selfLoopEdge));
        
        // Edge which connects a source on a higher level to a target on a lower level
        ElkEdge downLevelEdge = createEdgeWithoutContainment(node1,  node3_1);
        assertEquals(graph1, ElkGraphUtil.findBestEdgeContainment(downLevelEdge));
        
        // Edge which connects a source on a lower level to a target on a higher level
        ElkEdge upLevelEdge = createEdgeWithoutContainment(node3_1,  node2);
        assertEquals(graph1, ElkGraphUtil.findBestEdgeContainment(upLevelEdge));
        
        // Edge which connects a hierarchical node to its child
        ElkEdge toChildEdge = createEdgeWithoutContainment(node3,  node3_1);
        assertEquals(node3, ElkGraphUtil.findBestEdgeContainment(toChildEdge));
        
        // Edge which connects a node to its parent
        ElkEdge toParentEdge = createEdgeWithoutContainment(node3_1,  node3);
        assertEquals(node3, ElkGraphUtil.findBestEdgeContainment(toParentEdge));
        
        // Edge which connects a hierarchical node to its grand child
        ElkEdge toGrandChildEdge = createEdgeWithoutContainment(node3,  node3_1_1);
        assertEquals(node3, ElkGraphUtil.findBestEdgeContainment(toGrandChildEdge));
        
        // Edge which connects a node to its grand parent
        ElkEdge toGrandParentEdge = createEdgeWithoutContainment(node3_1_1,  node3);
        assertEquals(node3, ElkGraphUtil.findBestEdgeContainment(toGrandParentEdge));
        
        // Edge which connects to nodes on different branches of the inclusion tree
        ElkEdge crossHierarchyEdge = createEdgeWithoutContainment(node3_1,  node4_1);
        assertEquals(graph1, ElkGraphUtil.findBestEdgeContainment(crossHierarchyEdge));
        
        // Edge that connects nodes of different graphs
        ElkEdge crossGraphEdge = createEdgeWithoutContainment(node1, nodeA);
        assertNull(ElkGraphUtil.findBestEdgeContainment(crossGraphEdge));
        
        // Partially specified edges
        ElkEdge sourceMissingEdge = createEdgeWithoutContainment(null, node1);
        assertEquals(graph1, ElkGraphUtil.findBestEdgeContainment(sourceMissingEdge));

        ElkEdge targetMissingEdge = createEdgeWithoutContainment(node1, null);
        assertEquals(graph1, ElkGraphUtil.findBestEdgeContainment(targetMissingEdge));
    }
    
    @Test(expected = IllegalArgumentException.class)
    public void testFindBestEdgeContainmentWithUnconnectedEdge() {
        ElkEdge edge = ElkGraphUtil.createEdge(null);
        ElkGraphUtil.findBestEdgeContainment(edge);
    }
    
    @Test(expected = NullPointerException.class)
    public void testFindBestEdgeContainmentWithNullPointer() {
        ElkGraphUtil.findBestEdgeContainment(null);
    }
    
    /**
     * Creates an edge which connects the two shapes, but does not do so in a way that automatically sets the
     * edge's containment.
     */
    private ElkEdge createEdgeWithoutContainment(ElkConnectableShape source, ElkConnectableShape target) {
        ElkEdge edge = ElkGraphUtil.createEdge(null);
        
        if (source != null) {
            edge.getSources().add(source);
        }
        
        if (target != null) {
            edge.getTargets().add(target);
        }
        
        return edge;
    }

}
