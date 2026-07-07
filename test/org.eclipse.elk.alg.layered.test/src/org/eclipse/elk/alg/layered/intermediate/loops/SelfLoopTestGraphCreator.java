/*******************************************************************************
 * Copyright (c) 2019 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.layered.intermediate.loops;

import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.core.options.PortSide;

/**
 * Creates various test graphs with (and without) self loops.
 */
public final class SelfLoopTestGraphCreator {
    
    private SelfLoopTestGraphCreator() {
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Creation
    
    /**
     * Creates a graph with two nodes connected by an edge.
     */
    public static LGraph basicGraphWithoutSelfLoops() {
        LGraph lGraph = new LGraph();
        
        LNode n1 = node(lGraph);
        ports(n1, 0, 1, 0, 0);
        
        LNode n2 = node(lGraph);
        ports(n2, 0, 0, 0, 1);
        
        edge(n1.getPorts().get(0), n2.getPorts().get(0));
        
        return lGraph;
    }
    
    /**
     * Creates a graph with two nodes. The first is a regular node, the second contains self loops.
     */
    public static LGraph basicGraphWithSelfLoops() {
        LGraph lGraph = basicGraphWithoutSelfLoops();
        
        // Add more stuff to the second node
        LNode lNode = lGraph.getLayerlessNodes().get(1);
        ports(lNode, 2, 2, 2, 1);
        
        // There are always pairs of ports. Connect each pair
        for (int i = 1; i < lNode.getPorts().size(); i += 2) {
            edge(lNode.getPorts().get(i - 1), lNode.getPorts().get(i));
        }
        
        return lGraph;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Utilities
    
    public static LNode node(final LGraph lGraph) {
        LNode lNode = new LNode(lGraph);
        lGraph.getLayerlessNodes().add(lNode);
        return lNode;
    }
    
    public static LEdge edge(final LPort source, final LPort target) {
        LEdge lEdge = new LEdge();
        lEdge.setSource(source);
        lEdge.setTarget(target);
        return lEdge;
    }
    
    public static void ports(final LNode lNode, final int north, final int east, final int south, final int west) {
        for (int i = 0; i < north; i++) {
            LPort lPort = new LPort();
            lPort.setNode(lNode);
            lPort.setSide(PortSide.NORTH);
        }

        for (int i = 0; i < east; i++) {
            LPort lPort = new LPort();
            lPort.setNode(lNode);
            lPort.setSide(PortSide.EAST);
        }

        for (int i = 0; i < south; i++) {
            LPort lPort = new LPort();
            lPort.setNode(lNode);
            lPort.setSide(PortSide.SOUTH);
        }

        for (int i = 0; i < west; i++) {
            LPort lPort = new LPort();
            lPort.setNode(lNode);
            lPort.setSide(PortSide.WEST);
        }
    }

}
