/*******************************************************************************
 * Copyright (c) 2014, 2015 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.layered.compound;

import java.util.Comparator;

import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.options.PortType;

/**
 * Compares cross-hierarchy edge segments such that they can be sorted from the start to the end
 * segment.
 * 
 * @author msp
 */
final class CrossHierarchyEdgeComparator implements Comparator<CrossHierarchyEdge> {
    private final LGraph graph;
    
    /**
     * Creates a new comparator for sorting cross-hierarchy edge segments for the given top-level
     * compound graph.
     * 
     * @param graph the top-level compound graph.
     */
    CrossHierarchyEdgeComparator(final LGraph graph) {
        this.graph = graph;
    }
    
    @Override
    public int compare(final CrossHierarchyEdge edge1, final CrossHierarchyEdge edge2) {
        if (edge1.getType() == PortType.OUTPUT
                && edge2.getType() == PortType.INPUT) {
            return -1;
        } else if (edge1.getType() == PortType.INPUT
                && edge2.getType() == PortType.OUTPUT) {
            return 1;
        }
        int level1 = hierarchyLevel(edge1.getGraph(), graph);
        int level2 = hierarchyLevel(edge2.getGraph(), graph);
        if (edge1.getType() == PortType.OUTPUT) {
            // from deeper level to higher level
            return level2 - level1;
        } else {
            // from higher level to deeper level
            return level1 - level2;
        }
    }
    
    /**
     * Compute the hierarchy level of the given nested graph.
     * 
     * @param nestedGraph a nested graph
     * @param topLevelGraph the top-level graph
     * @return the hierarchy level (higher number means the node is nested deeper)
     */
    private static int hierarchyLevel(final LGraph nestedGraph, final LGraph topLevelGraph) {
        LGraph currentGraph = nestedGraph;
        int level = 0;
        do {
            if (currentGraph == topLevelGraph) {
                return level;
            }
            LNode currentNode = currentGraph.getParentNode();
            if (currentNode == null) {
                // the given node is not an ancestor of the graph node
                throw new IllegalArgumentException();
            }
            currentGraph = currentNode.getGraph();
            level++;
        } while (true);
    }
}