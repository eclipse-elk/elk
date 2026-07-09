/*******************************************************************************
 * Copyright (c) 2011, 2015 Kiel University and others.
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
package org.eclipse.elk.alg.force.graph;

import java.util.ArrayList;
import java.util.List;

import org.eclipse.elk.alg.force.options.ForceOptions;
import org.eclipse.elk.graph.properties.MapPropertyHolder;

import com.google.common.collect.Iterables;

/**
 * A graph for the force layouter.
 * 
 * @author owo
 * @author msp
 */
public final class FGraph extends MapPropertyHolder {
    
    /** the serial version UID. */
    private static final long serialVersionUID = -2396315570561498425L;
    
    /** All nodes of this graph. */
    private List<FNode> nodes = new ArrayList<>();
    /** All edges of this graph. */
    private List<FEdge> edges = new ArrayList<>();
    /** All labels of this graph. */
    private List<FLabel> labels = new ArrayList<>();
    /** All bend points of this graph. */
    private List<FBendpoint> bendPoints = new ArrayList<>();
    /** adjacency matrix of the graph. */
    private int[][] adjacency;

    /**
     * Returns the list of edges for this graph.
     * 
     * @return the edges
     */
    public List<FEdge> getEdges() {
        return edges;
    }

    /**
     * Returns the list of nodes for this graph.
     * 
     * @return the nodes
     */
    public List<FNode> getNodes() {
        return nodes;
    }

    /**
     * Returns the list of labels for this graph.
     * 
     * @return the labels
     */
    public List<FLabel> getLabels() {
        return labels;
    }
    
    /**
     * Returns the list of bend points for this graph.
     * 
     * @return the bend points
     */
    public List<FBendpoint> getBendpoints() {
        return bendPoints;
    }
    
    /**
     * Returns a list of all particles occurring in the graph: vertices, labels, and bend points.
     * 
     * @return iterable over all particles
     */
    public Iterable<FParticle> getParticles() {
        return Iterables.concat(nodes, labels, bendPoints);
    }
    
    /**
     * Determines the amount of connection between the two given particles, considering
     * the priority value of the respective edges.
     * 
     * @param particle1 first particle
     * @param particle2 second particle
     * @return the amount of connection
     */
    public int getConnection(final FParticle particle1, final FParticle particle2) {
        if (particle1 instanceof FNode && particle2 instanceof FNode) {
            FNode node1 = (FNode) particle1, node2 = (FNode) particle2;
            return adjacency[node1.id][node2.id] + adjacency[node2.id][node1.id];
        } else if (particle1 instanceof FBendpoint && particle2 instanceof FBendpoint) {
            FBendpoint bpoint1 = (FBendpoint) particle1, bpoint2 = (FBendpoint) particle2;
            if (bpoint1.getEdge() == bpoint2.getEdge()) {
                return bpoint2.getEdge().getProperty(ForceOptions.PRIORITY);
            }
        }
        return 0;
    }
    
    /**
     * Calculate the adjacency matrix of the graph.
     */
    public void calcAdjacency() {
        int n = nodes.size();
        adjacency = new int[n][n];
        for (FEdge edge : edges) {
            adjacency[edge.getSource().id][edge.getTarget().id] += edge.getProperty(ForceOptions.PRIORITY);
        }
    }

}
