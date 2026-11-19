/*******************************************************************************
 * Copyright (c) 2016 Kiel University and others.
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
package org.eclipse.elk.alg.force.stress;

import java.util.Arrays;
import java.util.List;
import java.util.PriorityQueue;
import java.util.Queue;

import org.eclipse.elk.alg.force.graph.FEdge;
import org.eclipse.elk.alg.force.graph.FGraph;
import org.eclipse.elk.alg.force.graph.FNode;
import org.eclipse.elk.alg.force.options.StressOptions;
import org.eclipse.elk.core.math.KVector;

import com.google.common.collect.LinkedListMultimap;
import com.google.common.collect.Multimap;

/**
 * Implementation of stress minimizing layout as described by Gansner, Koren, and North.
 * <ul><li>
 * Emden Gansner, Yehuda Koren, and Stephen North. Graph drawing by stress majorization. <em>Graph Drawing</em>, 2005.
 * </li></ul>
 * 
 * The implementation supports performing a layout in one dimension only, preserving the coordinates of the other
 * dimension. For this, set {@link StressOptions#DIMENSION} to either {@link Dimension#X} or {@link Dimension#Y}.
 * Furthermore, nodes can be fixed using the {@link StressOptions#FIXED} option.
 */
public class StressMajorization {

    /** The graph do be laid out, should be connected. */
    private FGraph graph;

    /** All pairs shortest path matrix. */
    private double[][] apsp;
    /** Weights for each pair of nodes. */
    private double[][] w;

    /** Common desired edge length, can be overridden by individual edges. */
    private double desiredEdgeLength;
    /** Dimensions to consider during layout. */
    private Dimension dim;
    /** Epsilon for terminating the stress minimizing process. */
    private double epsilon;
    /** Maximum number of iterations (overrides the {@link #epsilon}). */
    private int iterationLimit;

    private Multimap<FNode, FEdge> connectedEdges = LinkedListMultimap.create();


    /**
     * Initialize all internal structures that are required for the subsequent iterative procedure.. 
     * 
     * @param fgraph the graph to be laid out. 
     */
    public void initialize(final FGraph fgraph) {
        if (fgraph.getNodes().size() <= 1) {
            return;
        }

        this.graph = fgraph;

        this.dim = graph.getProperty(StressOptions.DIMENSION);
        this.iterationLimit = graph.getProperty(StressOptions.ITERATION_LIMIT);
        this.epsilon = graph.getProperty(StressOptions.EPSILON);
        this.desiredEdgeLength = graph.getProperty(StressOptions.DESIRED_EDGE_LENGTH);
        
        connectedEdges.clear();
        for (FEdge edge : graph.getEdges()) {
            connectedEdges.put(edge.getSource(), edge);
            connectedEdges.put(edge.getTarget(), edge);
        }

        // all pairs shortest path
        int n = graph.getNodes().size();
        apsp = new double[n][n];
        for (FNode source : graph.getNodes()) {
            dijkstra(source, apsp[source.id]);
        }

        // init weight matrix
        w = new double[n][n];
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                double dij = apsp[i][j];
                double wij = 1.0 / (dij * dij);
                w[i][j] = wij;
            }
        }
    }

    /**
     * Execute the stress-minimizing iteration until a termination criterion is reached. 
     */
    public void execute() {
        if (graph.getNodes().size() <= 1) {
            return;
        }
        
        int count = 0;
        double prevStress = computeStress();
        double curStress = Double.POSITIVE_INFINITY;

        do {
            if (count > 0) {
                prevStress = curStress;
            }

            for (FNode u : graph.getNodes()) {

                // note that we do not use 'NO_LAYOUT' here,
                // since that option results in the node already
                // being excluded by the layout engine
                if (u.getProperty(StressOptions.FIXED)) {
                    continue;
                }

                KVector newPos = computeNewPosition(u);
                u.getPosition().reset().add(newPos);
            }

            curStress = computeStress();
            
        } while (!done(count++, prevStress, curStress));
    }

    /**
     * Performs Dijkstra's all pairs shortest path algorithm.
     */
    private void dijkstra(final FNode source, final double[] dist) {
        Queue<FNode> nodes = new PriorityQueue<FNode>((n1, n2) -> Double.compare(dist[n1.id], dist[n2.id]));
        boolean[] mark = new boolean[graph.getNodes().size()];
        
        // init
        Arrays.fill(mark, false);
        dist[source.id] = 0;
        for (FNode n : graph.getNodes()) {
            if (n.id != source.id) {
                dist[n.id] = Integer.MAX_VALUE;
            }
            nodes.add(n);
        }

        // find shortest paths
        while (!nodes.isEmpty()) {
            FNode u = nodes.poll();
            mark[u.id] = true;

            for (FEdge e : connectedEdges.get(u)) {
                FNode v = getOther(e, u);
                if (mark[v.id]) {
                    continue;
                }
                // get e's desired length
                double el;
                if (e.hasProperty(StressOptions.DESIRED_EDGE_LENGTH)) {
                    el = e.getProperty(StressOptions.DESIRED_EDGE_LENGTH);
                } else { 
                    el = desiredEdgeLength;
                }
                double d = dist[u.id] + el;
                if (d < dist[v.id]) {
                    dist[v.id] = d;
                    nodes.remove(v);
                    nodes.add(v);
                }
            }
        }
    }

    /**
     * Done if either stress improvement is small than {@link StressOptions#EPSILON} or the
     * {@link StressOptions#ITERATION_LIMIT} is reached.
     */
    private boolean done(final int count, final double prevStress, final double curStress) {
        return prevStress == 0 
            || (((prevStress - curStress) / prevStress) < epsilon) 
            || (count >= iterationLimit);
    }

    /**
     * @return the stress value of the current node positioning.
     */
    private double computeStress() {
        double stress = 0;
        List<FNode> nodes = graph.getNodes();
        // we know 'nodes' is an arraylist
        for (int i = 0; i < nodes.size(); ++i) {
            FNode u = nodes.get(i);
            for (int j = i + 1; j < nodes.size(); ++j) {
                FNode v = nodes.get(j);
                double eucDist = u.getPosition().distance(v.getPosition());
                double eucDisplacement = eucDist - apsp[u.id][v.id];
                stress += w[u.id][v.id] * eucDisplacement * eucDisplacement;
            }
        }
        return stress;
    }

    /**
     * Computes a new position for the passed node. The procedure is described in
     * <em>Section 2.3 Localized optimization</em> of the paper.
     */
    private KVector computeNewPosition(final FNode u) {
        double weightSum = 0;
        double xDisp = 0;
        double yDisp = 0;
        
        // we need at least two nodes here, otherwise we would divide by zero below 
        assert graph.getNodes().size() > 1;

        for (FNode v : graph.getNodes()) {
            if (u == v) {
                continue;
            }

            double wij = w[u.id][v.id];
            weightSum += wij;

            double eucDist = u.getPosition().distance(v.getPosition());

            if (eucDist > 0 && dim != Dimension.Y) {
                xDisp += wij
                        * (v.getPosition().x + apsp[u.id][v.id] * (u.getPosition().x - v.getPosition().x) / eucDist);
            }

            if (eucDist > 0 && dim != Dimension.X) {
                yDisp += wij
                        * (v.getPosition().y + apsp[u.id][v.id] * (u.getPosition().y - v.getPosition().y) / eucDist);
            }
        }

        switch (dim) {
            case X:
                return new KVector(xDisp / weightSum, u.getPosition().y);
            case Y:
                return new KVector(u.getPosition().x, yDisp / weightSum);
            default:
                return new KVector(xDisp / weightSum, yDisp / weightSum);
        }
    }

    private FNode getOther(final FEdge edge, final FNode one) {
        if (edge.getSource() == one) {
            return edge.getTarget();
        } else if (edge.getTarget() == one) {
            return edge.getSource();
        } else {
            throw new IllegalArgumentException("Node 'one' must be either source or target of edge 'edge'.");
        }
    }
    
    /**
     * Dimensions in which nodes may be moved.
     */
    public enum Dimension {
        /** Both x and y allowed. */
        XY, 
        /** Nodes may only move in x. */
        X, 
        /** Nodes may only move in y. */
        Y
    }
}
