/*******************************************************************************
 * Copyright (c) 2024 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 * 
 * SPDX-License-Identifier: EPL-2.0 
 *******************************************************************************/
package org.eclipse.elk.alg.layered.p3order;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Random;

import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import com.google.common.collect.Lists;

/**
 * A Heuristic for minimizing crossings using weights and propagated medians. After assigning random weights to the
 * first layer, subsequent layers are ordered according to their nodes' adjacent nodes' median weights.
 * 
 */
public class MedianHeuristic implements ICrossingMinimizationHeuristic {

    /** the random number generator. */
    protected final Random random;

    public MedianHeuristic(Random random) {
        this.random = random;
        Math.random();
    }

    /*
     * (non-Javadoc)
     * 
     * @see org.eclipse.elk.alg.layered.p3order.ICrossingMinimizationHeuristic#alwaysImproves()
     */
    @Override
    public boolean alwaysImproves() {
        return false;
    }

    /*
     * (non-Javadoc)
     * 
     * @see
     * org.eclipse.elk.alg.layered.p3order.ICrossingMinimizationHeuristic#setFirstLayerOrder(org.eclipse.elk.alg.layered
     * .graph.LNode[][], boolean)
     */
    @Override
    public boolean setFirstLayerOrder(LNode[][] order, boolean forwardSweep) {
        // determine first index (
        int firstIndex = forwardSweep ? 0 : Math.max(0, order.length - 1);
        // extract firstLayer from 2D-array
        List<LNode> firstLayer = Lists.newArrayList(order[firstIndex]);
        // set random weights for nodes in firstLayer
        for (LNode node : firstLayer) {
            node.setProperty(InternalProperties.WEIGHT, random.nextDouble());
        }
        // sort firstLayer by their (randomized) weights
        // Collections.sort() is order-preserving (stable)
        Collections.sort(firstLayer, weightComparator);
        // insert nodes back into array
        int index = 0;
        for (LNode node : firstLayer) {
            order[firstIndex][index++] = node;
            // overwrite the weight with an integer from 1 to n
            node.setProperty(InternalProperties.WEIGHT, (double) index);
            // should the node's id be set as well?
        }

        return false;
    }

    /*
     * (non-Javadoc)
     * 
     * @see
     * org.eclipse.elk.alg.layered.p3order.ICrossingMinimizationHeuristic#minimizeCrossings(org.eclipse.elk.alg.layered.
     * graph.LNode[][], int, boolean, boolean)
     */
    @Override
    public boolean minimizeCrossings(LNode[][] order, int freeLayerIndex, boolean forwardSweep, boolean isFirstSweep) {
        List<LNode> freeLayer = Lists.newArrayList(order[freeLayerIndex]);
        // calculate Medians for the free Layer (does not sort the free layer)
        calculateMedians(freeLayer, forwardSweep ? freeLayerIndex - 1 : freeLayerIndex + 1);
        // sort the free Layer
        Collections.sort(freeLayer, weightComparator);
        // and insert the free Layer back into the array
        int index = 0;
        for (LNode node : freeLayer) {
            order[freeLayerIndex][index++] = node;
            // should the node's id be set as well?
        }
        return false;
    }

    /**
     * calculates the medians and writes them into the nodes as a property
     * 
     * @param nodes
     *            the list of nodes for which medians should be calculated, all in one layer
     * @param referenceLayer
     *            the Layer from which nodes should be taken into account when calculating medians
     */
    private void calculateMedians(List<LNode> nodes, int referenceLayer) {
        // minimum and maximum weight in free layer
        double minWeight = Double.MIN_VALUE;
        double maxWeight = Double.MAX_VALUE;
        // a list to be filled with nodes whose weights cannot be calculated from the referenceLayer alone.
        List<LNode> toRevisit = new ArrayList<>();

        // iterate through nodes
        for (LNode node : nodes) {
            List<LNode> connectedNodes = new ArrayList<>();
            // for every outgoing and incoming edge
            for (LEdge edge : node.getIncomingEdges()) {
                // add the adjacent node's weight to the weight list
                LNode target = edge.getTarget().getNode();
                LNode source = edge.getSource().getNode();
                if (target.getLayer().id == referenceLayer) {
                    connectedNodes.add(target);
                }
                if (source.getLayer().id == referenceLayer) {
                    connectedNodes.add(source);
                }
            }
            for (LEdge edge : node.getOutgoingEdges()) {
                LNode target = edge.getTarget().getNode();
                LNode source = edge.getSource().getNode();
                if (target.getLayer().id == referenceLayer) {
                    connectedNodes.add(target);
                }
                if (source.getLayer().id == referenceLayer) {
                    connectedNodes.add(source);
                }
            }
            // if the node's weight cannot be determined from referenceLayer (no connected nodes in that layer), save
            // this node in order to revisit it later
            if (connectedNodes.isEmpty()) {
                toRevisit.add(node);
            } else {
                // Collections.sort() is stable
                // therefore, no sort of normalization is needed
                // nodes will not keep switching places if they have the same weight
                Collections.sort(connectedNodes, weightComparator);
                // calculate weight from the median of connected nodes' weights
                double newWeight = connectedNodes.get(connectedNodes.size() / 2).getProperty(InternalProperties.WEIGHT);
                node.setProperty(InternalProperties.WEIGHT, newWeight);
                // update minWeight and maxWeight
                minWeight = Math.min(minWeight, newWeight);
                maxWeight = Math.max(maxWeight, newWeight);
            }
        }
        // if no nodes had any weight, avgWeight = 0
        double avgWeight = (maxWeight + minWeight) / 2.0;
        // go through yet unvisited nodes and set their weight to the layer's average
        for (LNode n : toRevisit) {
            n.setProperty(InternalProperties.WEIGHT, avgWeight);
        }
    }

    /*
     * (non-Javadoc)
     * 
     * @see org.eclipse.elk.alg.layered.p3order.ICrossingMinimizationHeuristic#isDeterministic()
     */
    @Override
    public boolean isDeterministic() {
        return true;
    }

    /**
     * Compares two {@link LNode}s based on their weights. Assume that both nodes have weights set.
     */
    protected Comparator<LNode> weightComparator = (n1, n2) -> {
        // even if weights are not set, return an appropriate value
        // if neither weight is set, the return value does not matter (here: return 1)
        if (!n1.hasProperty(InternalProperties.WEIGHT)) {
            return 1;
        }
        if (!n2.hasProperty(InternalProperties.WEIGHT)) {
            return -1;
        }
        Double w1 = n1.getProperty(InternalProperties.WEIGHT);
        Double w2 = n2.getProperty(InternalProperties.WEIGHT);
        if (w1 != null && w2 != null) {
            return w1.compareTo(w2);
        }
        // everything from here should not occur as weightComparator will only be called once weights have been set
        // will it??
        else if (w1 != null) {
            return -1;
        } else if (w2 != null) {
            return 1;
        }
        return 0;
    };

}
