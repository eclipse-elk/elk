/*******************************************************************************
 * Copyright (c) 2010, 2015 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.alg.layered.p3order;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.ListIterator;
import java.util.Random;

import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LNode.NodeType;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.layered.options.PortType;

import com.google.common.collect.Lists;

/**
 * Determines the node order of a given free layer. Uses heuristic methods to find an ordering that
 * minimizes edge crossings between the given free layer and a neighboring layer with fixed node
 * order. The barycenter heuristic is used here.
 * 
 * @author msp
 * @author cds
 * @author ima
 */
public class BarycenterHeuristic implements ICrossingMinimizationHeuristic {

    /** the array of port ranks. */
    protected float[] portRanks;
    /** the random number generator. */
    protected final Random random;
    /** the constraint resolver for ordering constraints. */
    protected ForsterConstraintResolver constraintResolver;
    /** the barycenter values of every node in the graph, indexed by layer.id and node.id. */
    protected BarycenterState[][] barycenterState;
    /** The Barycenter PortDistributor is used to ask for the port ranks.*/
    protected final AbstractBarycenterPortDistributor portDistributor;

    /**
     * Don't use!
     * Only public to be accessible by a test.
     */
    public void minimizeCrossings(final List<LNode> layer, final boolean preOrdered,
            final boolean randomize, final boolean forward) {

        if (randomize) {
            // Randomize barycenters (we don't need to update the edge count in this case;
            // there are no edges of interest since we're only concerned with one layer)
            // simply a permutation of nodes in layer
            randomizeBarycenters(layer);
        } else {
            // Calculate barycenters and assign barycenters to barycenterless node groups
            calculateBarycenters(layer, forward);
            fillInUnknownBarycenters(layer, preOrdered);
        }

        if (layer.size() > 1) {
            // Sort the vertices according to their barycenters
            if (layer.get(0).getGraph().getProperty(LayeredOptions.CROSSING_MINIMIZATION_FORCE_NODE_MODEL_ORDER)) {
                ModelOrderBarycenterHeuristic.insertionSort(layer, barycenterStateComparator,
                        (ModelOrderBarycenterHeuristic) this);
            } else {
                // use comparator based on "<"
                Collections.sort(layer, barycenterStateComparator);
            }

            // Resolve ordering constraints
            constraintResolver.processConstraints(layer);
        }
    }

    /**
     * Randomize the order of nodes for the given layer.
     * 
     * @param nodes
     *            a layer
     */
    protected void randomizeBarycenters(final List<LNode> nodes) {
        for (LNode node : nodes) {
            // Set barycenters only for nodeGroups containing a single node.
            stateOf(node).barycenter = random.nextDouble();
            stateOf(node).summedWeight = stateOf(node).barycenter;
            stateOf(node).degree = 1;
        }
    }

    /**
     * Try to find appropriate barycenter values for node groups whose barycenter is undefined.
     * 
     * @param nodes
     *            the nodeGroups to fill in barycenters for.
     * @param preOrdered
     *            whether the nodeGroups have been ordered in a previous run.
     */
    protected void fillInUnknownBarycenters(final List<LNode> nodes, final boolean preOrdered) {
        // Determine placements for nodes with undefined barycenter value
        if (preOrdered) {
            double lastValue = -1;
            ListIterator<LNode> nodesIterator = nodes.listIterator();

            while (nodesIterator.hasNext()) {
                LNode node = nodesIterator.next();
                Double value = stateOf(node).barycenter;

                if (value == null) {
                    // The barycenter is undefined - take the center of the previous value
                    // and the next defined value in the list
                    double nextValue = lastValue + 1;

                    ListIterator<LNode> nextNodeIterator =
                            nodes.listIterator(nodesIterator.nextIndex());
                    while (nextNodeIterator.hasNext()) {
                        Double x = stateOf(nextNodeIterator.next()).barycenter;
                        if (x != null) {
                            nextValue = x;
                            break;
                        }
                    }

                    value = (lastValue + nextValue) / 2;
                    stateOf(node).barycenter = value;
                    stateOf(node).summedWeight = value;
                    stateOf(node).degree = 1;
                }

                lastValue = value;
            }
        } else {
            // No previous ordering - determine random placement for new nodes
            double maxBary = 0;
            for (LNode node : nodes) {
                if (stateOf(node).barycenter != null) {
                    maxBary = Math.max(maxBary, stateOf(node).barycenter);
                }
            }

            maxBary += 2;
            for (LNode node : nodes) {
                if (stateOf(node).barycenter == null) {
                    double value = random.nextFloat() * maxBary - 1;
                    stateOf(node).barycenter = value;
                    stateOf(node).summedWeight = value;
                    stateOf(node).degree = 1;
                }
            }
        }
    }
    
    /**
     * Calculate the barycenters of the given node groups.
     * 
     * @param nodes
     *            the nodes
     * @param forward
     *            {@code true} if the current sweep moves forward
     */
    protected void calculateBarycenters(final List<LNode> nodes, final boolean forward) {
        // Set all visited flags to false
        for (LNode node : nodes) {
            stateOf(node).visited = false;
        }

        for (LNode node : nodes) {
            // Calculate the node groups's new barycenter (may be null)
            calculateBarycenter(node, forward);
        }
    }

    /**
     * Constructs a Barycenter heuristic for crossing minimization.
     * 
     * @param constraintResolver
     *            the constraint resolver
     * @param random
     *            the random number generator
     * @param portDistributor
     *            calculates the port ranks for the barycenter heuristic.
     * @param graph
     *            current node order
     */
    public BarycenterHeuristic(final ForsterConstraintResolver constraintResolver, final Random random,
            final AbstractBarycenterPortDistributor portDistributor, final LNode[][] graph) {
        this.constraintResolver = constraintResolver;
        this.random = random;
        this.portDistributor = portDistributor;
    }

    /** the amount of random value to add to each calculated barycenter. */
    private static final float RANDOM_AMOUNT = 0.07f;

    /**
     * Calculate the barycenter of the given single-node node group. This method is able to handle
     * in-layer edges, but it may give incorrect results if the in-layer edges form a cycle.
     * However, such cases do not occur in the present implementation.
     * 
     * @param node
     *            a node group consisting of a single node
     * @param forward
     *            {@code true} if the current sweep moves forward
     * @param portPos
     *            position array
     * @return a pair containing the summed port positions of the connected ports as the first, and
     *         the number of connected edges as the second entry.
     */
    private void calculateBarycenter(final LNode node, final boolean forward) {// javadoc for this method outdated

        // Check if the node group's barycenter was already computed
        if (stateOf(node).visited) {
            return;
        } else {
            stateOf(node).visited = true;
        }

        stateOf(node).degree = 0;
        stateOf(node).summedWeight = 0.0f;
        stateOf(node).barycenter = null;

        for (LPort freePort : node.getPorts()) {
            Iterable<LPort> portIterable =
                    forward ? freePort.getPredecessorPorts() : freePort.getSuccessorPorts();
            for (LPort fixedPort : portIterable) {
                // If the node the fixed port belongs to is part of the free layer (thus, if
                // we have an in-layer edge), use that node's barycenter calculation instead
                LNode fixedNode = fixedPort.getNode();

                if (fixedNode.getLayer() == node.getLayer()) {
                    // Self-loops are ignored
                    if (fixedNode != node) {
                        // Find the fixed node's node group and calculate its barycenter
                        calculateBarycenter(fixedNode, forward);

                        // Update this node group's values
                        stateOf(node).degree += stateOf(fixedNode).degree;
                        stateOf(node).summedWeight += stateOf(fixedNode).summedWeight;
                    } // else { }
                } else {
                    stateOf(node).summedWeight += portRanks[fixedPort.id];
                    stateOf(node).degree++;
                }
            }
        }

        // Iterate over the node's barycenter associates
        List<LNode> barycenterAssociates =
                node.getProperty(InternalProperties.BARYCENTER_ASSOCIATES);
        if (barycenterAssociates != null) {
            for (LNode associate : barycenterAssociates) {
                // Make sure the associate is in the same layer as this node
                if (node.getLayer() == associate.getLayer()) {
                    // Find the associate's node group and calculate its barycenter
                    calculateBarycenter(associate, forward);

                    // Update this vertex's values
                    stateOf(node).degree += stateOf(associate).degree;
                    stateOf(node).summedWeight += stateOf(associate).summedWeight;
                }
            }
        }

        if (stateOf(node).degree > 0) {
            // add a small random perturbation in order to increase diversity of solutions
            stateOf(node).summedWeight += random.nextFloat() * RANDOM_AMOUNT - RANDOM_AMOUNT / 2;
            stateOf(node).barycenter = stateOf(node).summedWeight / stateOf(node).degree;
        }
    }

    protected BarycenterState stateOf(final LNode node) {
        return barycenterState[node.getLayer().id][node.id];
    }
    
    /**
     * Represents the current barycenter state of a node.
     */
    public static final class BarycenterState {

        // SUPPRESS CHECKSTYLE NEXT 20 VisibilityModifier
        /** The node this state holds data for. */
        public LNode node;
        /** The sum of the node weights. Each node weight is the sum of the weights of the ports the
         *  node's ports are connected to. */
        public double summedWeight;
        /** The number of ports relevant to the barycenter calculation. */
        public int degree;
        /** This vertex' barycenter value. (summedWeight / degree) */
        public Double barycenter;
        /** Whether the node group has been visited in some traversing algorithm. */
        public boolean visited;

        /**
         * @param node
         *            the node this state holds data for.
         */
        public BarycenterState(final LNode node) {
            this.node = node;
        }

        @Override
        public String toString() {
            return "BarycenterState [node=" + node + ", summedWeight=" + summedWeight + ", degree=" + degree
                    + ", barycenter=" + barycenter + ", visited=" + visited + "]";
        }

    }
    
    /**
     * Compares two {@link LNode}s based on their barycenter values.
     */
    protected Comparator<LNode> barycenterStateComparator = 
        (n1, n2) -> {
            BarycenterState s1 = stateOf(n1);
            BarycenterState s2 = stateOf(n2);
                if (s1.barycenter != null && s2.barycenter != null) {
                    return s1.barycenter.compareTo(s2.barycenter);
                } else if (s1.barycenter != null) {
                    return -1;
                } else if (s2.barycenter != null) {
                    return 1;
                }
                return 0;
        };

    @Override
    public boolean minimizeCrossings(final LNode[][] order, final int freeLayerIndex,
            final boolean forwardSweep, final boolean isFirstSweep) {
        if (!isFirstLayer(order, freeLayerIndex, forwardSweep)) {
            LNode[] fixedLayer = order[freeLayerIndex - changeIndex(forwardSweep)];
            portDistributor.calculatePortRanks(fixedLayer, portTypeFor(forwardSweep));
        }

        LNode firstNodeInLayer = order[freeLayerIndex][0];
        boolean preOrdered = !isFirstSweep || isExternalPortDummy(firstNodeInLayer);

        List<LNode> nodes = Lists.newArrayList(order[freeLayerIndex]);
        minimizeCrossings(nodes, preOrdered, false, forwardSweep);
        // apply the new ordering
        int index = 0;
        for (LNode nodeGroup : nodes) {
                order[freeLayerIndex][index++] = nodeGroup;
        }

        return false; // Does not always improve.
    }

    @Override
    public boolean setFirstLayerOrder(final LNode[][] order, final boolean isForwardSweep) {
        // if sweeping forward, startIndex = 0, else the last existing element of order
        int startIndex = startIndex(isForwardSweep, order.length);
        // extract first layer into List<LNode>
        List<LNode> nodes = Lists.newArrayList(
                order[startIndex]);
        // randomize nodes' barycenters
        minimizeCrossings(nodes, false, true, isForwardSweep);
        // fill first layer with nodes
        int index = 0;
        for (LNode nodeGroup : nodes) {
                order[startIndex][index++] = nodeGroup;
        }

        return false; // Does not always improve
    }

    private boolean isExternalPortDummy(final LNode firstNode) {
        return firstNode.getType() == NodeType.EXTERNAL_PORT;
    }

    private int changeIndex(final boolean dir) {
        return dir ? 1 : -1;
    }

    private PortType portTypeFor(final boolean direction) {
        return direction ? PortType.OUTPUT : PortType.INPUT;
    }

    private int startIndex(final boolean dir, final int length) {
        return dir ? 0 : Math.max(0, length - 1);
    }

    private boolean isFirstLayer(final LNode[][] nodeOrder, final int currentIndex,
            final boolean forwardSweep) {
        return currentIndex == startIndex(forwardSweep, nodeOrder.length);
    }

    @Override
    public boolean alwaysImproves() {
        return false;
    }

    @Override
    public boolean isDeterministic() {
        return false;
    }

    @Override
    public void initAfterTraversal() {
        barycenterState = constraintResolver.getBarycenterStates();
        portRanks = portDistributor.getPortRanks();
    }

    @Override
    public void initAtLayerLevel(final int l, final LNode[][] nodeOrder) {
        nodeOrder[l][0].getLayer().id = l;
    }
}
