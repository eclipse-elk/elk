/*******************************************************************************
 * Copyright (c) 2014, 2015 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.alg.layered.intermediate.greedyswitch;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.ArrayList;
import java.util.HashMap;

import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.p3order.counting.CrossMinUtil;
import org.eclipse.elk.core.options.PortSide;


/**
 * Calculates the number of crossings for edges incident to two nodes. In the case where there is
 * free port order and two edges go into one port, this crossing counter can in some cases count to
 * few crossings. See the ignored test in the test class.
 */
public final class BetweenLayerEdgeTwoNodeCrossingsCounter {
    private int upperLowerCrossings;
    private int lowerUpperCrossings;
    private AdjacencyList upperAdjacencies;
    private AdjacencyList lowerAdjacencies;
    private final LNode[][] currentNodeOrder;
    private final int freeLayerIndex;
    private final Map<LPort, Integer> portPositions;
    private final Map<LNode, AdjacencyList> easternAdjacencies;
    private final Map<LNode, AdjacencyList> westernAdjacencies;

    /**
     * Create {@link BetweenLayerEdgeTwoNodeCrossingsCounter}. Naming assumes a left-right layer
     * ordering.
     * 
     * @param currentNodeOrder
     *            Currently considered node ordering.
     * @param freeLayerIndex
     *            Index of free layer.
     */
    public BetweenLayerEdgeTwoNodeCrossingsCounter(final LNode[][] currentNodeOrder, final int freeLayerIndex) {
        portPositions = new HashMap<>();
        easternAdjacencies = new HashMap<>();
        westernAdjacencies = new HashMap<>();
        this.currentNodeOrder = currentNodeOrder;
        this.freeLayerIndex = freeLayerIndex;
        setPortPositionsForNeighbouringLayers();
    }

    private void setPortPositionsForNeighbouringLayers() {
        if (freeLayerIsNotFirstLayer()) {
            setPortPositionsForLayer(freeLayerIndex - 1, PortSide.EAST);
        }
        if (freeLayerIsNotLastLayer()) {
            setPortPositionsForLayer(freeLayerIndex + 1, PortSide.WEST);
        }
    }

    private boolean freeLayerIsNotFirstLayer() {
        return freeLayerIndex > 0;
    }

    private boolean freeLayerIsNotLastLayer() {
        return freeLayerIndex < currentNodeOrder.length - 1;
    }

    private void setPortPositionsForLayer(final int layerIndex, final PortSide portSide) {
        int portId = 0;
        for (LNode node : currentNodeOrder[layerIndex]) {
            Iterable<LPort> ports = CrossMinUtil.inNorthSouthEastWestOrder(node, portSide);
            for (LPort port : ports) {
                portPositions.put(port, portId++);
            }
        }
    }

    /**
     * Calculates the number of crossings for incident edges coming from the west to the nodes. The
     * crossing numbers can be received with getCrossingForOrderUpperLower and
     * getCrossingForOrderLowerUpper for the order upperNode -> lowerNode or lowerNode -> upperNode
     * respectively.
     * 
     * @param upperNode
     *            Upper node assuming left-right layout.
     * @param lowerNode
     *            Lower node assuming left-right layout
     */
    public void countEasternEdgeCrossings(final LNode upperNode, final LNode lowerNode) {
        resetCrossingCount();
        if (upperNode.equals(lowerNode)) {
            return;
        }
        addEasternCrossings(upperNode, lowerNode);
    }

    /**
     * Calculates the number of crossings for incident edges coming from the east to the nodes. The
     * crossing numbers can be received with getCrossingForOrderUpperLower and
     * getCrossingForOrderLowerUpper for the order upperNode -> lowerNode or lowerNode -> upperNode
     * respectively.
     * 
     * @param upperNode
     *            Upper node assuming left-right layout.
     * @param lowerNode
     *            Lower node assuming left-right layout
     */
    public void countWesternEdgeCrossings(final LNode upperNode, final LNode lowerNode) {
        resetCrossingCount();
        if (upperNode.equals(lowerNode)) {
            return;
        }
        addWesternCrossings(upperNode, lowerNode);
    }

    /**
     * Calculates the number of crossings for incident edges coming from the both sides to the
     * nodes. The crossing numbers can be received with getCrossingForOrderUpperLower and
     * getCrossingForOrderLowerUpper for the order upperNode -> lowerNode or lowerNode -> upperNode
     * respectively.
     * 
     * @param upperNode
     *            Upper node assuming left-right layout.
     * @param lowerNode
     *            Lower node assuming left-right layout
     */
    public void countBothSideCrossings(final LNode upperNode, final LNode lowerNode) {
        resetCrossingCount();
        if (upperNode.equals(lowerNode)) {
            return;
        }
        addWesternCrossings(upperNode, lowerNode);
        addEasternCrossings(upperNode, lowerNode);
    }

    private void resetCrossingCount() {
        upperLowerCrossings = 0;
        lowerUpperCrossings = 0;
    }

    private void addEasternCrossings(final LNode upperNode, final LNode lowerNode) {
        upperAdjacencies = getAdjacencyFor(upperNode, PortSide.EAST, easternAdjacencies);
        lowerAdjacencies = getAdjacencyFor(lowerNode, PortSide.EAST, easternAdjacencies);
        if (upperAdjacencies.size() == 0 || lowerAdjacencies.size() == 0) {
            return;
        }
        countCrossingsByMergingAdjacencyLists();
    }

    /**
     * Since calculating adjacencies is a little expensive, it is only done once for each
     * configuration and the sorted adjacencies saved in Maps.
     */
    private AdjacencyList getAdjacencyFor(final LNode node, final PortSide side,
            final Map<LNode, AdjacencyList> adjacencies) {
        
        if (adjacencies.isEmpty()) {
            for (LNode n : currentNodeOrder[freeLayerIndex]) {
                adjacencies.put(n, new AdjacencyList(n, side));
            }
        }
        AdjacencyList aL = adjacencies.get(node);
        aL.reset();
        return aL;
    }

    private void addWesternCrossings(final LNode upperNode, final LNode lowerNode) {
        upperAdjacencies = getAdjacencyFor(upperNode, PortSide.WEST, westernAdjacencies);
        lowerAdjacencies = getAdjacencyFor(lowerNode, PortSide.WEST, westernAdjacencies);
        if (upperAdjacencies.size() == 0 || lowerAdjacencies.size() == 0) {
            return;
        }
        countCrossingsByMergingAdjacencyLists();
    }

    /**
     * The main algorithm. Adjacency Lists are lists of ports connected to a node. If a connected
     * node has no fixed port ordering all ports have the same position value. By merging adjacency
     * lists, both the number of between-layer crossings for the order upper - lower and for the
     * opposite order can be found. Consider:
     * 
     * <pre>
     * A   p0 
     *  \\/
     *  /\\
     * B   p1
     * </pre>
     * <ul>
     * <li>The adjacency list La for A is: p1, p1
     * <li>The adjacency list Lb for B is: p0
     * </ul>
     * Since p1 is below p0, edge (B, p0) is crossed by all edges from A and can p0 can be removed.
     * If the next adjacencies in both lists have the same position value p we add to
     * upperLowerCrossings the number of remaining adjacencies upperAdjacencies below the current
     * node. We do the same for lowerUpperCrossings.
     */
    private void countCrossingsByMergingAdjacencyLists() {
        while (!upperAdjacencies.isEmpty() && !lowerAdjacencies.isEmpty()) {
            if (isBelow(upperAdjacencies.first(), lowerAdjacencies.first())) {
                upperLowerCrossings += upperAdjacencies.size();
                lowerAdjacencies.removeFirst();
            } else if (isBelow(lowerAdjacencies.first(), upperAdjacencies.first())) {
                lowerUpperCrossings += lowerAdjacencies.size();
                upperAdjacencies.removeFirst();
            } else {
                upperLowerCrossings += upperAdjacencies.countAdjacenciesBelowNodeOfFirstPort();
                lowerUpperCrossings += lowerAdjacencies.countAdjacenciesBelowNodeOfFirstPort();
                upperAdjacencies.removeFirst();
                lowerAdjacencies.removeFirst();
            }
        }
    }

    private boolean isBelow(final int firstPort, final int secondPort) {
        return firstPort > secondPort;
    }

    /**
     * The adjacency list of a node holds the position of connected ports in a neighboring layer on
     * the given side. Since we want to save it for further use, the remove operation does not
     * actually delete the entries in the adjacency list. Instead we use currentIndex, currentSize
     * and currentCardinality (in the inner class) to show the current state of the list. Use
     * reset() to reset to the original state.
     */
    private class AdjacencyList {
        private final LNode node;
        private final List<Adjacency> adjacencyList;
        private final PortSide side;
        private int size;
        private int currentSize;
        private int currentIndex;
        
        AdjacencyList(final LNode node, final PortSide side) {
            this.node = node;
            this.side = side;
            adjacencyList = new ArrayList<>();
            getAdjacenciesSortedByPosition();
        }

        private void getAdjacenciesSortedByPosition() {
            iterateTroughEdgesCollectingAdjacencies();

            Collections.sort(adjacencyList);
        }

        private void iterateTroughEdgesCollectingAdjacencies() {
            Iterable<LPort> ports = CrossMinUtil.inNorthSouthEastWestOrder(node, side);
            for (LPort port : ports) {
                List<LEdge> edges = getEdgesConnectedTo(port);
                for (LEdge edge : edges) {
                    if (!edge.isSelfLoop() && isNotInLayer(edge)) {
                        addAdjacencyOf(edge);
                        size++;
                        currentSize++;
                    }
                }
            }
        }

        private List<LEdge> getEdgesConnectedTo(final LPort port) {
            return side == PortSide.WEST ? port.getIncomingEdges() : port.getOutgoingEdges();
        }

        private boolean isNotInLayer(final LEdge edge) {
            return edge.getSource().getNode().getLayer() != edge.getTarget().getNode().getLayer();
        }

        private void addAdjacencyOf(final LEdge edge) {
            LPort adjacentPort = adjacentPortOf(edge, side);
            int adjacentPortPosition = portPositions.get(adjacentPort);
            int lastIndex = adjacencyList.size() - 1;
            if (!adjacencyList.isEmpty()
                    && adjacencyList.get(lastIndex).position == adjacentPortPosition) {
                adjacencyList.get(lastIndex).cardinality++;
                adjacencyList.get(lastIndex).currentCardinality++;
            } else {
                adjacencyList.add(new Adjacency(adjacentPortPosition, adjacentPort));
            }
        }

        private LPort adjacentPortOf(final LEdge e, final PortSide s) {
            return s == PortSide.WEST ? e.getSource() : e.getTarget();
        }

        public void reset() {
            currentIndex = 0;
            currentSize = size;
            if (!isEmpty()) {
                currentAdjacency().reset();
            }
        }

        public int countAdjacenciesBelowNodeOfFirstPort() {
            return currentSize - currentAdjacency().currentCardinality;
        }

        public void removeFirst() {
            if (isEmpty()) {
                return;
            }
            Adjacency currentEntry = currentAdjacency();
            if (currentEntry.currentCardinality == 1) {
                incrementCurrentIndex();
            } else {
                currentEntry.currentCardinality--;
            }

            currentSize--;
        }

        private void incrementCurrentIndex() {
            currentIndex++;
            // reset Adjacency for reuse
            if (currentIndex < adjacencyList.size()) {
                currentAdjacency().reset();
            }
        }

        public boolean isEmpty() {
            return currentSize == 0;
        }

        public int first() {
            return currentAdjacency().position;
        }

        public int size() {
            return currentSize;
        }

        private Adjacency currentAdjacency() {
            return adjacencyList.get(currentIndex);
        }

        @Override
        public String toString() {
            return "AdjacencyList [node=" + node + ", adjacencies= " + adjacencyList + "]";
        }

        /**
         * Adjacency containing only the position and number of ports with the same position.
         */
        private class Adjacency implements Comparable<Adjacency> {
            /** The position of the port. */
            private final int position;
            /** The number of adjacencies with the same position. */
            private int cardinality;
            /** The current number of adjacencies with the same position. */
            private int currentCardinality;

            Adjacency(final int adjacentPortPosition, final LPort port) {
                position = adjacentPortPosition;
                cardinality = 1;
                currentCardinality = 1;
            }

            public void reset() {
                currentCardinality = cardinality;
            }

            @Override
            public int compareTo(final Adjacency o) {
                return (position < o.position) ? -1 : ((position == o.position) ? 0 : 1);
            }

            @Override
            public String toString() {
                return "Adjacency [position=" + position + ", cardinality=" + cardinality
                        + ", currentCardinality=" + currentCardinality + "]";
            }
        }
    }

    /**
     * @return the upperLowerCrossings
     */
    public int getUpperLowerCrossings() {
        return upperLowerCrossings;
    }

    /**
     * @return the lowerUpperCrossings
     */
    public int getLowerUpperCrossings() {
        return lowerUpperCrossings;
    }
}
