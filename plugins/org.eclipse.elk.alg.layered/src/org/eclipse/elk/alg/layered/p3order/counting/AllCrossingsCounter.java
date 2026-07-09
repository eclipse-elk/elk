/*******************************************************************************
 * Copyright (c) 2014, 2015 Kiel University and others.
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
package org.eclipse.elk.alg.layered.p3order.counting;

import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LNode.NodeType;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.core.options.PortSide;

/**
 * Counts all crossings in a graph. Must be initialized using {@link IInitializable#init(java.util.List)}!
 *
 * <p>
 * Must be initialized using {@link IInitializable#init(java.util.List)}!
 * </p>
 */
public final class AllCrossingsCounter implements IInitializable {

    private CrossingsCounter crossingCounter;
    private boolean[] hasHyperEdgesEastOfIndex;
    private HyperedgeCrossingsCounter hyperedgeCrossingsCounter;

    /** For Initialization. */
    private int[] inLayerEdgeCounts;
    private boolean[] hasNorthSouthPorts;
    private int nPorts;

    /**
     * Returns crossings counter. Must be initialized using {@link IInitializable#init(java.util.List)}!
     * 
     * @param graph
     *            the current node order
     */
    public AllCrossingsCounter(final LNode[][] graph) {
        inLayerEdgeCounts = new int[graph.length];
        hasNorthSouthPorts = new boolean[graph.length];
        hasHyperEdgesEastOfIndex = new boolean[graph.length];
        nPorts = 0;
    }

    /**
     * Count all crossings.
     * 
     * @param currentOrder
     *            the current node order
     * @return the number of crossings in the graph
     */
    public int countAllCrossings(final LNode[][] currentOrder) {
        if (currentOrder.length == 0) {
            return 0;
        }
        int crossings = crossingCounter.countInLayerCrossingsOnSide(currentOrder[0], PortSide.WEST);
        crossings += crossingCounter.countInLayerCrossingsOnSide(currentOrder[currentOrder.length - 1], PortSide.EAST);
        for (int layerIndex = 0; layerIndex < currentOrder.length; layerIndex++) {
            crossings += countCrossingsAt(layerIndex, currentOrder);
        }
        return crossings;
    }

    private int countCrossingsAt(final int layerIndex, final LNode[][] currentOrder) {
        int totalCrossings = 0;
        LNode[] leftLayer = currentOrder[layerIndex];
        if (layerIndex < currentOrder.length - 1) {
            LNode[] rightLayer = currentOrder[layerIndex + 1];
            if (hasHyperEdgesEastOfIndex[layerIndex]) {
                totalCrossings = hyperedgeCrossingsCounter.countCrossings(leftLayer, rightLayer);
                totalCrossings += crossingCounter.countInLayerCrossingsOnSide(leftLayer, PortSide.EAST);
                totalCrossings += crossingCounter.countInLayerCrossingsOnSide(rightLayer, PortSide.WEST);
            } else {
                totalCrossings = crossingCounter.countCrossingsBetweenLayers(leftLayer, rightLayer);
            }
        }

        if (hasNorthSouthPorts[layerIndex]) {
            totalCrossings += crossingCounter.countNorthSouthPortCrossingsInLayer(leftLayer);
        }
        
        return totalCrossings;
    }

    @Override
    public void initAtNodeLevel(final int l, final int n, final LNode[][] nodeOrder) {
        LNode node = nodeOrder[l][n];
        hasNorthSouthPorts[l] |= node.getType() == NodeType.NORTH_SOUTH_PORT;
    }

    @Override
    public void initAtPortLevel(final int l, final int n, final int p, final LNode[][] nodeOrder) {
        LPort port = nodeOrder[l][n].getPorts().get(p);
        port.id = nPorts++;
        if (port.getOutgoingEdges().size() + port.getIncomingEdges().size() > 1) {
            if (port.getSide() == PortSide.EAST) {
                hasHyperEdgesEastOfIndex[l] = true;
            } else if (port.getSide() == PortSide.WEST && l > 0) {
                hasHyperEdgesEastOfIndex[l - 1] = true;
            }
        }
    }

    @Override
    public void initAtEdgeLevel(final int l, final int n, final int p, final int e, final LEdge edge,
            final LNode[][] nodeOrder) {
        LPort port = nodeOrder[l][n].getPorts().get(p);
        if (edge.getSource() == port
                && edge.getSource().getNode().getLayer() == edge.getTarget().getNode().getLayer()) {
            inLayerEdgeCounts[l]++;
        }
    }

    @Override
    public void initAfterTraversal() {
        int[] portPos = new int[nPorts];
        hyperedgeCrossingsCounter = new HyperedgeCrossingsCounter(inLayerEdgeCounts, hasNorthSouthPorts, portPos);
        crossingCounter = new CrossingsCounter(portPos);
    }
}
