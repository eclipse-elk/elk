/*******************************************************************************
 * Copyright (c) 2016 Kiel University.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.alg.layered.p3order;

import java.util.List;

import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.intermediate.greedyswitch.BetweenLayerEdgeTwoNodeCrossingsCounter;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.p3order.counting.CrossingsCounter;
import org.eclipse.elk.alg.layered.p3order.counting.IInitializable;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.core.options.PortSide;
import org.eclipse.elk.core.util.Pair;
import org.eclipse.elk.alg.layered.graph.LGraphUtil;


/**
 * Distribute ports greedily on a single node.
 * 
 * @author alan
 */
public class GreedyPortDistributor implements ISweepPortDistributor, IInitializable {

    private CrossingsCounter crossingsCounter;
    private int nPorts;
    private int[] portPos;
    private BetweenLayerEdgeTwoNodeCrossingsCounter hierarchicalCrossingsCounter;

    @Override
    public boolean distributePortsWhileSweeping(final LNode[][] nodeOrder, final int currentIndex,
            final boolean isForwardSweep) {
        initialize(nodeOrder, currentIndex, isForwardSweep);

        return distributePortsInLayer(nodeOrder, currentIndex, isForwardSweep);
    }

    private boolean distributePortsInLayer(final LNode[][] nodeOrder, final int currentIndex,
            final boolean isForwardSweep) {
        PortSide side = isForwardSweep ? PortSide.WEST : PortSide.EAST;
        boolean improved = false;
        for (LNode node : nodeOrder[currentIndex]) {
            if (node.getProperty(LayeredOptions.PORT_CONSTRAINTS).isOrderFixed()) {
                continue;
            }
            LGraph nestedGraph = node.getNestedGraph();
            boolean useHierarchicalCrossCounter = !node.getPortSideView(side).isEmpty() && nestedGraph != null;
            if (useHierarchicalCrossCounter) {
                LNode[][] innerGraph = nestedGraph.toNodeArray();
                hierarchicalCrossingsCounter = new BetweenLayerEdgeTwoNodeCrossingsCounter(innerGraph,
                        isForwardSweep ? 0 : innerGraph.length - 1);
            }
            improved |= distributePortsOnNode(node, side, useHierarchicalCrossCounter);
        }
        return improved;
    }

    /**
     * Distribute ports greedily on a single node.
     */
    private boolean distributePortsOnNode(final LNode node, final PortSide side,
            final boolean useHierarchicalCrosscounter) {

        List<LPort> ports = node.getPortSideView(side);
        if (side == PortSide.SOUTH || side == PortSide.WEST) {
            ports = LGraphUtil.reversed(ports);
        }
        boolean improved = false;
        boolean continueSwitching;
        do {
            continueSwitching = false;
            for (int i = 0; i < ports.size() - 1; i++) {
                LPort upperPort = ports.get(i);
                LPort lowerPort = ports.get(i + 1);
                if (switchingDecreasesCrossings(upperPort, lowerPort, node, useHierarchicalCrosscounter)) {
                    improved = true;
                    switchPorts(ports, node, i, i + 1);
                    continueSwitching = true;
                }
            }
        } while (continueSwitching);
        return improved;
    }

    /**
     * Initialize crossings counter for given layers and side.
     */
    private void initForLayers(final LNode[] leftLayer, final LNode[] rightLayer) {
        crossingsCounter.initForCountingBetween(leftLayer, rightLayer);
    }

    private boolean switchingDecreasesCrossings(final LPort upperPort, final LPort lowerPort, final LNode node,
            final boolean useHierarchicalCrosscounter) {
        Pair<Integer, Integer> originalNSwitchedCrossings =
                crossingsCounter.countCrossingsBetweenPortsInBothOrders(upperPort, lowerPort);
        int upperLowerCrossings = originalNSwitchedCrossings.getFirst();
        int lowerUpperCrossings = originalNSwitchedCrossings.getSecond();
        if (useHierarchicalCrosscounter) {
            LNode upperNode = upperPort.getProperty(InternalProperties.PORT_DUMMY);
            LNode lowerNode = lowerPort.getProperty(InternalProperties.PORT_DUMMY);
            if (upperNode != null && lowerNode != null) {
                hierarchicalCrossingsCounter.countBothSideCrossings(upperNode, lowerNode);
                upperLowerCrossings += hierarchicalCrossingsCounter.getUpperLowerCrossings();
                lowerUpperCrossings += hierarchicalCrossingsCounter.getLowerUpperCrossings();
            }
        }
        return upperLowerCrossings > lowerUpperCrossings;
    }

    private void switchPorts(final List<LPort> ports, final LNode node, final int topPort, final int bottomPort) {
        crossingsCounter.switchPorts(ports.get(topPort), ports.get(bottomPort));
        LPort lower = ports.get(bottomPort);
        ports.set(bottomPort, ports.get(topPort));
        ports.set(topPort, lower);
    }

    private void initialize(final LNode[][] nodeOrder, final int currentIndex, final boolean isForwardSweep) {
        if (isForwardSweep && currentIndex > 0) {
            initForLayers(nodeOrder[currentIndex - 1], nodeOrder[currentIndex]);
        } else if (!isForwardSweep && currentIndex < nodeOrder.length - 1) {
            initForLayers(nodeOrder[currentIndex], nodeOrder[currentIndex + 1]);
        } else {
            crossingsCounter.initPortPositionsForInLayerCrossings(nodeOrder[currentIndex],
                    isForwardSweep ? PortSide.WEST : PortSide.EAST);
        }
    }

    @Override
    public void initAtNodeLevel(final int l, final int n, final LNode[][] nodeOrder) {
        LNode node = nodeOrder[l][n];
        nPorts += node.getPorts().size();
    }

    @Override
    public void initAfterTraversal() {
        portPos = new int[nPorts];
        crossingsCounter = new CrossingsCounter(portPos);
    }

}
