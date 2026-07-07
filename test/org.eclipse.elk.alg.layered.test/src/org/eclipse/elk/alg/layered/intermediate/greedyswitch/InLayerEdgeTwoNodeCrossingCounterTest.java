/*******************************************************************************
 * Copyright (c) 2016 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.layered.intermediate.greedyswitch;

import static org.hamcrest.CoreMatchers.is;
import static org.junit.Assert.assertThat;

import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.p3order.counting.CrossingsCounter;
import org.eclipse.elk.core.options.PortSide;
import org.eclipse.elk.core.util.Pair;
import org.junit.Ignore;
import org.junit.Test;

/**
 * Tests {@link InLayerEdgeTwoNodeCrossingCounter}.
 *
 * @author alan
 */
public class InLayerEdgeTwoNodeCrossingCounterTest extends InLayerEdgeTestGraphCreator {
    private CrossingsCounter leftCounter;
    private CrossingsCounter rightCounter;
    private int lowerUpperCrossings;
    private int upperLowerCrossings;
    private LNode[] nodeOrder;

    //CHECKSTYLEOFF javadoc
    //CHECKSTYLEOFF MagicNumber
    //CHECKSTYLEOFF MethodName

    @Test
    public void ignoresInBetweenLayerEdges() {
        getCrossFormedGraph();
        countCrossingsInLayerForUpperNodeLowerNode(1, 0, 1);

        assertThat("upperLowerCrossings", upperLowerCrossings, is(0));
        assertThat("lowerUpperCrossings", lowerUpperCrossings, is(0));
    }

    @Test
    public void countInLayerEdgeWithNormalEdgeCrossing() {
        getInLayerEdgesGraph();
        countCrossingsInLayerForUpperNodeLowerNode(1, 0, 1);

        assertThat("upperLowerCrossings", upperLowerCrossings, is(1));
        assertThat("lowerUpperCrossings", lowerUpperCrossings, is(0));
    }

    @Test
    public void crossingsWhenSwitched() {
        getInLayerEdgesGraphWhichResultsInCrossingsWhenSwitched();

        countCrossingsInLayerForUpperNodeLowerNode(1, 1, 2);

        assertThat("upperLowerCrossings", upperLowerCrossings, is(0));
        assertThat("lowerUpperCrossings", lowerUpperCrossings, is(1));
    }

    @Test
    public void inLayerEdgeOnLowerNode() {
        getInLayerEdgesGraph();

        countCrossingsInLayerForUpperNodeLowerNode(1, 0, 1);

        assertThat("upperLowerCrossings", upperLowerCrossings, is(1));
        assertThat("lowerUpperCrossings", lowerUpperCrossings, is(0));
    }

    @Test
    public void switchNodeOrder() {
        getInLayerEdgesGraph();

        initCrossingCounterForLayerIndex(1);
        switchOrderAndNotifyCounter(1, 2);

        countCrossings(0, 1);

        assertThat("upperLowerCrossings", upperLowerCrossings, is(0));
        assertThat("lowerUpperCrossings", lowerUpperCrossings, is(0));
    }

    @Test
    public void fixedPortOrderCrossingToInBetweenLayerEdge() {
        getInLayerEdgesGraphWithCrossingsToBetweenLayerEdgeWithFixedPortOrder();
        countCrossingsInLayerForUpperNodeLowerNode(1, 0, 1);
        assertThat("upperLowerCrossings", upperLowerCrossings, is(1));
        assertThat("lowerUpperCrossings", lowerUpperCrossings, is(2));
        switchOrderAndNotifyCounter(0, 1);

        countCrossings(0, 1);

        assertThat("upperLowerCrossings", upperLowerCrossings, is(2));
        assertThat("lowerUpperCrossings", lowerUpperCrossings, is(1));
    }

    @Test
    public void fixedPortOrderCrossingsAndNormalEdgeCrossings() {
        getInLayerEdgesWithFixedPortOrderAndNormalEdgeCrossings();

        countCrossingsInLayerForUpperNodeLowerNode(1, 0, 1);

        assertThat("upperLowerCrossings", upperLowerCrossings, is(2));
        assertThat("lowerUpperCrossings", lowerUpperCrossings, is(1));

        switchOrderAndNotifyCounter(0, 1);

        countCrossings(0, 1);

        assertThat("upperLowerCrossings", upperLowerCrossings, is(1));
        assertThat("lowerUpperCrossings", lowerUpperCrossings, is(2));
    }

    @Test
    public void ignoresSelfLoops() {
        getCrossWithManySelfLoopsGraph();

        countCrossingsInLayerForUpperNodeLowerNode(1, 0, 1);

        assertThat("upperLowerCrossings", upperLowerCrossings, is(0));
        assertThat("lowerUpperCrossings", lowerUpperCrossings, is(0));
    }


    @Test
    public void crossingsOnBothSides() {
        getInLayerCrossingsOnBothSides();

        countCrossingsInLayerForUpperNodeLowerNode(1, 0, 1);

        assertThat("upperLowerCrossings", upperLowerCrossings, is(2));
        assertThat("lowerUpperCrossings", lowerUpperCrossings, is(0));
    }

    @Test
    public void fixedPortOrderInLayerNoCrossings() {
        getFixedPortOrderInLayerEdgesDontCrossEachOther();

        countCrossingsInLayerForUpperNodeLowerNode(0, 0, 1);

        assertThat("upperLowerCrossings", upperLowerCrossings, is(0));
        assertThat("lowerUpperCrossings", lowerUpperCrossings, is(0));
    }

    @Test
    public void fixedPortOrderInLayerWithAlwaysRemaingCrossingsAreNotCounted() {
        getFixedPortOrderInLayerEdgesWithCrossings();

        countCrossingsInLayerForUpperNodeLowerNode(0, 0, 1);

        assertThat("upperLowerCrossings", upperLowerCrossings, is(1));
        assertThat("lowerUpperCrossings", lowerUpperCrossings, is(1));
    }

    @Test
    public void oneNode() {
        getOneNodeGraph();

        countCrossingsInLayerForUpperNodeLowerNode(0, 0, 0);

        assertThat("upperLowerCrossings", upperLowerCrossings, is(0));
        assertThat("lowerUpperCrossings", lowerUpperCrossings, is(0));
    }

    @Test
    public void moreComplex() {
        getMoreComplexInLayerGraph();

        countCrossingsInLayerForUpperNodeLowerNode(1, 0, 1);

        assertThat("upperLowerCrossings", upperLowerCrossings, is(6));
        assertThat("lowerUpperCrossings", lowerUpperCrossings, is(6));
    }

    @Test
    public void downwardInLayerEdgesOnLowerNode() {
        getInLayerEdgesFixedPortOrderInLayerAndInBetweenLayerCrossing();
        countCrossingsInLayerForUpperNodeLowerNode(1, 0, 1);

        assertThat("upperLowerCrossings", upperLowerCrossings, is(2));
        assertThat("lowerUpperCrossings", lowerUpperCrossings, is(2));
    }

    @Test
    public void oneLayerInLayerCrossingShouldDisappearAfterAnySwitch() {
        getOneLayerWithInLayerCrossings();

        countCrossingsInLayerForUpperNodeLowerNode(0, 0, 1);

        assertThat("upperLowerCrossings", upperLowerCrossings, is(1));
        assertThat("lowerUpperCrossings", lowerUpperCrossings, is(0));

        countCrossingsInLayerForUpperNodeLowerNode(0, 1, 2);

        assertThat("upperLowerCrossings", upperLowerCrossings, is(1));
        assertThat("lowerUpperCrossings", lowerUpperCrossings, is(0));

        countCrossingsInLayerForUpperNodeLowerNode(0, 2, 3);

        assertThat("upperLowerCrossings", upperLowerCrossings, is(1));
        assertThat("lowerUpperCrossings", lowerUpperCrossings, is(0));

        switchOrderAndNotifyCounter(0, 1);

        countCrossings(0, 1);

        assertThat("upperLowerCrossings", upperLowerCrossings, is(0));
        assertThat("lowerUpperCrossings", lowerUpperCrossings, is(1));

        switchOrderAndNotifyCounter(0, 1);
        switchOrderAndNotifyCounter(1, 2);

        countCrossings(1, 2);

        assertThat("upperLowerCrossings", upperLowerCrossings, is(0));
        assertThat("lowerUpperCrossings", lowerUpperCrossings, is(1));

        switchOrderAndNotifyCounter(1, 2);
        switchOrderAndNotifyCounter(2, 3);

        countCrossings(2, 3);

        assertThat("upperLowerCrossings", upperLowerCrossings, is(0));
        assertThat("lowerUpperCrossings", lowerUpperCrossings, is(1));
    }

    @Test
    public void moreThanOneEdgeIntoAPort() {
        getInLayerEdgesMultipleEdgesIntoSinglePort();

        countCrossingsInLayerForUpperNodeLowerNode(1, 1, 2);

        assertThat("upperLowerCrossings", upperLowerCrossings, is(2));
        assertThat("lowerUpperCrossings", lowerUpperCrossings, is(0));
    }

    @Test
    public void inBetweenLayerEdgesIntoNodeWithNoFixedPortOrderCauseCrossings() {
        multipleInBetweenLayerEdgesIntoNodeWithNoFixedPortOrderCauseCrossings();

        countCrossingsInLayerForUpperNodeLowerNode(1, 0, 1);

        assertThat("upperLowerCrossings", upperLowerCrossings, is(2));
        assertThat("lowerUpperCrossings", lowerUpperCrossings, is(0));

        countCrossingsInLayerForUpperNodeLowerNode(1, 1, 2);

        assertThat("upperLowerCrossings", upperLowerCrossings, is(2));
        assertThat("lowerUpperCrossings", lowerUpperCrossings, is(0));
    }

    @Test
    public void inLayerEdgesPassEachOther() {
        getInLayerOneLayerNoCrossings();
        countCrossingsInLayerForUpperNodeLowerNode(0, 0, 1);
        assertThat("upperLowerCrossings", upperLowerCrossings, is(0));
        assertThat("lowerUpperCrossings", lowerUpperCrossings, is(1));
        countCrossingsInLayerForUpperNodeLowerNode(0, 1, 2);
        assertThat("upperLowerCrossings", upperLowerCrossings, is(0));
        assertThat("lowerUpperCrossings", lowerUpperCrossings, is(0));
        countCrossingsInLayerForUpperNodeLowerNode(0, 2, 3);
        assertThat("upperLowerCrossings", upperLowerCrossings, is(0));
        assertThat("lowerUpperCrossings", lowerUpperCrossings, is(1));
    }

    @Test
    public void fixedPortOrderCrossingToInLayerEdge() {
        getInLayerEdgesFixedPortOrderInLayerCrossing();

        countCrossingsInLayerForUpperNodeLowerNode(0, 1, 2);

        assertThat("upperLowerCrossings", upperLowerCrossings, is(1));
        assertThat("lowerUpperCrossings", lowerUpperCrossings, is(0));
    }

    @Test
    public void fixedPortOrderTwoInLayerEdgesCrossEachOther() {
        getFixedPortOrderTwoInLayerEdgesCrossEachOther();

        countCrossingsInLayerForUpperNodeLowerNode(0, 0, 1);

        assertThat("upperLowerCrossings", upperLowerCrossings, is(1));
        assertThat("lowerUpperCrossings", lowerUpperCrossings, is(0));
    }

    /**
     * <pre>
     * *------
     * *---\ |
     * ___ | |
     * | |--/
     * |_|
     * </pre>
     */
    @Test
    public void multipleEdgesIntoOnePort_ShouldNotCauseCrossing() {
        LNode[] nodes = addNodesToLayer(3, makeLayer(getGraph()));
        final PortSide portSide = PortSide.EAST;
        LPort portOne = addPortOnSide(nodes[0], portSide);
        LPort portTwo = addPortOnSide(nodes[1], portSide);
        LPort portThree = addPortOnSide(nodes[2], portSide);
        addEdgeBetweenPorts(portOne, portThree);
        addEdgeBetweenPorts(portTwo, portThree);

        countCrossingsInLayerForUpperNodeLowerNode(0, 0, 1);

        assertThat("upperLowerCrossings", upperLowerCrossings, is(0));
        assertThat("lowerUpperCrossings", lowerUpperCrossings, is(0));
    }

    @Ignore
    // this case is not supported by the current algorithm
    public void multipleEdgesIntoOnePortAndFreePortOrderCausesCrossings() {
        multipleEdgesIntoOnePortAndFreePortOrder();

        countCrossingsInLayerForUpperNodeLowerNode(0, 0, 1);

        assertThat("upperLowerCrossings", upperLowerCrossings, is(1));
        assertThat("lowerUpperCrossings", lowerUpperCrossings, is(0));
    }

    private void countCrossingsInLayerForUpperNodeLowerNode(final int layerIndex, final int upperNodeIndex,
            final int lowerNodeIndex) {

        initCrossingCounterForLayerIndex(layerIndex);

        countCrossings(upperNodeIndex, lowerNodeIndex);
    }

    private void countCrossings(final int upperNodeIndex, final int lowerNodeIndex) {
        Pair<Integer, Integer> leftCrossings = leftCounter.countInLayerCrossingsBetweenNodesInBothOrders(
                nodeOrder[upperNodeIndex], nodeOrder[lowerNodeIndex],
                PortSide.WEST);
        Pair<Integer, Integer> rightCrossings = rightCounter.countInLayerCrossingsBetweenNodesInBothOrders(
                nodeOrder[upperNodeIndex], nodeOrder[lowerNodeIndex],
                PortSide.EAST);
        upperLowerCrossings = leftCrossings.getFirst() + rightCrossings.getFirst();
        lowerUpperCrossings = leftCrossings.getSecond() + rightCrossings.getSecond();
    }

    /**
     * Initializes Crossing counter, sets nodeOrder with the nodes from the
     * layer being considered and numbers its in ascending form as required by
     * the counter.
     *
     * @param layerIndex
     */
    private void initCrossingCounterForLayerIndex(final int layerIndex) {
        LNode[][] currentOrder = getGraph().toNodeArray();
        nodeOrder = currentOrder[layerIndex];
        numberIdsAscendinglyIn(nodeOrder);
        leftCounter = new CrossingsCounter(new int[getNPorts(currentOrder)]);
        rightCounter =  new CrossingsCounter(new int[getNPorts(currentOrder)]);
        leftCounter.initPortPositionsForInLayerCrossings(nodeOrder, PortSide.WEST);
        rightCounter.initPortPositionsForInLayerCrossings(nodeOrder, PortSide.EAST);
    }

    private int getNPorts(final LNode[][] currentOrder) {
        int nPorts = 0;
        for (LNode[] lNodes : currentOrder) {
            for (LNode n : lNodes) {
                nPorts += n.getPorts().size();
            }
        }
        return nPorts;
    }

    private void numberIdsAscendinglyIn(final LNode[] nodes) {
        for (int i = 0; i < nodes.length; i++) {
            nodes[i].id = i;
        }
    }

    private void switchOrderAndNotifyCounter(final int indexOne, final int indexTwo) {
        leftCounter.switchNodes(nodeOrder[indexOne], nodeOrder[indexTwo], PortSide.WEST);
        rightCounter.switchNodes(nodeOrder[indexOne], nodeOrder[indexTwo], PortSide.EAST);
        LNode one = nodeOrder[indexOne];
        nodeOrder[indexOne] = nodeOrder[indexTwo];
        nodeOrder[indexTwo] = one;
    }
}
