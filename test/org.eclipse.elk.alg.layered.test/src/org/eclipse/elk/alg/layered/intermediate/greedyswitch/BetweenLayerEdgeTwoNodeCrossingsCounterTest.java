/*******************************************************************************
 * Copyright (c) 2016 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.alg.layered.intermediate.greedyswitch;

import static org.hamcrest.CoreMatchers.is;
import static org.junit.Assert.assertThat;

import java.util.List;

import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.core.options.PortSide;
import org.junit.Ignore;
import org.junit.Test;

/**
 * Tests counting crosses of edges in between layers for neighbouring nodes.
 *
 * @author alan
 */
public class BetweenLayerEdgeTwoNodeCrossingsCounterTest extends TestGraphCreator {
    private BetweenLayerEdgeTwoNodeCrossingsCounter crossingCounter;
    private LNode upperNode;
    private LNode lowerNode;
    private Layer layerToCountIn;
    private LNode[][] nodeOrder;

    // CHECKSTYLEOFF javadoc
    // CHECKSTYLEOFF MagicNumber

    @Test
    public void twoNodeNoEdges() {
        getTwoNodesNoConnectionGraph();

        nodeOrder = getGraph().toNodeArray();
        layerToCountIn = getGraph().getLayers().get(0);
        crossingCounter = new BetweenLayerEdgeTwoNodeCrossingsCounter(nodeOrder, 0);
        setUpperNode(0);
        setLowerNode(1);

        assertBothSideUpperLowerCrossingsIs(0);
        assertBothSideLowerUpperCrossingsIs(0);
        assertWesternSideUpperLowerCrossingsIs(0);
        assertWesternSideLowerUpperCrossingsIs(0);
        assertEasternSideUpperLowerCrossingsIs(0);
        assertEasternSideLowerUpperCrossingsIs(0);
    }

    @Test
    public void crossFormed() {
        getCrossFormedGraph();

        nodeOrder = getGraph().toNodeArray();
        layerToCountIn = getGraph().getLayers().get(1);
        crossingCounter = new BetweenLayerEdgeTwoNodeCrossingsCounter(nodeOrder, 1);
        setUpperNode(0);
        setLowerNode(1);

        assertBothSideUpperLowerCrossingsIs(1);
        assertBothSideLowerUpperCrossingsIs(0);
        assertWesternSideUpperLowerCrossingsIs(1);
        assertWesternSideLowerUpperCrossingsIs(0);
        assertEasternSideUpperLowerCrossingsIs(0);
        assertEasternSideLowerUpperCrossingsIs(0);
    }

    @Test
    public void oneNode() {
        getOneNodeGraph();

        nodeOrder = getGraph().toNodeArray();
        layerToCountIn = getGraph().getLayers().get(0);
        crossingCounter = new BetweenLayerEdgeTwoNodeCrossingsCounter(nodeOrder, 0);
        setUpperNode(0);
        setLowerNode(0);

        assertBothSideUpperLowerCrossingsIs(0);
        assertBothSideLowerUpperCrossingsIs(0);
        assertWesternSideUpperLowerCrossingsIs(0);
        assertWesternSideLowerUpperCrossingsIs(0);
        assertEasternSideUpperLowerCrossingsIs(0);
        assertEasternSideLowerUpperCrossingsIs(0);
    }

    @Test
    public void crossFormedMultipleEdgesBetweenSameNodes() {
        getMultipleEdgesBetweenSameNodesGraph();

        nodeOrder = getGraph().toNodeArray();
        layerToCountIn = getGraph().getLayers().get(1);
        crossingCounter = new BetweenLayerEdgeTwoNodeCrossingsCounter(nodeOrder, 1);
        setUpperNode(0);
        setLowerNode(1);

        assertBothSideUpperLowerCrossingsIs(4);
        assertBothSideLowerUpperCrossingsIs(0);
        assertWesternSideUpperLowerCrossingsIs(4);
        assertWesternSideLowerUpperCrossingsIs(0);
        assertEasternSideUpperLowerCrossingsIs(0);
        assertEasternSideLowerUpperCrossingsIs(0);
    }

    @Test
    public void crossWithExtraEdgeInBetween() {
        getCrossWithExtraEdgeInBetweenGraph();

        nodeOrder = getGraph().toNodeArray();
        layerToCountIn = getGraph().getLayers().get(1);
        crossingCounter = new BetweenLayerEdgeTwoNodeCrossingsCounter(nodeOrder, 1);
        setUpperNode(0);
        setLowerNode(2);

        assertBothSideUpperLowerCrossingsIs(1);
        assertBothSideLowerUpperCrossingsIs(0);
        assertWesternSideUpperLowerCrossingsIs(1);
        assertWesternSideLowerUpperCrossingsIs(0);
        assertEasternSideUpperLowerCrossingsIs(0);
        assertEasternSideLowerUpperCrossingsIs(0);
    }

    @Test
    public void ignoreInLayerEdges() {
        getInLayerEdgesGraph();

        nodeOrder = getGraph().toNodeArray();
        layerToCountIn = getGraph().getLayers().get(1);
        crossingCounter = new BetweenLayerEdgeTwoNodeCrossingsCounter(nodeOrder, 1);
        setUpperNode(0);
        setLowerNode(2);

        assertBothSideUpperLowerCrossingsIs(0);
        assertBothSideLowerUpperCrossingsIs(0);
        assertWesternSideUpperLowerCrossingsIs(0);
        assertWesternSideLowerUpperCrossingsIs(0);
        assertEasternSideUpperLowerCrossingsIs(0);
        assertEasternSideLowerUpperCrossingsIs(0);
    }

    @Test
    public void ignoreSelfLoops() {
        getCrossWithManySelfLoopsGraph();

        nodeOrder = getGraph().toNodeArray();
        layerToCountIn = getGraph().getLayers().get(1);
        crossingCounter = new BetweenLayerEdgeTwoNodeCrossingsCounter(nodeOrder, 1);
        setUpperNode(0);
        setLowerNode(1);

        assertBothSideUpperLowerCrossingsIs(1);
        assertBothSideLowerUpperCrossingsIs(0);
        assertWesternSideUpperLowerCrossingsIs(1);
        assertWesternSideLowerUpperCrossingsIs(0);
        assertEasternSideUpperLowerCrossingsIs(0);
        assertEasternSideLowerUpperCrossingsIs(0);
    }

    @Test
    public void moreComplexThreeLayerGraph() {
        getMoreComplexThreeLayerGraph();

        nodeOrder = getGraph().toNodeArray();
        layerToCountIn = getGraph().getLayers().get(1);
        crossingCounter = new BetweenLayerEdgeTwoNodeCrossingsCounter(nodeOrder, 1);
        setUpperNode(0);
        setLowerNode(1);

        assertWesternSideUpperLowerCrossingsIs(1);
        assertWesternSideLowerUpperCrossingsIs(1);
        assertEasternSideUpperLowerCrossingsIs(2);
        assertEasternSideLowerUpperCrossingsIs(3);
        assertBothSideUpperLowerCrossingsIs(3);
        assertBothSideLowerUpperCrossingsIs(4);
    }

    @Test
    public void fixedPortOrder() {
        getFixedPortOrderGraph();

        nodeOrder = getGraph().toNodeArray();
        layerToCountIn = getGraph().getLayers().get(1);
        crossingCounter = new BetweenLayerEdgeTwoNodeCrossingsCounter(nodeOrder, 1);
        setUpperNode(0);
        setLowerNode(1);

        assertEasternSideUpperLowerCrossingsIs(0);
        assertEasternSideLowerUpperCrossingsIs(0);
        assertWesternSideUpperLowerCrossingsIs(1);
        assertWesternSideLowerUpperCrossingsIs(0);
        assertBothSideUpperLowerCrossingsIs(1);
        assertBothSideLowerUpperCrossingsIs(0);
    }

    @Test
    public void switchThreeTimes() {
        LNode[] leftNodes = addNodesToLayer(2, makeLayer(getGraph()));
        LNode[] rightNodes = addNodesToLayer(4, makeLayer(getGraph()));
        LPort leftTopPort = addPortOnSide(leftNodes[0], PortSide.EAST);
        LPort leftLowerPort = addPortOnSide(leftNodes[1], PortSide.EAST);
        LPort rightTopPort = addPortOnSide(rightNodes[0], PortSide.WEST);

        addEdgeBetweenPorts(leftLowerPort, rightTopPort);
        eastWestEdgeFromTo(leftLowerPort, rightNodes[2]);

        addEdgeBetweenPorts(leftTopPort, rightTopPort);
        eastWestEdgeFromTo(leftTopPort, rightNodes[1]);
        eastWestEdgeFromTo(leftTopPort, rightNodes[3]);

        setUpIds();
        nodeOrder = getGraph().toNodeArray();
        layerToCountIn = getGraph().getLayers().get(0);
        crossingCounter = new BetweenLayerEdgeTwoNodeCrossingsCounter(nodeOrder, 0);
        setUpperNode(0);
        setLowerNode(1);

        assertEasternSideUpperLowerCrossingsIs(3);
        assertEasternSideLowerUpperCrossingsIs(2);
        assertWesternSideUpperLowerCrossingsIs(0);
        assertWesternSideLowerUpperCrossingsIs(0);
        assertBothSideUpperLowerCrossingsIs(3);
        assertBothSideLowerUpperCrossingsIs(2);
    }

    @Test
    public void intoSamePort() {
        twoEdgesIntoSamePort();
        nodeOrder = getGraph().toNodeArray();
        layerToCountIn = getGraph().getLayers().get(1);
        crossingCounter = new BetweenLayerEdgeTwoNodeCrossingsCounter(nodeOrder, 1);
        setUpperNode(0);
        setLowerNode(1);

        assertEasternSideUpperLowerCrossingsIs(0);
        assertEasternSideLowerUpperCrossingsIs(0);
        assertWesternSideUpperLowerCrossingsIs(2);
        assertWesternSideLowerUpperCrossingsIs(0);
        assertBothSideUpperLowerCrossingsIs(2);
        assertBothSideLowerUpperCrossingsIs(0);
    }

    @Test
    public void intoSamePortCausesCrossingsOnSwitch() {
        twoEdgesIntoSamePortCrossesWhenSwitched();
        nodeOrder = getGraph().toNodeArray();
        layerToCountIn = getGraph().getLayers().get(0);
        crossingCounter = new BetweenLayerEdgeTwoNodeCrossingsCounter(nodeOrder, 0);
        setUpperNode(0);
        setLowerNode(1);

        assertEasternSideUpperLowerCrossingsIs(0);
        assertEasternSideLowerUpperCrossingsIs(1);
        assertWesternSideUpperLowerCrossingsIs(0);
        assertWesternSideLowerUpperCrossingsIs(0);
        assertBothSideUpperLowerCrossingsIs(0);
        assertBothSideLowerUpperCrossingsIs(1);
    }

    @Test
    public void intoSamePortReducesCrossingsOnSwitch() {
        twoEdgesIntoSamePortResolvesCrossingWhenSwitched();
        nodeOrder = getGraph().toNodeArray();
        layerToCountIn = getGraph().getLayers().get(0);
        crossingCounter = new BetweenLayerEdgeTwoNodeCrossingsCounter(nodeOrder, 0);
        setUpperNode(0);
        setLowerNode(1);

        assertEasternSideUpperLowerCrossingsIs(1);
        assertEasternSideLowerUpperCrossingsIs(0);
        assertWesternSideUpperLowerCrossingsIs(0);
        assertWesternSideLowerUpperCrossingsIs(0);
        assertBothSideUpperLowerCrossingsIs(1);
        assertBothSideLowerUpperCrossingsIs(0);
    }

    @Test
    public void intoSamePortFromEastSwitchWithFixedPortOrder() {
        twoEdgesIntoSamePortFromEastWithFixedPortOrder();
        nodeOrder = getGraph().toNodeArray();
        layerToCountIn = getGraph().getLayers().get(0);
        crossingCounter = new BetweenLayerEdgeTwoNodeCrossingsCounter(nodeOrder, 0);
        setUpperNode(0);
        setLowerNode(1);

        assertEasternSideUpperLowerCrossingsIs(0);
        assertEasternSideLowerUpperCrossingsIs(1);
        assertWesternSideUpperLowerCrossingsIs(0);
        assertWesternSideLowerUpperCrossingsIs(0);
        assertBothSideUpperLowerCrossingsIs(0);
        assertBothSideLowerUpperCrossingsIs(1);
    }

    @Test
    public void multipleEdgesIntoSamePort_causesNoCrossings() {
        Layer leftLayer = makeLayer(graph);
        Layer rightLayer = makeLayer(graph);
        
        LNode topLeft = addNodeToLayer(leftLayer);
        LNode bottomLeft = addNodeToLayer(leftLayer);
        LNode bottomRight = addNodeToLayer(rightLayer);
        
        LPort bottomRightPort = addPortOnSide(bottomRight, PortSide.WEST);

        eastWestEdgeFromTo(topLeft, bottomRightPort);
        eastWestEdgeFromTo(topLeft, bottomRightPort);
        eastWestEdgeFromTo(bottomLeft, bottomRightPort);
        setUpIds();
        nodeOrder = getGraph().toNodeArray();
        layerToCountIn = getGraph().getLayers().get(0);
        crossingCounter = new BetweenLayerEdgeTwoNodeCrossingsCounter(nodeOrder, 0);
        setUpperNode(0);
        setLowerNode(1);

        assertEasternSideUpperLowerCrossingsIs(0);
        assertEasternSideUpperLowerCrossingsIs(0);
    }

    /**
     * This behavior is not supported by the
     * BetweenLayerEdgeTwoNodeCrossingsCounter.
     *
     * <pre>
     * ___   ____
     * | o---o  |
     * | |\  |  |
     * |_o-+-o  |
     *     | |__|
     *      \
     *       *
     * </pre>
     *
     * Right layer has fixed port order, left has free.
     */
    @Ignore
    public void partlyFixedPartlyUnfixedPortOrder() {
        LNode leftNode = addNodeToLayer(makeLayer(getGraph()));
        LNode[] rightNodes = addNodesToLayer(2, makeLayer(getGraph()));
        setFixedOrderConstraint(rightNodes[0]);
        setFixedOrderConstraint(rightNodes[1]);
        LPort leftTopPort = addPortOnSide(leftNode, PortSide.EAST);
        LPort leftLowerPort = addPortOnSide(leftNode, PortSide.EAST);
        LPort rightUpperNodeLowerPort = addPortOnSide(rightNodes[0], PortSide.WEST);
        LPort rightUpperNodeTopPort = addPortOnSide(rightNodes[0], PortSide.WEST);
        LPort rightLowerPort = addPortOnSide(rightNodes[1], PortSide.WEST);
        addEdgeBetweenPorts(leftTopPort, rightUpperNodeTopPort);
        addEdgeBetweenPorts(leftTopPort, rightLowerPort);
        addEdgeBetweenPorts(leftLowerPort, rightUpperNodeLowerPort);

        nodeOrder = getGraph().toNodeArray();
        layerToCountIn = getGraph().getLayers().get(1);
        crossingCounter = new BetweenLayerEdgeTwoNodeCrossingsCounter(nodeOrder, 1);
        setUpperNode(0);
        setLowerNode(1);

        assertBothSideUpperLowerCrossingsIs(1);
    }

    private void assertEasternSideLowerUpperCrossingsIs(final int expectedCrossings) {
        crossingCounter.countEasternEdgeCrossings(upperNode, lowerNode);
        int crossings = crossingCounter.getLowerUpperCrossings();
        assertThat(failMessage(" east, lower upper "), crossings, is(expectedCrossings));
    }

    private void assertEasternSideUpperLowerCrossingsIs(final int expectedCrossings) {
        crossingCounter.countEasternEdgeCrossings(upperNode, lowerNode);
        int crossings = crossingCounter.getUpperLowerCrossings();
        assertThat(failMessage(" east, upper lower "), crossings, is(expectedCrossings));
    }

    private void assertWesternSideLowerUpperCrossingsIs(final int expectedCrossings) {
        crossingCounter.countWesternEdgeCrossings(upperNode, lowerNode);
        int crossings = crossingCounter.getLowerUpperCrossings();
        assertThat(failMessage(" west, lower upper "), crossings, is(expectedCrossings));
    }

    private void assertWesternSideUpperLowerCrossingsIs(final int expectedCrossings) {
        crossingCounter.countWesternEdgeCrossings(upperNode, lowerNode);
        int crossings = crossingCounter.getUpperLowerCrossings();
        assertThat(failMessage(" west, upper lower "), crossings, is(expectedCrossings));
    }

    private void assertBothSideLowerUpperCrossingsIs(final int expectedCrossings) {
        crossingCounter.countBothSideCrossings(upperNode, lowerNode);
        int crossings = crossingCounter.getLowerUpperCrossings();
        assertThat(failMessage(" both, lower upper "), crossings, is(expectedCrossings));
    }

    private void assertBothSideUpperLowerCrossingsIs(final int expectedCrossings) {
        crossingCounter.countBothSideCrossings(upperNode, lowerNode);
        int crossings = crossingCounter.getUpperLowerCrossings();
        assertThat(failMessage(" both, upper lower "), crossings, is(expectedCrossings));
    }

    private void setUpperNode(final int nodeIndex) {
        List<LNode> nodes = layerToCountIn.getNodes();
        upperNode = nodes.get(nodeIndex);
    }

    private void setLowerNode(final int nodeIndex) {
        List<LNode> nodes = layerToCountIn.getNodes();
        lowerNode = nodes.get(nodeIndex);
    }

    private String failMessage(final String configuration) {
        return "With configuration: " + configuration + " in Layer " + layerToCountIn
                + " between the nodes: Upper node " + upperNode.getIndex() + " and lower Node " + lowerNode.getIndex();
    }
}
