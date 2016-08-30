/*******************************************************************************
 * Copyright (c) 2011, 2015 Kiel University and others.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the Eclipse Public License v1.0
 * which accompanies this distribution, and is available at
 * http://www.eclipse.org/legal/epl-v10.html
 *
 * Contributors:
 *     Kiel University - initial API and implementation
 *******************************************************************************/
package org.eclipse.elk.alg.test.layered.intermediate.greedyswitch;

import static org.hamcrest.CoreMatchers.is;
import static org.junit.Assert.assertThat;

import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.intermediate.greedyswitch.NorthSouthEdgeNeighbouringNodeCrossingsCounter;
import org.eclipse.elk.alg.layered.properties.InternalProperties;
import org.eclipse.elk.alg.layered.properties.LayeredOptions;
import org.eclipse.elk.core.options.EdgeRouting;
import org.eclipse.elk.core.options.PortSide;
import org.junit.Ignore;
import org.junit.Test;

/**
 * Tests counting crosses generated by the ordering of north south ports.
 *
 * @author alan
 *
 */
public class NorthSouthEdgeNeighbouringNodeCrossingsCounterTest extends NorthSouthEdgeTestGraphCreator {
    private NorthSouthEdgeNeighbouringNodeCrossingsCounter counter;
    private LNode[] layer;

    // CHECKSTYLEOFF javadoc
    // CHECKSTYLEOFF MagicNumber
    // CHECKSTYLEOFF MethodName

    @Test
    public void noNorthSouthNode() {
        getCrossFormedGraph();
        countCrossingsInLayerBetweenNodes(0, 0, 1);
        assertThat(counter.getUpperLowerCrossings(), is(0));
        assertThat(counter.getLowerUpperCrossings(), is(0));
    }

    @Test
    public void southernNorthSouthNodeCrossing() {
        getNorthSouthDownwardCrossingGraph();
        countCrossingsInLayerBetweenNodes(0, 1, 2);
        assertThat(counter.getUpperLowerCrossings(), is(1));
        assertThat(counter.getLowerUpperCrossings(), is(0));
    }

    @Test
    public void northernNorthSouthNodeCrossings() {
        getNorthSouthUpwardCrossingGraph();

        countCrossingsInLayerBetweenNodes(0, 0, 1);

        assertThat(counter.getUpperLowerCrossings(), is(1));
        assertThat(counter.getLowerUpperCrossings(), is(0));
    }

    @Test
    public void oneNodeIsLongEdgeDummy() {
        getSouthernNorthSouthDummyEdgeCrossingGraph();
        countCrossingsInLayerBetweenNodes(1, 1, 2);
        assertThat(counter.getUpperLowerCrossings(), is(1));
        assertThat(counter.getLowerUpperCrossings(), is(0));

        switchNodes(1, 2);

        countCrossingsInLayerBetweenNodes(1, 1, 2);
        assertThat(counter.getUpperLowerCrossings(), is(0));
        assertThat(counter.getLowerUpperCrossings(), is(1));
    }

    @Test
    public void oneNodeIsLongEdgeDummyNorthern() {
        getNorthernNorthSouthDummyEdgeCrossingGraph();
        countCrossingsInLayerBetweenNodes(1, 0, 1);
        assertThat(counter.getUpperLowerCrossings(), is(1));
        assertThat(counter.getLowerUpperCrossings(), is(0));

        switchNodes(0, 1);

        countCrossingsInLayerBetweenNodes(1, 0, 1);
        assertThat(counter.getUpperLowerCrossings(), is(0));
        assertThat(counter.getLowerUpperCrossings(), is(1));
    }

    @Test
    public void withNormalNode() {
        getNorthSouthDownwardCrossingGraph();
        countCrossingsInLayerBetweenNodes(0, 0, 1);
        assertThat(counter.getUpperLowerCrossings(), is(0));
        assertThat(counter.getLowerUpperCrossings(), is(0));
    }

    @Test
    public void northSouthEdgesComeFromBothSidesDontCross() {
        getSouthernNorthSouthGraphEdgesFromEastAndWestNoCrossings();
        countCrossingsInLayerBetweenNodes(1, 1, 2);
        assertThat(counter.getUpperLowerCrossings(), is(0));
        assertThat(counter.getLowerUpperCrossings(), is(0));

        getNorthernNorthSouthGraphEdgesFromEastAndWestNoCrossings();
        countCrossingsInLayerBetweenNodes(1, 0, 1);
        assertThat(counter.getUpperLowerCrossings(), is(0));
        assertThat(counter.getLowerUpperCrossings(), is(0));
    }

    @Test
    public void southernNorthSouthEdgesBothToEast() {
        getSouthernNorthSouthEdgesBothToEast();
        countCrossingsInLayerBetweenNodes(0, 1, 2);
        assertThat(counter.getUpperLowerCrossings(), is(0));
        assertThat(counter.getLowerUpperCrossings(), is(1));
    }

    @Test
    public void crossingsWithNorthSouthPortsBelongingToDifferentNodesShouldNotBeCounted() {
        getGraphWhereLayoutUnitPreventsSwitch();
        countCrossingsInLayerBetweenNodes(0, 1, 2);
        assertThat(counter.getUpperLowerCrossings(), is(0));
        assertThat(counter.getLowerUpperCrossings(), is(0));
    }

    @Test
    public void northSouthEdgesComeFromBothSidesDoCross() {
        getNorthSouthEdgesFromEastAndWestAndCross();
        countCrossingsInLayerBetweenNodes(1, 1, 2);
        assertThat(counter.getUpperLowerCrossings(), is(1));
        assertThat(counter.getLowerUpperCrossings(), is(1));

    }

    @Test
    public void switchNodesAndRecount() {
        getNorthSouthUpwardCrossingGraph();
        countCrossingsInLayerBetweenNodes(0, 0, 1);
        assertThat(counter.getUpperLowerCrossings(), is(1));
        assertThat(counter.getLowerUpperCrossings(), is(0));
        switchAndRecount(0, 1);
        assertThat(counter.getUpperLowerCrossings(), is(0));
        assertThat(counter.getLowerUpperCrossings(), is(1));
    }

    @Test
    public void southPortOndNormalNodeBelowLongEdgeDummy() {
        getSouthPortOnNormalNodeBelowLongEdgeDummy();

        countCrossingsInLayerBetweenNodes(1, 0, 1);
        assertThat(counter.getUpperLowerCrossings(), is(0));
        assertThat(counter.getLowerUpperCrossings(), is(1));
        switchAndRecount(0, 1);
        assertThat(counter.getUpperLowerCrossings(), is(1));
        assertThat(counter.getLowerUpperCrossings(), is(0));
    }

    @Test
    public void northPortOndNormalNodeAboveLongEdgeDummy() {
        getNorthPortOndNormalNodeAboveLongEdgeDummy();

        countCrossingsInLayerBetweenNodes(1, 1, 2);
        assertThat(counter.getUpperLowerCrossings(), is(0));
        assertThat(counter.getLowerUpperCrossings(), is(1));
        switchAndRecount(1, 2);
        assertThat(counter.getUpperLowerCrossings(), is(1));
        assertThat(counter.getLowerUpperCrossings(), is(0));
    }

    private void switchAndRecount(final int upperNodeIndex, final int lowerNodeIndex) {
        switchNodes(upperNodeIndex, lowerNodeIndex);
        counter.countCrossings(layer[upperNodeIndex], layer[lowerNodeIndex]);
    }

    @Test
    public void southernTwoWesternEdges() {
        getNorthSouthSouthernTwoWesternEdges();
        countCrossingsInLayerBetweenNodes(1, 1, 2);
        assertThat(counter.getUpperLowerCrossings(), is(1));
        assertThat(counter.getLowerUpperCrossings(), is(0));
        switchAndRecount(1, 2);
        assertThat(counter.getUpperLowerCrossings(), is(0));
        assertThat(counter.getLowerUpperCrossings(), is(1));
    }

    @Test
    public void southernWesternPortToEastAndEasternPortToWest() {
        getNorthSouthSouthernWesternPortToEastAndEasternPortToWest();
        countCrossingsInLayerBetweenNodes(1, 1, 2);
        assertThat(counter.getUpperLowerCrossings(), is(1));
        assertThat(counter.getLowerUpperCrossings(), is(1));
        switchAndRecount(1, 2);
        assertThat(counter.getUpperLowerCrossings(), is(1));
        assertThat(counter.getLowerUpperCrossings(), is(1));
    }

    @Test
    public void northernBothEdgesWestern() {
        getNorthSouthNorthernWesternEdges();
        countCrossingsInLayerBetweenNodes(1, 0, 1);
        assertThat(counter.getUpperLowerCrossings(), is(0));
        assertThat(counter.getLowerUpperCrossings(), is(1));
        switchAndRecount(0, 1);
        assertThat(counter.getUpperLowerCrossings(), is(1));
        assertThat(counter.getLowerUpperCrossings(), is(0));
    }

    @Test
    public void northernEasternPortToWestWesternPortToEast() {
        getNorthSouthNorthernEasternPortToWestWesternPortToEast();
        countCrossingsInLayerBetweenNodes(1, 0, 1);
        assertThat(counter.getUpperLowerCrossings(), is(1));
        assertThat(counter.getLowerUpperCrossings(), is(1));
        switchAndRecount(0, 1);
        assertThat(counter.getUpperLowerCrossings(), is(1));
        assertThat(counter.getLowerUpperCrossings(), is(1));
    }

    @Test
    public void normalNodesNorthSouthEdgesHaveCrossingsToLongEdgeDummy() {
        getNorthernNorthSouthDummyEdgeCrossingGraph();

        countCrossingsInLayerBetweenNodes(1, 0, 1);

        assertThat(counter.getUpperLowerCrossings(), is(1));
        assertThat(counter.getLowerUpperCrossings(), is(0));

        countCrossingsInLayerBetweenNodes(1, 1, 2);

        assertThat(counter.getUpperLowerCrossings(), is(1));
        assertThat(counter.getLowerUpperCrossings(), is(0));

        getSouthernNorthSouthDummyEdgeCrossingGraph();

        countCrossingsInLayerBetweenNodes(1, 0, 1);

        assertThat(counter.getUpperLowerCrossings(), is(1));
        assertThat(counter.getLowerUpperCrossings(), is(0));

        countCrossingsInLayerBetweenNodes(1, 1, 2);

        assertThat(counter.getUpperLowerCrossings(), is(1));
        assertThat(counter.getLowerUpperCrossings(), is(0));
    }

    @Test
    public void normalNodesNorthSouthEdgesHaveCrossingsToLongEdgeDummyOnBothSides() {
        getMultipleNorthSouthAndLongEdgeDummiesOnBothSides();
        countCrossingsInLayerBetweenNodes(1, 2, 3);

        assertThat(counter.getUpperLowerCrossings(), is(2));
        assertThat(counter.getLowerUpperCrossings(), is(2));
    }

    @Test
    public void ignoresUnconnectedPortsForNormalNodeAndLongEdgeDummies() {
        getLongEdgeDummyAndNormalNodeWithUnusedPortsOnSouthernSide();
        countCrossingsInLayerBetweenNodes(1, 0, 1);

        assertThat(counter.getUpperLowerCrossings(), is(0));

        getLongEdgeDummyAndNormalNodeWithUnusedPortsOnNorthernSide();
        countCrossingsInLayerBetweenNodes(1, 0, 1);

        assertThat(counter.getUpperLowerCrossings(), is(0));
    }

    @Test
    public void oneEdgeWestOneEdgeEastDontCross() {
        getNorthernNorthSouthGraphEdgesFromEastAndWestNoCrossings();
        countCrossingsInLayerBetweenNodes(1, 0, 1);

        assertThat(counter.getUpperLowerCrossings(), is(0));
        assertThat(counter.getLowerUpperCrossings(), is(0));
    }

    @Test
    public void oneEdgeEastOneEdgeWestDontCross() {
        getNorthernNorthSouthGraphEdgesFromEastAndWestNoCrossingsUpperEdgeEast();

        countCrossingsInLayerBetweenNodes(1, 0, 1);

        assertThat(counter.getUpperLowerCrossings(), is(0));
        assertThat(counter.getLowerUpperCrossings(), is(0));
    }

    /**
     * <pre>
     *
     * *----
     *    /+--*
     *   --+--*
     *   | |
     *  _|_|_
     *  |   |
     *  |___|
     *  .
     * </pre>
     *
     */
    @Ignore // TODO-alan ask is this possible
    public void givenPolylineRoutingWhenMoreThanOneEdgeIntoNSNode_countsTheseToo() {
        LNode leftNode = addNodeToLayer(makeLayer(getGraph()));
        LNode[] middleNodes = addNodesToLayer(3, makeLayer(getGraph()));
        LNode[] rightNodes = addNodesToLayer(2, makeLayer(getGraph()));

        setFixedOrderConstraint(middleNodes[2]);

        // ports are added in clockwise fashion!
        addNorthSouthEdge(PortSide.NORTH, middleNodes[2], middleNodes[1], rightNodes[0], false);
        addNorthSouthEdge(PortSide.NORTH, middleNodes[2], middleNodes[0], leftNode, true);
        // second edge on middle node
        LPort middleNodePort = middleNodes[1].getPorts().get(0);
        eastWestEdgeFromTo(middleNodePort, rightNodes[1]);
        getGraph().setProperty(LayeredOptions.EDGE_ROUTING, EdgeRouting.POLYLINE);

        countCrossingsInLayerBetweenNodes(1, 0, 1);

        assertThat(counter.getUpperLowerCrossings(), is(2));
        assertThat(counter.getLowerUpperCrossings(), is(1));
    }

    /**
     * <pre>
     *     .--*
     *     |
     * *-.-+--*
     *   | |
     *  _|_|_
     *  |   |
     *  |___|
     * </pre>
     */
    @Test
    public void givenMultipleEdgesInOneNSNodeCountsCrossings() {
        LNode leftNode = addNodeToLayer(makeLayer());
        LNode[] middleLayer = addNodesToLayer(3, makeLayer());
        LNode[] rightLayer = addNodesToLayer(2, makeLayer());

        setFixedOrderConstraint(middleLayer[2]);

        addNorthSouthEdge(PortSide.NORTH, middleLayer[2], middleLayer[1], leftNode, true);

        LPort normalNodePort = addPortOnSide(rightLayer[1], PortSide.WEST);
        LPort dummyNodePort = addPortOnSide(middleLayer[1], PortSide.EAST);
        addEdgeBetweenPorts(dummyNodePort, normalNodePort);
        LPort originPort = middleLayer[2].getPorts().get(0);
        dummyNodePort.setProperty(InternalProperties.ORIGIN, originPort);

        addNorthSouthEdge(PortSide.NORTH, middleLayer[2], middleLayer[0], rightLayer[0], false);

        countCrossingsInLayerBetweenNodes(1, 0, 1);

        assertThat(counter.getUpperLowerCrossings(), is(1));
        assertThat(counter.getLowerUpperCrossings(), is(1));
    }

    /**
     * <pre>
     * *---.--*
     *     |
     * *-.-+--*
     *   | |
     *  _|_|_
     *  |   |
     *  |___|
     * </pre>
     */
    @Test
    public void multipleEdgesInBithNSNode() {
        LNode[] leftLayer = addNodesToLayer(2, makeLayer());
        LNode[] middleLayer = addNodesToLayer(3, makeLayer());
        LNode[] rightLayer = addNodesToLayer(2, makeLayer());

        setFixedOrderConstraint(middleLayer[2]);

        addNorthSouthEdge(PortSide.NORTH, middleLayer[2], middleLayer[1], leftLayer[1], true);

        LPort normalNodePort = addPortOnSide(rightLayer[1], PortSide.WEST);
        LPort dummyNodePort = addPortOnSide(middleLayer[1], PortSide.EAST);
        addEdgeBetweenPorts(dummyNodePort, normalNodePort);
        LPort originPort = middleLayer[2].getPorts().get(0);
        dummyNodePort.setProperty(InternalProperties.ORIGIN, originPort);

        addNorthSouthEdge(PortSide.NORTH, middleLayer[2], middleLayer[0], leftLayer[0], true);

        normalNodePort = addPortOnSide(rightLayer[0], PortSide.WEST);
        dummyNodePort = addPortOnSide(middleLayer[0], PortSide.EAST);
        addEdgeBetweenPorts(dummyNodePort, normalNodePort);
        originPort = middleLayer[2].getPorts().get(1);
        dummyNodePort.setProperty(InternalProperties.ORIGIN, originPort);

        countCrossingsInLayerBetweenNodes(1, 0, 1);

        assertThat(counter.getUpperLowerCrossings(), is(1));
        // Note that the two other crossings which result from this switch are
        // counted by a
        // different crossings counter.
        assertThat(counter.getLowerUpperCrossings(), is(1));
    }

    private void switchNodes(final int upper, final int lower) {
        LNode upperNode = layer[upper];
        LNode lowerNode = layer[lower];
        layer[upper] = lowerNode;
        layer[lower] = upperNode;
    }

    private void countCrossingsInLayerBetweenNodes(final int layerIndex, final int upperNodeIndex,
            final int lowerNodeIndex) {
        setLayerIfNotSet(layerIndex);
        counter = new NorthSouthEdgeNeighbouringNodeCrossingsCounter(layer);
        counter.countCrossings(layer[upperNodeIndex], layer[lowerNodeIndex]);

    }

    private void setLayerIfNotSet(final int layerIndex) {
        if (layer == null) {
            layer = getGraph().toNodeArray()[layerIndex];
        }
    }
}
