/*******************************************************************************
 * Copyright (c) 2016 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.alg.layered.p3order.counting;

import static org.hamcrest.CoreMatchers.is;
import static org.junit.Assert.assertThat;

import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.intermediate.greedyswitch.NorthSouthEdgeTestGraphCreator;
import org.eclipse.elk.core.data.LayoutMetaDataService;
import org.eclipse.elk.core.options.PortSide;
import org.junit.Before;
import org.junit.Test;

/**
 * Tests the north/south edge crossing counting of the {@link CrossingsCounter}.
 */
@SuppressWarnings("restriction")
public class NorthSouthEdgeAllCrossingsCounterTest extends NorthSouthEdgeTestGraphCreator {

    // CHECKSTYLEOFF Javadoc
    // CHECKSTYLEOFF MagicNumber
    // CHECKSTYLEOFF MethodName
    
    @Before
    public void setup() {
        LayoutMetaDataService.getInstance();
    }

    @Test
    public void northernNorthSouthNodeSingleCrossing() {
        getNorthSouthUpwardCrossingGraph();
        int crossingCount = initCounterForLayerWithIndexAndCountInLayer(0);
        assertThat(crossingCount, is(1));
    }

    @Test
    public void northernNorthSouthNodeMultipleCrossings() {
        getNorthSouthUpwardMultipleCrossingGraph();
        int crossingCount = initCounterForLayerWithIndexAndCountInLayer(0);
        assertThat(crossingCount, is(3));
    }

    @Test
    public void southernTowEdgeEastCrossing() {
        getNorthSouthDownwardCrossingGraph();
        int crossingCount = initCounterForLayerWithIndexAndCountInLayer(0);
        assertThat(crossingCount, is(1));
    }

    @Test
    public void southernNorthSouthMultipleNodeCrossing() {
        getNorthSouthDownwardMultipleCrossingGraph();
        int crossingCount = initCounterForLayerWithIndexAndCountInLayer(0);
        assertThat(crossingCount, is(3));
    }

    @Test
    public void southernTwoWesternEdges() {
        getNorthSouthSouthernTwoWesternEdges();
        int crossingCount = initCounterForLayerWithIndexAndCountInLayer(1);
        assertThat(crossingCount, is(1));
    }

    @Test
    public void southernThreeWesternEdges() {
        getNorthSouthSouthernThreeWesternEdges();
        int crossingCount = initCounterForLayerWithIndexAndCountInLayer(1);
        assertThat(crossingCount, is(3));
    }

    @Test
    public void northSouthEdgesComeFromBothSidesDontCross() {
        getSouthernNorthSouthGraphEdgesFromEastAndWestNoCrossings();
        int crossingCount = initCounterForLayerWithIndexAndCountInLayer(1);
        assertThat(crossingCount, is(0));
    }

    @Test
    public void southernNorthSouthEdgesBothToEastDontCross() {
        getSouthernNorthSouthEdgesBothToEast();
        int crossingCount = initCounterForLayerWithIndexAndCountInLayer(0);
        assertThat(crossingCount, is(0));
    }

    @Test
    public void northSouthEdgesComeFromBothSidesDoCross() {
        getNorthSouthEdgesFromEastAndWestAndCross();

        int crossingCount = initCounterForLayerWithIndexAndCountInLayer(1);

        assertThat(crossingCount, is(1));
    }

    @Test
    public void northernBothEdgesWestern() {
        getNorthSouthNorthernWesternEdges();
        int crossingCount = initCounterForLayerWithIndexAndCountInLayer(0);
        assertThat(crossingCount, is(0));
    }

    @Test
    public void northernEasternPortToWestWesternPortToEast() {
        getNorthSouthNorthernEasternPortToWestWesternPortToEast();
        int crossingCount = initCounterForLayerWithIndexAndCountInLayer(1);

        assertThat(crossingCount, is(1));
    }

    @Test
    public void allSidesMultipleCrossings() {
        getNorthSouthAllSidesMultipleCrossings();
        int crossingCount = initCounterForLayerWithIndexAndCountInLayer(1);
        assertThat(crossingCount, is(4));
    }

    @Test
    public void oneEdgeDummyIsCrossedByOneSouthernNorthSouthPortEdge() {
        getSouthernNorthSouthDummyEdgeCrossingGraph();
        int crossingCount = initCounterForLayerWithIndexAndCountInLayer(1);
        assertThat(crossingCount, is(1));
    }

    @Test
    public void oneEdgeDummyIsCrossedByTwoSouthernNorthSouthPortEdges() {
        getSouthernNorthSouthDummyEdgeTwoCrossingGraph();
        int crossingCount = initCounterForLayerWithIndexAndCountInLayer(1);
        assertThat(crossingCount, is(2));
    }

    @Test
    public void southernTwoDummyEdgeAndTwoNorthSouthShouldCrossFourTimes() {
        getSouthernTwoDummyEdgeAndNorthSouthCrossingGraph();
        int crossingCount = initCounterForLayerWithIndexAndCountInLayer(1);
        assertThat(crossingCount, is(4));
    }

    @Test
    public void normalNodesNorthSouthEdgesHaveCrossingsToLongEdgeDummyOnBothSides() {
        getMultipleNorthSouthAndLongEdgeDummiesOnBothSides();
        int crossingCount = initCounterForLayerWithIndexAndCountInLayer(1);
        assertThat(crossingCount, is(4));
    }

    @Test
    public void ignoresUnconnectedPortsForNormalNodeAndLongEdgeDummies() {
        getLongEdgeDummyAndNormalNodeWithUnusedPortsOnSouthernSide();
        int crossingCount = initCounterForLayerWithIndexAndCountInLayer(1);
        assertThat(crossingCount, is(0));
    }

    @Test
    public void noNorthSouthNode() {
        getCrossFormedGraph();
        int crossingCount = initCounterForLayerWithIndexAndCountInLayer(0);
        assertThat(crossingCount, is(0));
    }

    /**
     * <pre>
     * *---o
     *     +--*  
     *    /| 
     *   o-+--*
     *   | |
     *  _|_|_
     *  |   |
     *  |___|
     *  .
     * </pre>
     *
     */
    @Test
    public void moreThanOneEdgeIntoNSNode_countsTheseToo() {
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

        int crossingCount = initCounterForLayerWithIndexAndCountInLayer(1);

        assertThat(crossingCount, is(2));
    }

    /**
     * *--------o
     * *------o |
     * *----o | |
     * *--o | | |
     *    | | | |
     *   _|_|_|_|_ 
     *  |         |
     */
    @Test 
    public void theOneThatFailedWithTheOldCounting() {
        LNode[] leftNodes = addNodesToLayer(4, makeLayer(getGraph()));
        LNode[] middleNodes = addNodesToLayer(5, makeLayer(getGraph()));
        
        setFixedOrderConstraint(middleNodes[4]);
        
        for (int i = 3; i >= 0; --i) {
            addNorthSouthEdge(PortSide.NORTH, middleNodes[4], middleNodes[i], leftNodes[i], false);
        }
        
        int crossingCount = initCounterForLayerWithIndexAndCountInLayer(1);
        
        assertThat(crossingCount, is(0));
    }
    
    private int initCounterForLayerWithIndexAndCountInLayer(final int layerIndex) {
        setUpIds();
        int numPorts = 0;
        for (Layer l : getGraph()) {
            for (LNode lNode : l) {
                numPorts += lNode.getPorts().size();
            }
        }
        CrossingsCounter counter = new CrossingsCounter(new int[numPorts]);
        return counter.countNorthSouthPortCrossingsInLayer(getGraph().toNodeArray()[layerIndex]);
    }

}
