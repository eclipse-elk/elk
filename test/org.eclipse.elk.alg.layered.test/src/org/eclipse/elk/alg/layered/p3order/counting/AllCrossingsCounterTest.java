/*******************************************************************************
 * Copyright (c) 2016 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.layered.p3order.counting;

import static org.hamcrest.CoreMatchers.is;
import static org.junit.Assert.assertThat;

import java.util.List;

import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.intermediate.greedyswitch.InLayerEdgeTestGraphCreator;
import org.eclipse.elk.alg.layered.intermediate.greedyswitch.NorthSouthEdgeTestGraphCreator;
import org.eclipse.elk.alg.layered.intermediate.greedyswitch.TestGraphCreator;
import org.eclipse.elk.alg.layered.p3order.GraphInfoHolder;
import org.eclipse.elk.alg.layered.p3order.LayerSweepCrossingMinimizer.CrossMinType;
import org.eclipse.elk.alg.test.PlainJavaInitialization;
import org.eclipse.elk.core.options.PortSide;
import org.junit.Before;
import org.junit.Test;

/**
 * Test for extracted and modified AllCrossingsCounter.
 *
 */
public class AllCrossingsCounterTest extends TestGraphCreator {

    // CHECKSTYLEOFF Javadoc
    // CHECKSTYLEOFF MagicNumber
    
    @Before
    public void initPlainJavaLayout() {
        // Necessary as this is a plain junit test for which the LayoutMetaDataService is not initialized automatically
        PlainJavaInitialization.initializePlainJavaLayout();
    }

    @Test
    public void countOneCrossing() {
        graph = getCrossFormedGraph();
        int amountOfCrossings = allCrossings();
        assertThat(amountOfCrossings, is(1));
    }

    @Test
    public void countInLayerCrossing() {
        graph = getInLayerEdgesGraph();
        int amountOfCrossings = allCrossings();
        assertThat(amountOfCrossings, is(1));
    }

    @Test
    public void countInLayerCrossingAndSwitch() {
        graph = getInLayerEdgesGraph();
        int amountOfCrossings = allCrossings();
        assertThat(amountOfCrossings, is(1));
    }

    @Test
    public void countNorthSouthCrossing() {
        graph = new NorthSouthEdgeTestGraphCreator().getNorthSouthDownwardCrossingGraph();
        int amountOfCrossings = allCrossings();
        assertThat(amountOfCrossings, is(1));
    }

    @Test
    public void countNorthernNorthSouthCrossing() {
        graph = new NorthSouthEdgeTestGraphCreator().getNorthSouthUpwardCrossingGraph();
        int amountOfCrossings = allCrossings();
        assertThat(amountOfCrossings, is(1));
    }

    @Test
    public void northSouthDummyEdgeCrossing() {
        graph = new NorthSouthEdgeTestGraphCreator().getSouthernNorthSouthDummyEdgeCrossingGraph();
        int amountOfCrossings = allCrossings();
        assertThat(amountOfCrossings, is(1));
    }

    @Test
    public void switchAndCountTwice() {
        graph = getCrossFormedGraph();
        int amountOfCrossings = allCrossings();
        assertThat(amountOfCrossings, is(1));
        switchNodesInLayer(0, 1, 1);
        amountOfCrossings = allCrossings();
        assertThat(amountOfCrossings, is(0));
    }

    private LNode[][] switchNodesInLayer(final int upperNodeIndex, final int lowerNodeIndex, final int layerIndex) {
        Layer layer = graph.getLayers().get(layerIndex);
        List<LNode> nodes = layer.getNodes();
        LNode upperNode = nodes.get(upperNodeIndex);
        // counter.notifyOfSwitch(upperNode, nodes.get(lowerNodeIndex),
        // layerIndex);

        nodes.set(upperNodeIndex, nodes.get(lowerNodeIndex));
        nodes.set(lowerNodeIndex, upperNode);
        return getCurrentOrder(graph);
    }

    @Test
    public void tooManyInLayerCrossingsWithTheOldMethod() {
        graph = new InLayerEdgeTestGraphCreator().getInLayerOneLayerNoCrossings();
        int amountOfCrossings = allCrossings();
        assertThat(amountOfCrossings, is(0));
    }

    /**
     * <pre>
     * *    *
     *  \\//
     *  //\\
     * *    *
     * </pre>
     */
    @Test
    public void countCrossingsWithMultipleEdgesBetweenSameNodes() {
        LNode[] left = addNodesToLayer(2, makeLayer());
        LNode[] right = addNodesToLayer(2, makeLayer());
        
        LPort[] rightLowerPorts = addPortsOnSide(2, right[1], PortSide.WEST);
        eastWestEdgeFromTo(left[0], rightLowerPorts[1]);
        eastWestEdgeFromTo(left[0], rightLowerPorts[0]);
        LPort[] rightUpperPorts = addPortsOnSide(2, right[0], PortSide.WEST);
        eastWestEdgeFromTo(left[1], rightUpperPorts[1]);
        eastWestEdgeFromTo(left[1], rightUpperPorts[0]);
        
        assertThat(allCrossings(), is(4));
    }

    @Test
    public void countCrossingsInEmptyGraph() {
        graph = getEmptyGraph();
        int amountOfCrossings = allCrossings();
        assertThat(amountOfCrossings, is(0));
    }

    @Test
    public void oneNodeIsLongEdgeDummy() {
        graph = new NorthSouthEdgeTestGraphCreator().getSouthernNorthSouthDummyEdgeCrossingGraph();
        int amountOfCrossings = allCrossings();
        assertThat(amountOfCrossings, is(1));
    }

    @Test
    public void oneNodeIsLongEdgeDummyNorthern() {
        graph = new NorthSouthEdgeTestGraphCreator().getNorthernNorthSouthDummyEdgeCrossingGraph();
        int amountOfCrossings = allCrossings();
        assertThat(amountOfCrossings, is(1));
    }

    /**
     * <pre>
     *   ----*
     *   |---*
     *   ||
     * *-++--*
     *   ||
     *  ----
     *  |__|
     *   ||
     * *-++--*
     *   ||--*
     *   |---*
     * </pre>
     */
    @Test
    public void multipleNorthSouthAndLongEdgeDummiesOnBothSides() {
        LNode[] leftNodes = addNodesToLayer(2, makeLayer());
        LNode[] middleNodes = addNodesToLayer(7, makeLayer());
        LNode[] rightNodes = addNodesToLayer(6, makeLayer());
        
        eastWestEdgeFromTo(leftNodes[0], middleNodes[2]);
        eastWestEdgeFromTo(middleNodes[2], rightNodes[2]);
        eastWestEdgeFromTo(leftNodes[1], middleNodes[4]);
        eastWestEdgeFromTo(middleNodes[4], rightNodes[3]);
        
        setAsLongEdgeDummy(middleNodes[2]);
        setAsLongEdgeDummy(middleNodes[4]);
        
        addNorthSouthEdge(PortSide.NORTH, middleNodes[3], middleNodes[0], rightNodes[0], false);
        addNorthSouthEdge(PortSide.NORTH, middleNodes[3], middleNodes[1], rightNodes[1], false);
        addNorthSouthEdge(PortSide.SOUTH, middleNodes[3], middleNodes[5], rightNodes[4], false);
        addNorthSouthEdge(PortSide.SOUTH, middleNodes[3], middleNodes[6], rightNodes[5], false);

        int amountOfCrossings = allCrossings();
        assertThat(amountOfCrossings, is(4));
    }

    /**
     * <pre>
    *   --*
    *   | ____
    *   |/|  |
    *   /\|  |
    *   | |  |
    *   | |__|
    *   |
    *    \
    *     *
     * </pre>
     */
    @Test
    public void inLayerCrossingsOnFarLeft() {
        LNode[] nodes = addNodesToLayer(3, makeLayer(graph));

        setFixedOrderConstraint(nodes[1]);

        addInLayerEdge(nodes[0], nodes[1], PortSide.WEST);
        addInLayerEdge(nodes[1], nodes[2], PortSide.WEST);

        int amountOfCrossings = allCrossings();
        assertThat(amountOfCrossings, is(1));
    }

    private int allCrossings() {
        LNode[][] nodeArray = graph.toNodeArray();
        int portId = 0;
        for (LNode[] lNodes : nodeArray) {
            for (LNode lNode : lNodes) {
                for (LPort port : lNode.getPorts()) {
                    port.id = portId;
                    portId++;
                }
            }
        }
        GraphInfoHolder gd = new GraphInfoHolder(graph, CrossMinType.BARYCENTER, null);
        return gd.crossCounter().countAllCrossings(nodeArray);
    }
}
