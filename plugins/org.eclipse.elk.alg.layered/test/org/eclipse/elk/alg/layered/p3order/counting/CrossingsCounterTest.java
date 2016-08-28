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
package org.eclipse.elk.alg.layered.p3order.counting;

import static org.hamcrest.CoreMatchers.is;
import static org.junit.Assert.assertThat;

import java.util.Arrays;
import java.util.Date;
import java.util.Random;

import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.intermediate.greedyswitch.InLayerEdgeTestGraphCreator;
import org.eclipse.elk.alg.layered.p3order.GraphData;
import org.eclipse.elk.alg.layered.p3order.LayerSweepCrossingMinimizer.CrossMinType;
import org.eclipse.elk.core.options.PortSide;
import org.junit.Ignore;
import org.junit.Test;

// CHECKSTYLEOFF MagicNumber
// CHECKSTYLEOFF javadoc
// CHECKSTYLEOFF MethodName
public class CrossingsCounterTest extends InLayerEdgeTestGraphCreator {
    private CrossingsCounter counter;
    private LNode[][] order;


    /**
     * <pre>
     * ___  ___
     * | |\/| |
     * |_|/\|_|
     * </pre>
     */
    @Test
    public void countCrossingsBetweenLayers_fixedPortOrderCrossingOnTwoNodes() {
        LNode left = addNodeToLayer(makeLayer(getGraph()));
        LNode right = addNodeToLayer(makeLayer(getGraph()));
        eastWestEdgeFromTo(left, right);
        eastWestEdgeFromTo(left, right);

        counter = new CrossingsCounter(new int[getNumPorts(order())]);

        assertThat(counter.countCrossingsBetweenLayers(order()[0], order()[1]), is(1));
    }

    /**
     * <pre>
     * *
     *  \
     *  /
     * *
     *  \
     * *+--
     *  | |
     * * /
     * *
     * </pre>
     * 
     * @throws Exception
     */
    @Test
    public void longInLayerCrossings() throws Exception {
        LNode[] nodes = addNodesToLayer(5, makeLayer());
        addInLayerEdge(nodes[0], nodes[1], PortSide.EAST);
        addInLayerEdge(nodes[1], nodes[3], PortSide.EAST);
        addInLayerEdge(nodes[2], nodes[4], PortSide.EAST);

        counter = new CrossingsCounter(new int[getNumPorts(order())]);

        assertThat(counter.countInLayerCrossingsOnSide(order()[0], PortSide.EAST), is(1));
    }

    @Test
    public void countCrossingsBetweenLayers_crossFormed() {
        getCrossFormedGraph();

        counter = new CrossingsCounter(new int[getNumPorts(order())]);

        assertThat(counter.countCrossingsBetweenLayers(order()[0], order()[1]), is(1));
    }

    private LNode[][] order() {
        return getGraph().toNodeArray();
    }

    /**
     * Constructs a cross formed graph with two edges between the corners
     *
     * <pre>
     * *    *
     *  \\//
     *  //\\
     * *    *
     * .
     * </pre>
     */
    @Test
    public void countCrossingsBetweenLayers_crossFormedMultipleEdgesBetweenSameNodes() {
        Layer leftLayer = makeLayer(graph);
        Layer rightLayer = makeLayer(graph);

        LNode topLeft = addNodeToLayer(leftLayer);
        LNode bottomLeft = addNodeToLayer(leftLayer);
        LNode topRight = addNodeToLayer(rightLayer);
        LNode bottomRight = addNodeToLayer(rightLayer);

        LPort topLeftTopPort = addPortOnSide(topLeft, PortSide.EAST);
        LPort topLeftBottomPort = addPortOnSide(topLeft, PortSide.EAST);
        LPort bottomRightBottomPort = addPortOnSide(bottomRight, PortSide.WEST);
        LPort bottomRightTopPort = addPortOnSide(bottomRight, PortSide.WEST);
        addEdgeBetweenPorts(topLeftTopPort, bottomRightTopPort);
        addEdgeBetweenPorts(topLeftBottomPort, bottomRightBottomPort);

        LPort bottomLeftTopPort = addPortOnSide(bottomLeft, PortSide.EAST);
        LPort bottomLeftBottomPort = addPortOnSide(bottomLeft, PortSide.EAST);
        LPort topRightBottomPort = addPortOnSide(topRight, PortSide.WEST);
        LPort topRightTopPort = addPortOnSide(topRight, PortSide.WEST);
        addEdgeBetweenPorts(bottomLeftTopPort, topRightTopPort);
        addEdgeBetweenPorts(bottomLeftBottomPort, topRightBottomPort);

        GraphData gd = new GraphData(graph, CrossMinType.BARYCENTER, null);
        gd.portDistributor().distributePortsWhileSweeping(order(), 1, true);

        counter = new CrossingsCounter(new int[getNumPorts(order())]);

        assertThat(counter.countCrossingsBetweenLayers(order()[0], order()[1]), is(4));
    }

    @Test
    public void countCrossingsBetweenLayers_crossWithExtraEdgeInBetween() {
        getCrossWithExtraEdgeInBetweenGraph();

        counter = new CrossingsCounter(new int[getNumPorts(order())]);

        assertThat(counter.countCrossingsBetweenLayers(order()[0], order()[1]), is(3));
    }

    @Test
    public void countCrossingsBetweenLayers_ignoreSelfLoops() {
        getCrossWithManySelfLoopsGraph();

        counter = new CrossingsCounter(new int[getNumPorts(order())]);

        assertThat(counter.countCrossingsBetweenLayers(order()[0], order()[1]), is(1));
    }

    @Test
    public void countCrossingsBetweenLayers_moreComplexThreeLayerGraph() {
        getMoreComplexThreeLayerGraph();
        GraphData gd = new GraphData(graph, CrossMinType.BARYCENTER, null);
        gd.portDistributor().distributePortsWhileSweeping(order(), 1, true);
        counter = new CrossingsCounter(new int[getNumPorts(order())]);
        assertThat(counter.countCrossingsBetweenLayers(order()[0], order()[1]), is(1));
    }

    @Test
    public void countCrossingsBetweenLayers_fixedPortOrder() {
        getFixedPortOrderGraph();

        counter = new CrossingsCounter(new int[getNumPorts(order())]);

        assertThat(counter.countCrossingsBetweenLayers(order()[0], order()[1]), is(1));
    }

    /**
     * <pre>
     * *   *<- Into same port
     *  \//
     *  //\
     * *   *
     * </pre>
     */
    @Test
    public void countCrossingsBetweenLayers_intoSamePort() {
        Layer leftLayer = makeLayer(graph);
        Layer rightLayer = makeLayer(graph);
        
        LNode topLeft = addNodeToLayer(leftLayer);
        LNode bottomLeft = addNodeToLayer(leftLayer);
        LNode topRight = addNodeToLayer(rightLayer);
        LNode bottomRight = addNodeToLayer(rightLayer);
        
        eastWestEdgeFromTo(topLeft, bottomRight);
        LPort bottomLeftFirstPort = addPortOnSide(bottomLeft, PortSide.EAST);
        LPort bottomLeftSecondPort = addPortOnSide(bottomLeft, PortSide.EAST);
        LPort topRightFirstPort = addPortOnSide(topRight, PortSide.WEST);
        
        addEdgeBetweenPorts(bottomLeftFirstPort, topRightFirstPort);
        addEdgeBetweenPorts(bottomLeftSecondPort, topRightFirstPort);
        setUpIds();
        counter = new CrossingsCounter(new int[getNumPorts(order())]);

        assertThat(counter.countCrossingsBetweenLayers(order()[0], order()[1]), is(2));
    }

    /**
     * <pre>
     * *   /*
     * |  /
     * \ /____
     *  x/|  |
     * |/\|  |
     * *  |__|
     * </pre>
     */
    @Test
    public void countCrossingsBetweenPorts_givenWesternCrossings_OnlyCountsForGivenPorts() {
        LNode[] leftNodes = addNodesToLayer(2, makeLayer(getGraph()));
        LNode[] rightNodes = addNodesToLayer(2, makeLayer(getGraph()));
        eastWestEdgeFromTo(leftNodes[0], rightNodes[1]);
        eastWestEdgeFromTo(leftNodes[1], rightNodes[1]);
        eastWestEdgeFromTo(leftNodes[1], rightNodes[0]);

        counter = new CrossingsCounter(new int[getNumPorts(order())]);
        counter.initForCountingBetweenOnSide(leftNodes, rightNodes, PortSide.WEST);
        assertThat(counter.countCrossingsBetweenPorts(rightNodes[1].getPorts().get(1), rightNodes[1].getPorts().get(0)),
                is(1));
    }

    /**
     * <pre>
     * ___
     * | |\/*
     * |_|/\*
     * </pre>
     */
    @Test
    public void countCrossingsBetweenPorts_GivenCrossingsOnEasternSide_() throws Exception {
        LNode[] leftNodes = addNodesToLayer(1, makeLayer());
        LNode[] rightNodes = addNodesToLayer(2, makeLayer());
        eastWestEdgeFromTo(leftNodes[0], rightNodes[1]);
        eastWestEdgeFromTo(leftNodes[0], rightNodes[0]);

        counter = new CrossingsCounter(new int[getNumPorts(order())]);
        counter.initForCountingBetweenOnSide(leftNodes, rightNodes, PortSide.EAST);
        assertThat(counter.countCrossingsBetweenPorts(leftNodes[0].getPorts().get(0), leftNodes[0].getPorts().get(1)),
                is(1));
    }

    /**
     * <pre>
     * *   *
     *  \//
     *  //\
     * *   *
     * ^Into same port
     * </pre>
     */
    @Test
    public void countCrossingsBetweenPorts_twoEdgesIntoSamePort() {
        Layer leftLayer = makeLayer();
        Layer rightLayer = makeLayer();
        
        LNode topLeft = addNodeToLayer(leftLayer);
        LNode bottomLeft = addNodeToLayer(leftLayer);
        LNode topRight = addNodeToLayer(rightLayer);
        LNode bottomRight = addNodeToLayer(rightLayer);
        
        eastWestEdgeFromTo(topLeft, bottomRight);
        LPort bottomLeftPort = addPortOnSide(bottomLeft, PortSide.EAST);
        LPort topRightPort = addPortOnSide(topRight, PortSide.WEST);
        
        addEdgeBetweenPorts(bottomLeftPort, topRightPort);
        addEdgeBetweenPorts(bottomLeftPort, topRightPort);
        setUpIds();
        
        counter = new CrossingsCounter(new int[getNumPorts(order())]);
        counter.initForCountingBetweenOnSide(order()[0], order()[1], PortSide.WEST);
        
        assertThat(counter.countCrossingsBetweenPorts(bottomLeftPort, topLeft.getPorts().get(0)), is(2));
    }

    @Ignore
    // @Test
    public void benchmark() {
        makeTwoLayerRandomGraphWithNodesPerLayer(6000, 6);

        counter = new CrossingsCounter(new int[getNumPorts(order())]);
        System.out.println("Starting");
        int length = 400;
        long[] times = new long[length];
        for (int i = 0; i < length; i++) {
            long tick = new Date().getTime();
            counter.countCrossingsBetweenLayers(order()[0], order()[1]);
            times[i] = new Date().getTime() - tick;
        }
        System.out.println(Arrays.stream(times).min());
    }

    private void makeTwoLayerRandomGraphWithNodesPerLayer(final int numNodes, final int edgesPerNode) {
        LNode[] leftNodes = addNodesToLayer(numNodes, makeLayer());
        LNode[] rightNodes = addNodesToLayer(numNodes, makeLayer());
        Random random = new Random(0);
        for (int i = 0; i < edgesPerNode * numNodes; i++) {
            if (random.nextBoolean()) {
                LNode left = leftNodes[random.nextInt(numNodes)];
                LNode right = rightNodes[random.nextInt(numNodes)];
                eastWestEdgeFromTo(left, right);
            } else {
                addInLayerEdge(leftNodes[random.nextInt(numNodes)], leftNodes[random.nextInt(numNodes)], PortSide.EAST);
            }
        }
        for (LNode node : rightNodes) {
            node.cachePortSides();
        }
        for (LNode node : leftNodes) {
            node.cachePortSides();
        }
    }

    private int getNumPorts(final LNode[][] currentOrder) {
        int numPorts = 0;
        for (LNode[] lNodes : currentOrder) {
            for (LNode node : lNodes) {
                numPorts += node.getPorts().size();
            }
        }
        return numPorts;
    }

}
