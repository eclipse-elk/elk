/*******************************************************************************
 * Copyright (c) 2016 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.alg.layered.p3order;

import static org.hamcrest.CoreMatchers.is;
import static org.junit.Assert.assertThat;

import java.util.Arrays;
import java.util.List;

import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.intermediate.greedyswitch.TestGraphCreator;
import org.eclipse.elk.alg.layered.p3order.LayerSweepCrossingMinimizer.CrossMinType;
import org.eclipse.elk.alg.layered.p3order.counting.IInitializable;
import org.eclipse.elk.core.data.LayoutMetaDataService;
import org.eclipse.elk.core.options.PortSide;
import org.junit.BeforeClass;
import org.junit.Ignore;
import org.junit.Test;

import com.google.common.collect.Lists;

@SuppressWarnings("restriction")
public class AbstractBarycenterPortDistributorTest extends TestGraphCreator {
    
    @BeforeClass
    public static void initialize() {
        LayoutMetaDataService.getInstance();
    }
    
    // CHECKSTYLEOFF javadoc
    // CHECKSTYLEOFF MagicNumber
    // CHECKSTYLEOFF MethodName
    /**
     * <pre>
     * *  ___
     *  \/| |
     *  /\| |
     * *  |_|
     * </pre>
     */
    @Test
    public void distributePortsOnSide_GivenCrossOnWesternSide_ShouldRemoveCrossing() {
        LNode[] leftNodes = addNodesToLayer(2, makeLayer(getGraph()));
        LNode rightNode = addNodeToLayer(makeLayer(getGraph()));
        eastWestEdgeFromTo(leftNodes[0], rightNode);
        eastWestEdgeFromTo(leftNodes[1], rightNode);

        List<LPort> expectedPortOrderRightNode =
                Lists.newArrayList(rightNode.getPorts().get(1), rightNode.getPorts().get(0));

        distributePortsInCompleteGraph(4);

        assertThat(rightNode.getPorts(), is(expectedPortOrderRightNode));
    }

    /**
     * <pre>
     * *  ___  *
     *  \/| |\/
     *  /\| |/\
     * *  |_|  *
     * </pre>
     *
     */
    @Test
    public void distributePortsOfGraph_GivenCrossOnBothSides_ShouldRemoveCrossin() {
        LNode[] leftNodes = addNodesToLayer(2, makeLayer(getGraph()));
        LNode middleNode = addNodeToLayer(makeLayer(getGraph()));
        LNode[] rightNodes = addNodesToLayer(2, makeLayer(getGraph()));
        eastWestEdgeFromTo(middleNode, rightNodes[1]);
        eastWestEdgeFromTo(middleNode, rightNodes[0]);
        eastWestEdgeFromTo(leftNodes[0], middleNode);
        eastWestEdgeFromTo(leftNodes[1], middleNode);
        setUpIds();
        List<LPort> expectedPortOrderMiddleNode = copyPortsInIndexOrder(middleNode, 1, 0, 3, 2);

        distributePortsInCompleteGraph(8);

        assertThat(middleNode.getPorts(), is(expectedPortOrderMiddleNode));
    }

    private void distributePortsInCompleteGraph(final int numberOfPorts) {
        GraphInfoHolder gd = new GraphInfoHolder(graph, CrossMinType.BARYCENTER, null);
        LNode[][] nodes = graph.toNodeArray();
        for (int i = 0; i < nodes.length; i++) {
            gd.portDistributor().distributePortsWhileSweeping(nodes, i, true);
        }
        for (int i = nodes.length - 1; i >= 0; i--) {
            gd.portDistributor().distributePortsWhileSweeping(nodes, i, false);
        }
    }

    /**
     * <pre>
     * ___
     * | |\ /-*
     * | | x
     * |_|/ \-*
     * </pre>
     */
    @Test
    public void distributePortsOfGraph_GivenCrossOnEasternSide_ShouldRemoveCrossing() {
        LNode leftNode = addNodeToLayer(makeLayer(getGraph()));
        LNode[] rightNodes = addNodesToLayer(2, makeLayer(getGraph()));
        eastWestEdgeFromTo(leftNode, rightNodes[1]);
        eastWestEdgeFromTo(leftNode, rightNodes[0]);

        List<LPort> expectedPortOrderLeftNode = copyPortsInIndexOrder(leftNode, 1, 0);

        distributePortsInCompleteGraph(4);

        assertThat(leftNode.getPorts(), is(expectedPortOrderLeftNode));
    }

    /**
     * <pre>
     *     *-----
     *     *-\  |
     *   ____ | |
     * * |  |-+--
     *   |__|-|
     * </pre>
     */
    @Test
    public void distributePortsOfGraph_GivenInLayerEdgePortOrderCrossing_ShouldRemoveIt() {
        addNodeToLayer(makeLayer());
        LNode[] nodes = addNodesToLayer(3, makeLayer());
        addInLayerEdge(nodes[0], nodes[2], PortSide.EAST);
        addInLayerEdge(nodes[1], nodes[2], PortSide.EAST);

        List<LPort> expectedPortOrderLowerNode = copyPortsInIndexOrder(nodes[2], 1, 0);

        distributePortsInCompleteGraph(4);

        assertThat(nodes[2].getPorts(), is(expectedPortOrderLowerNode));
    }

    /**
     * <pre>
     *     *-->*
     *     |
     *   *-+-->*
     *   | |
     *  _|_|_
     *  |   |
     *  |___|
     *  .
     * </pre>
     */
    @Test
    public void distributePortsOfGraph_GivenNorthSouthPortOrderCrossing_ShouldSwitchPortOrder() {
        LNode[] leftNodes = addNodesToLayer(3, makeLayer(getGraph()));
        LNode[] rightNodes = addNodesToLayer(2, makeLayer(getGraph()));

        addNorthSouthEdge(PortSide.NORTH, leftNodes[2], leftNodes[1], rightNodes[1], false);
        addNorthSouthEdge(PortSide.NORTH, leftNodes[2], leftNodes[0], rightNodes[0], false);

        List<LPort> expectedPortOrderLowerNode =
                Lists.newArrayList(leftNodes[2].getPorts().get(1), leftNodes[2].getPorts().get(0));

        distributePortsInCompleteGraph(6);

        assertThat(leftNodes[2].getPorts(), is(expectedPortOrderLowerNode));
    }

    /**
     * <pre>
     * ___  ____
     * | |\/|  |
     * |_|/\|  |
     *      |--|
     * </pre>
     */
    @Test
    public void distributePortsWhileSweeping_givenSimpleCross_ShouldRemoveCrossing() {
        LNode leftNode = addNodeToLayer(makeLayer(getGraph()));
        LNode rightNode = addNodeToLayer(makeLayer(getGraph()));
        eastWestEdgeFromTo(leftNode, rightNode);
        eastWestEdgeFromTo(leftNode, rightNode);
        List<LPort> expectedPortRightNode = copyPortsInIndexOrder(rightNode, 1, 0);
        setUpIds();
        LNode[][] nodeArray = graph.toNodeArray();
        ISweepPortDistributor portDist = new LayerTotalPortDistributor(nodeArray.length);
        IInitializable.init(Arrays.asList(portDist), nodeArray);
        portDist.distributePortsWhileSweeping(nodeArray, 1, true);

        assertThat(rightNode.getPorts(), is(expectedPortRightNode));
    }

    /**
     * <pre>
     * ____
     * | *+--  *
     * |  |  \/
     * |  |\ /\
     * | *+-x  *
     * |__|  \
     *        -*
     * </pre>
     */
    // TODO this is a problem which currently cannot be solved by our algorithm :-(
    @Ignore
    public void distributePortsOnSide_partlyCrossHierarchicalEdges_CrossHierarchyStaysOuterChanges() {
        LNode leftOuterNode = addNodeToLayer(makeLayer());
        LNode[] rightNodes = addNodesToLayer(3, makeLayer());
        LPort[] leftOuterPorts = addPortsOnSide(3, leftOuterNode, PortSide.EAST);
        LGraph leftInnerGraph = nestedGraph(leftOuterNode);
        LNode[] leftInnerNodes = addNodesToLayer(2, makeLayer(leftInnerGraph));
        LNode[] leftInnerDummyNodes = new LNode[2];
        Layer dummyLayer = makeLayer();
        leftInnerDummyNodes[0] = addExternalPortDummyNodeToLayer(dummyLayer, leftOuterPorts[0]);
        leftInnerDummyNodes[1] = addExternalPortDummyNodeToLayer(dummyLayer, leftOuterPorts[2]);
        eastWestEdgeFromTo(leftInnerNodes[0], leftInnerDummyNodes[0]);
        eastWestEdgeFromTo(leftInnerNodes[1], leftInnerDummyNodes[1]);
        eastWestEdgeFromTo(leftOuterPorts[0], rightNodes[1]);
        eastWestEdgeFromTo(leftOuterPorts[1], rightNodes[2]);
        eastWestEdgeFromTo(leftOuterPorts[2], rightNodes[0]);
        // leftOuterNode.setProperty(InternalProperties.HAS_HIERARCHICAL_AND_NORMAL_PORTS, true);
        setPortOrderFixed(leftOuterNode);

        setUpIds();

        List<LPort> expectedOrder = Lists.newArrayList(switchOrderInArray(1, 2, leftOuterPorts));

        distributePortsInCompleteGraph(8);

        assertThat(leftOuterNode.getPorts(), is(expectedOrder));
    }

}
