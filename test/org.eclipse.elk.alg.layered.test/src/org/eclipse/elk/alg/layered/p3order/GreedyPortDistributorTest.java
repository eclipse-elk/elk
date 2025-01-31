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
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertThat;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.intermediate.greedyswitch.TestGraphCreator;
import org.eclipse.elk.alg.layered.p3order.counting.IInitializable;
import org.eclipse.elk.core.options.PortSide;
import org.junit.Test;

@SuppressWarnings("restriction")
public class GreedyPortDistributorTest extends TestGraphCreator {
    private GreedyPortDistributor portDist;

    // CHECKSTYLEOFF javadoc
    // CHECKSTYLEOFF MagicNumber
    // CHECKSTYLEOFF MethodName

    private void setUpDistributor() {
        portDist = new GreedyPortDistributor();
        LNode[][] nodeArray = graph.toNodeArray();
        IInitializable.init(Arrays.asList(portDist), nodeArray);
    }

    /**
     * <pre>
     * *  ___
     *  \/| |
     *  /\| |
     * *  |_|
     * </pre>
     */
    @Test
    public void distributePorts_GivenCrossOnWesternSide_RemoveCrossing() {
        LNode[] leftNodes = addNodesToLayer(2, makeLayer(getGraph()));
        LNode rightNode = addNodeToLayer(makeLayer(getGraph()));
        eastWestEdgeFromTo(leftNodes[0], rightNode);
        eastWestEdgeFromTo(leftNodes[1], rightNode);

        List<LPort> expectedPortOrderRightNode = portsOrderedAs(rightNode, 1, 0);

        setUpDistributor();
        boolean improved = portDist.distributePortsWhileSweeping(graph.toNodeArray(), 1, true);

        assertTrue(improved);
        assertThat(rightNode.getPorts(), is(expectedPortOrderRightNode));
    }

    /**
     * <pre>
     * *  ___
     *  \/| | *
     *  /\| | *
     * *  |_|
     * </pre>
     */
    @Test
    public void distributePorts_GivenNoPortsOnRightSide_NothingHappens() {
        LNode[] leftNodes = addNodesToLayer(2, makeLayer(getGraph()));
        LNode middleNode = addNodeToLayer(makeLayer(getGraph()));
        addNodesToLayer(2, makeLayer(getGraph()));
        eastWestEdgeFromTo(leftNodes[0], middleNode);
        eastWestEdgeFromTo(leftNodes[1], middleNode);

        List<LPort> expectedPortOrderRightNode = portsOrderedAs(middleNode, 0, 1);

        setUpDistributor();
        portDist.distributePortsWhileSweeping(graph.toNodeArray(), 1, false);

        assertThat(middleNode.getPorts(), is(expectedPortOrderRightNode));
    }

    /**
     * <pre>
     * *    ___
     *  \/--| |
     *  /\ /| |
     * *  x | |
     * *-/ \|_|
     * </pre>
     */
    @Test
    public void distributePorts_GivenMultipleCrossingsOnWesternSide_RemoveCrossing() {
        LNode[] leftNodes = addNodesToLayer(3, makeLayer(getGraph()));
        LNode rightNode = addNodeToLayer(makeLayer(getGraph()));
        eastWestEdgeFromTo(leftNodes[0], rightNode);
        eastWestEdgeFromTo(leftNodes[2], rightNode);
        eastWestEdgeFromTo(leftNodes[1], rightNode);

        List<LPort> expectedPortOrderRightNode = portsOrderedAs(rightNode, 1, 2, 0);

        setUpDistributor();

        portDist.distributePortsWhileSweeping(graph.toNodeArray(), 1, true);

        assertThat(rightNode.getPorts(), is(expectedPortOrderRightNode));
    }

    /**
     * <pre>
     * ___
     * | |\/*
     * |_|/\*
     * </pre>
     */
    @Test
    public void distributePorts_GivenCrossingsOnEasternSide_RemoveThem() throws Exception {
        LNode[] leftNodes = addNodesToLayer(1, makeLayer());
        LNode[] rightNodes = addNodesToLayer(2, makeLayer());
        eastWestEdgeFromTo(leftNodes[0], rightNodes[1]);
        eastWestEdgeFromTo(leftNodes[0], rightNodes[0]);

        List<LPort> expectedPortOrderLeftNode = portsOrderedAs(leftNodes[0], 1, 0);

        setUpDistributor();
        portDist.distributePortsWhileSweeping(graph.toNodeArray(), 0, false);

        assertThat(leftNodes[0].getPorts(), is(expectedPortOrderLeftNode));
    }

    /**
     * <pre>
     * ___
     * | |\/*
     * |_|/\*
     * </pre>
     */
    @Test
    public void distributePorts_fixedPortOrder_NoChange() throws Exception {
        LNode[] leftNodes = addNodesToLayer(1, makeLayer());
        LNode[] rightNodes = addNodesToLayer(2, makeLayer());
        eastWestEdgeFromTo(leftNodes[0], rightNodes[1]);
        eastWestEdgeFromTo(leftNodes[0], rightNodes[0]);
        setFixedOrderConstraint(leftNodes[0]);

        List<LPort> expectedPortOrderLeftNode = portsOrderedAs(leftNodes[0], 0, 1);

        setUpDistributor();
        portDist.distributePortsWhileSweeping(graph.toNodeArray(), 0, false);

        assertThat(leftNodes[0].getPorts(), is(expectedPortOrderLeftNode));
    }

    /**
     * <pre>
     * ____
     * |*-+   *
     * |  |\\/
     * |*-+/\\
     * |--|   *
     * </pre>
     */
    @Test
    public void givenDoubleCrossBetweenCompoundAndNonCompoundNodes_SwitchesPorts() {
        LNode leftOuterNode = addNodeToLayer(makeLayer(getGraph()));
        LNode[] rightNodes = addNodesToLayer(2, makeLayer(getGraph()));
        LPort[] leftOuterPorts = addPortsOnSide(2, leftOuterNode, PortSide.EAST);
        eastWestEdgeFromTo(leftOuterPorts[0], rightNodes[1]);
        eastWestEdgeFromTo(leftOuterPorts[0], rightNodes[1]);
        eastWestEdgeFromTo(leftOuterPorts[1], rightNodes[0]);
        LGraph leftInnerGraph = nestedGraph(leftOuterNode);
        LNode[] leftInnerNodes = addNodesToLayer(2, makeLayer(leftInnerGraph));
        LNode[] leftInnerDummyNodes = addExternalPortDummiesToLayer(makeLayer(leftInnerGraph), leftOuterPorts);
        eastWestEdgeFromTo(leftInnerNodes[0], leftInnerDummyNodes[0]);
        eastWestEdgeFromTo(leftInnerNodes[1], leftInnerDummyNodes[1]);
        setUpIds();

        List<LPort> expectedPortOrderLeftNode = portsOrderedAs(leftOuterNode, 1, 0);

        setUpDistributor();
        portDist.distributePortsWhileSweeping(graph.toNodeArray(), 0, false);

        assertThat(leftOuterNode.getPorts(), is(expectedPortOrderLeftNode));
    }

    /**
     * <pre>
     * ____
     * |*-+  *
     * |  |\/
     * |*-+/\
     * |--|  *
     * </pre>
     */
    @Test
    public void givenSingleCrossBetweenCompoundAndNonCompoundNodes_DoesNotSwitchPorts() {
        LNode leftOuterNode = addNodeToLayer(makeLayer(getGraph()));
        LNode[] rightNodes = addNodesToLayer(2, makeLayer(getGraph()));
        LPort[] leftOuterPorts = addPortsOnSide(2, leftOuterNode, PortSide.EAST);
        eastWestEdgeFromTo(leftOuterPorts[0], rightNodes[1]);
        eastWestEdgeFromTo(leftOuterPorts[1], rightNodes[0]);
        LGraph leftInnerGraph = nestedGraph(leftOuterNode);
        LNode[] leftInnerNodes = addNodesToLayer(2, makeLayer(leftInnerGraph));
        LNode[] leftInnerDummyNodes = addExternalPortDummiesToLayer(makeLayer(leftInnerGraph), leftOuterPorts);
        eastWestEdgeFromTo(leftInnerNodes[0], leftInnerDummyNodes[0]);
        eastWestEdgeFromTo(leftInnerNodes[1], leftInnerDummyNodes[1]);

        // Order stays the same.
        List<LPort> expectedPortOrderLeftNode = portsOrderedAs(leftOuterNode, 0, 1);

        setUpDistributor();
        portDist.distributePortsWhileSweeping(graph.toNodeArray(), 0, false);

        assertThat(leftOuterNode.getPorts(), is(expectedPortOrderLeftNode));
    }

    /**
     * <pre>
     * ____
     * |*-+  *
     * |  |\/
     * |*-+/\
     * |  |  *
     * |*-+--*
     * |--|
     * </pre>
     */
    @Test
    public void givenMoreHierarchicalNodes_DoesNotSwitchPorts() {
        LNode leftOuterNode = addNodeToLayer(makeLayer(getGraph()));
        LNode[] rightNodes = addNodesToLayer(3, makeLayer(getGraph()));
        LPort[] leftOuterPorts = addPortsOnSide(3, leftOuterNode, PortSide.EAST);
        eastWestEdgeFromTo(leftOuterPorts[0], rightNodes[1]);
        eastWestEdgeFromTo(leftOuterPorts[1], rightNodes[0]);
        eastWestEdgeFromTo(leftOuterPorts[2], rightNodes[2]);
        LGraph leftInnerGraph = nestedGraph(leftOuterNode);
        LNode[] leftInnerNodes = addNodesToLayer(3, makeLayer(leftInnerGraph));
        LNode[] leftInnerDummyNodes = addExternalPortDummiesToLayer(makeLayer(leftInnerGraph), leftOuterPorts);
        eastWestEdgeFromTo(leftInnerNodes[0], leftInnerDummyNodes[0]);
        eastWestEdgeFromTo(leftInnerNodes[1], leftInnerDummyNodes[1]);
        eastWestEdgeFromTo(leftInnerNodes[2], leftInnerDummyNodes[2]);

        // Order stays the same.
        List<LPort> expectedPortOrderLeftNode = portsOrderedAs(leftOuterNode, 0, 1, 2);

        setUpDistributor();
        portDist.distributePortsWhileSweeping(graph.toNodeArray(), 0, false);

        assertThat(leftOuterNode.getPorts(), is(expectedPortOrderLeftNode));
    }

    /**
     * <pre>
     * _______
     * |   *-+  *
     * |     |\/
     * |   *-+/\
     * |     |  *
     * |*--*-+--*
     * |-----|
     * </pre>
     */
    @Test
    public void givenMoreHierarchicalNodes2_DoesNotSwitchPorts() {
        LNode leftOuterNode = addNodeToLayer(makeLayer(getGraph()));
        LNode[] rightNodes = addNodesToLayer(3, makeLayer(getGraph()));
        LPort[] leftOuterPorts = addPortsOnSide(3, leftOuterNode, PortSide.EAST);
        eastWestEdgeFromTo(leftOuterPorts[0], rightNodes[1]);
        eastWestEdgeFromTo(leftOuterPorts[1], rightNodes[0]);
        eastWestEdgeFromTo(leftOuterPorts[2], rightNodes[2]);
        LGraph leftInnerGraph = nestedGraph(leftOuterNode);
        LNode leftInnerNode = addNodeToLayer(makeLayer(leftInnerGraph));
        LNode[] rightInnerNodes = addNodesToLayer(3, makeLayer(leftInnerGraph));
        LNode[] dummyNodes = addExternalPortDummiesToLayer(makeLayer(leftInnerGraph), leftOuterPorts);
        eastWestEdgeFromTo(leftInnerNode, rightInnerNodes[2]);
        eastWestEdgeFromTo(rightInnerNodes[0], dummyNodes[0]);
        eastWestEdgeFromTo(rightInnerNodes[1], dummyNodes[1]);
        eastWestEdgeFromTo(rightInnerNodes[2], dummyNodes[2]);

        // Order stays the same.
        List<LPort> expectedPortOrderLeftNode = portsOrderedAs(leftOuterNode, 0, 1, 2);

        setUpDistributor();
        portDist.distributePortsWhileSweeping(graph.toNodeArray(), 0, false);

        assertThat(leftOuterNode.getPorts(), is(expectedPortOrderLeftNode));
    }

    /**
     * <pre>
     * ___
     * | |\/*
     * |_|/\*
     * ___
     * | |\/*
     * |_|/\*
     * </pre>
     */
    @Test
    public void distributePortsWhileSweeping_givenTwoHierarchicalNodesInOneLayer() throws Exception {
        Layer leftLayer = makeLayer(getGraph());
        Layer rightLayer = makeLayer(getGraph());

        for (int i = 0; i < 2; i++) {
            LNode[] leftNodes = addNodesToLayer(1, leftLayer);
            LNode[] rightNodes = addNodesToLayer(2, rightLayer);
            eastWestEdgeFromTo(leftNodes[0], rightNodes[1]);
            eastWestEdgeFromTo(leftNodes[0], rightNodes[0]);
        }

        // Order stays the same.
        LNode leftOuterNode = leftLayer.getNodes().get(0);
        List<LPort> expectedPortOrderLeftNode = portsOrderedAs(leftOuterNode, 1, 0);

        LNode[][] nodeOrder = getGraph().toNodeArray();
        setUpDistributor();
        portDist.distributePortsWhileSweeping(nodeOrder, 0, false);

        assertThat(leftOuterNode.getPorts(), is(expectedPortOrderLeftNode));
    }

    /**
     * <pre>
     * ___
     * | |--*
     * |_|--*
     * </pre>
     * 
     * @throws Exception
     */
    @Test
    public void noChange() throws Exception {
        Layer leftLayer = makeLayer(getGraph());
        Layer rightLayer = makeLayer(getGraph());
        LNode[] leftNodes = addNodesToLayer(1, leftLayer);
        LNode[] rightNodes = addNodesToLayer(2, rightLayer);
        eastWestEdgeFromTo(leftNodes[0], rightNodes[0]);
        eastWestEdgeFromTo(leftNodes[0], rightNodes[1]);
        LNode[][] nodeOrder = getGraph().toNodeArray();
        setUpDistributor();
        boolean improved = portDist.distributePortsWhileSweeping(nodeOrder, 0, false);

        assertFalse(improved);
    }

    private List<LPort> portsOrderedAs(final LNode rightNode, final int... indices) {
        List<LPort> ordered = new ArrayList<>();
        for (int i : indices) {
            ordered.add(rightNode.getPorts().get(i));
        }
        return ordered;
    }
}
