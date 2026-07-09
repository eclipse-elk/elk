/*******************************************************************************
 * Copyright (c) 2016 Kiel University and others.
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
package org.eclipse.elk.alg.layered.p3order;

import static org.hamcrest.CoreMatchers.is;
import static org.junit.Assert.assertThat;

import java.util.Arrays;
import java.util.List;

import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.intermediate.greedyswitch.TestGraphCreator;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.PortType;
import org.eclipse.elk.alg.layered.p3order.counting.IInitializable;
import org.eclipse.elk.core.data.LayoutMetaDataService;
import org.eclipse.elk.core.options.PortSide;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import com.google.common.collect.LinkedListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;

// CHECKSTYLEOFF javadoc
// CHECKSTYLEOFF MagicNumber
// CHECKSTYLEOFF MethodName
@SuppressWarnings("restriction")
public class BarycenterHeuristicTest extends TestGraphCreator {
    private MockRandom random;
    @Before
    public void setUp() {
        random = new MockRandom();

    }

    @BeforeClass
    public static void setup() {
        LayoutMetaDataService.initElkReflect();
    }

    /**
     * <pre>
     * *  *
     *  \/
     *  /\
     * *  *
     * .
     * </pre>
     *
     * @return Graph of the form above.
     */
    @Test
    public void minimizeCrossings_removesCrossingInSimpleCross() {
        LNode[] leftNodes = addNodesToLayer(2, makeLayer());
        LNode[] rightNodes = addNodesToLayer(2, makeLayer());
        eastWestEdgeFromTo(leftNodes[0], rightNodes[1]);
        eastWestEdgeFromTo(leftNodes[1], rightNodes[0]);
        setUpIds();
        LNode[][] nodes = getGraph().toNodeArray();

        NodeRelativePortDistributor portDist = new NodeRelativePortDistributor(nodes.length);
        ForsterConstraintResolver constraintResolver = new ForsterConstraintResolver(nodes);
        IInitializable.init(Arrays.asList(portDist, constraintResolver), nodes);
        IInitializable.init(Arrays.asList(portDist, constraintResolver), nodes);

        portDist.calculatePortRanks(nodes[0], PortType.OUTPUT);
        BarycenterHeuristic crossMin =
                new BarycenterHeuristic(constraintResolver, random, portDist, nodes);
        IInitializable.init(Arrays.asList(crossMin), nodes);

        LNode[] expectedOrder = switchOrderInArray(0, 1, nodes[1]);

        minimizeCrossings(crossMin, nodes[1], false, false, true);

        assertThat(expectedOrder, is(nodes[1]));
    }

    /**
     * <pre>
     *  
     * *  *
     *  \/
     *  /\
     * *  *
     * .
     * </pre>
     *
     * @return Graph of the form above.
     */
    @Test
    public void mockRandomizeFirstLayer() {
        LNode[] leftNodes = addNodesToLayer(2, makeLayer());
        LNode[] rightNodes = addNodesToLayer(2, makeLayer());
        eastWestEdgeFromTo(leftNodes[0], rightNodes[1]);
        eastWestEdgeFromTo(leftNodes[1], rightNodes[0]);
        setUpIds();

        LNode[][] nodes = getGraph().toNodeArray();
        NodeRelativePortDistributor portDist = new NodeRelativePortDistributor(nodes.length);
        ForsterConstraintResolver constraintResolver = new ForsterConstraintResolver(nodes);
        IInitializable.init(Arrays.asList(portDist, constraintResolver), nodes);
        portDist.calculatePortRanks(nodes[0], PortType.OUTPUT);
        BarycenterHeuristic crossMin =
                new BarycenterHeuristic(constraintResolver, random, portDist, nodes);
        IInitializable.init(Arrays.asList(crossMin), nodes);

        LNode[] expectedOrder = Arrays.copyOf(nodes[0], nodes[0].length);
        LNode[] expectedSwitchedOrder = switchOrderInArray(0, 1, nodes[0]);

        minimizeCrossings(crossMin, nodes[0], false, true, true);
        assertThat(expectedOrder, is(nodes[0]));

        random.setChangeBy(-0.01);
        minimizeCrossings(crossMin, nodes[0], false, true, true);
        assertThat(nodes[0], is(expectedSwitchedOrder));
    }

    /**
     * <pre>
     *   *  *
     *    \/
     *    /\
     * *-*  *
     * .
     * </pre>
     *
     * @return Graph of the form above.
     */
    @Test
    public void fillingInUnknownBarycenters() {
        LNode leftNode = addNodeToLayer(makeLayer());
        LNode[] middleNodes = addNodesToLayer(2, makeLayer());
        LNode[] rightNodes = addNodesToLayer(2, makeLayer());
        eastWestEdgeFromTo(middleNodes[0], rightNodes[1]);
        eastWestEdgeFromTo(middleNodes[1], rightNodes[0]);
        eastWestEdgeFromTo(leftNode, middleNodes[1]);
        setUpIds();

        LNode[][] nodes = getGraph().toNodeArray();
        LNode[] expectedSwitchedOrder = switchOrderInArray(0, 1, nodes[2]);
        LNode[] expectedOrderSecondLayer = Arrays.copyOf(nodes[1], nodes[1].length);

        NodeRelativePortDistributor portDist = new NodeRelativePortDistributor(nodes.length);
        ForsterConstraintResolver constraintResolver = new ForsterConstraintResolver(nodes);
        IInitializable.init(Arrays.asList(portDist, constraintResolver), nodes);

        BarycenterHeuristic crossMin =
                new BarycenterHeuristic(constraintResolver, random, portDist, nodes);
        IInitializable.init(Arrays.asList(crossMin), nodes);
        portDist.calculatePortRanks(nodes[0], PortType.OUTPUT);
        minimizeCrossings(crossMin, nodes[0], false, true, true);

        portDist.calculatePortRanks(nodes[1], PortType.OUTPUT);
        minimizeCrossings(crossMin, nodes[1], false, false, true);
        assertThat(nodes[1], is(expectedOrderSecondLayer));

        minimizeCrossings(crossMin, nodes[2], false, false, true);

        assertThat(nodes[2], is(expectedSwitchedOrder));
    }


    /**
     * <pre>
     * ____  *
     * |  |\/
     * |__|/\
     *       *
     * .
     * </pre>
     *
     * @return Graph of the form above.
     */
    @Test
    public void assumingFixedPortOrder_givenSimplePortOrderCross_removesCrossingIndependentOfRandom() {
        Layer leftLayer = makeLayer(getGraph());
        Layer rightLayer = makeLayer(getGraph());

        LNode leftNode = addNodeToLayer(leftLayer);

        LNode rightTopNode = addNodeToLayer(rightLayer);
        LNode rightBottomNode = addNodeToLayer(rightLayer);

        eastWestEdgeFromTo(leftNode, rightBottomNode);
        eastWestEdgeFromTo(leftNode, rightTopNode);

        setFixedOrderConstraint(leftNode);

        setUpIds();

        LNode[][] nodes = getGraph().toNodeArray();

        NodeRelativePortDistributor portDist = new NodeRelativePortDistributor(nodes.length);
        ForsterConstraintResolver constraintResolver = new ForsterConstraintResolver(nodes);
        IInitializable.init(Arrays.asList(portDist, constraintResolver), nodes);

        portDist.calculatePortRanks(nodes[0], PortType.OUTPUT);
        BarycenterHeuristic crossMin =
                new BarycenterHeuristic(constraintResolver, random, portDist, nodes);
        IInitializable.init(Arrays.asList(crossMin), nodes);

        LNode[] expectedOrder = switchOrderInArray(0, 1, nodes[1]);

        minimizeCrossings(crossMin, nodes[1], false, false, true);
        assertThat(nodes[1], is(expectedOrder));

        random.setChangeBy(-0.1);
        random.setNextBoolean(false);
        minimizeCrossings(crossMin, nodes[1], false, false, true);

        assertThat(nodes[1], is(expectedOrder));
    }

    /**
     * <pre>
     *
     * *  ___
     *  \/| |
     *  /\|_|
     * *
     * </pre>
     *
     */
    @Test
    public void assumingFixedPortOrder_givenSimplePortOrderCross_removesCrossingBackwards() {
        LNode[] leftNodes = addNodesToLayer(2, makeLayer(getGraph()));
        LNode rightNode = addNodeToLayer(makeLayer(getGraph()));
        eastWestEdgeFromTo(leftNodes[0], rightNode);
        eastWestEdgeFromTo(leftNodes[1], rightNode);
        setFixedOrderConstraint(rightNode);
        setUpIds();

        LNode[][] nodes = getGraph().toNodeArray();

        NodeRelativePortDistributor portDist = new NodeRelativePortDistributor(nodes.length);
        ForsterConstraintResolver constraintResolver = new ForsterConstraintResolver(nodes);
        IInitializable.init(Arrays.asList(portDist, constraintResolver), nodes);

        portDist.calculatePortRanks(nodes[1], PortType.INPUT);
        BarycenterHeuristic crossMin =
                new BarycenterHeuristic(constraintResolver, random, portDist, nodes);
        IInitializable.init(Arrays.asList(crossMin), nodes);

        LNode[] expectedOrder = switchOrderInArray(0, 1, nodes[0]);

        minimizeCrossings(crossMin, nodes[0], false, false, false);

        assertThat(expectedOrder, is(nodes[0]));
    }

    /**
     * <pre>
     *       ___
     *    ---| |
     *    |  | |  <- switch this
     * ---+--|_|
     * |  |
     * *--|--*  <- with this
     *    |
     *    ---*
     * .
     * </pre>
     *
     * With fixed Port PortOrder.
     *
     * @return Graph of the form above.
     */
    @Test
    public void inLayerEdges() {
        LNode leftNode = addNodeToLayer(makeLayer());
        LNode[] rightNodes = addNodesToLayer(3, makeLayer());
        setFixedOrderConstraint(rightNodes[0]);
        eastWestEdgeFromTo(leftNode, rightNodes[0]);
        addInLayerEdge(rightNodes[0], rightNodes[2], PortSide.WEST);
        eastWestEdgeFromTo(leftNode, rightNodes[1]);
        setUpIds();
        LNode[][] nodes = getGraph().toNodeArray();

        NodeRelativePortDistributor portDist = new NodeRelativePortDistributor(nodes.length);
        ForsterConstraintResolver constraintResolver = new ForsterConstraintResolver(nodes);
        IInitializable.init(Arrays.asList(portDist, constraintResolver), nodes);

        portDist.calculatePortRanks(nodes[0], PortType.INPUT);
        BarycenterHeuristic crossMin =
                new BarycenterHeuristic(constraintResolver, random, portDist, nodes);
        IInitializable.init(Arrays.asList(crossMin), nodes);

        LNode[] expectedOrder = getArrayInIndexOrder(nodes[1], 2, 0, 1);

        minimizeCrossings(crossMin, nodes[1], false, false, true);

        assertThat(nodes[1], is(expectedOrder));
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
     * </pre>
     * @return Graph of the form above.
     */
    @Test
    public void northSouthEdges() {
        LNode[] leftNodes = addNodesToLayer(1, makeLayer());
        LNode[] middleNodes = addNodesToLayer(4, makeLayer());
        LNode[] rightNodes = addNodesToLayer(3, makeLayer());
        eastWestEdgeFromTo(leftNodes[0], middleNodes[2]);
        eastWestEdgeFromTo(middleNodes[2], rightNodes[2]);
        setAsLongEdgeDummy(middleNodes[2]);
        addNorthSouthEdge(PortSide.NORTH, middleNodes[3], middleNodes[0], rightNodes[0], false);
        addNorthSouthEdge(PortSide.NORTH, middleNodes[3], middleNodes[1], rightNodes[1], false);
        setUpIds();

        Multimap<LNode, LNode> layoutUnits = LinkedListMultimap.create();
        LNode[][] nodes = getGraph().toNodeArray();
        for (LNode[] list : nodes) {
            for (LNode node : list) {
                LNode layoutUnit = node.getProperty(InternalProperties.IN_LAYER_LAYOUT_UNIT);
                if (layoutUnit != null) {
                    layoutUnits.put(layoutUnit, node);
                }
            }
        }

        NodeRelativePortDistributor portDist = new NodeRelativePortDistributor(nodes.length);
        ForsterConstraintResolver constraintResolver = new ForsterConstraintResolver(nodes);
        IInitializable.init(Arrays.asList(portDist, constraintResolver), nodes);

        portDist.calculatePortRanks(nodes[0], PortType.INPUT);

        BarycenterHeuristic crossMin =
                new BarycenterHeuristic(constraintResolver, random, portDist, nodes);
        IInitializable.init(Arrays.asList(crossMin), nodes);
        LNode[] expectedOrder = getArrayInIndexOrder(nodes[1], 1, 0, 3, 2);

        random.setChangeBy(-0.01);
        minimizeCrossings(crossMin, nodes[1], false, false, true);

        assertThat(nodes[1], is(expectedOrder));

    }

    /**
     * Helper method that wraps the nodes array in a list and applies the order back to the array.
     */
    private void minimizeCrossings(final BarycenterHeuristic crossMin, final LNode[] nodes,
            final boolean preOrdered, final boolean randomized, final boolean forward) {
        List<LNode> nodeList = Lists.newArrayList(nodes);
        crossMin.minimizeCrossings(nodeList, preOrdered, randomized, forward);
        int i = 0;
        for (LNode node : nodeList) {
            nodes[i++] = node;
        }
    }
}
