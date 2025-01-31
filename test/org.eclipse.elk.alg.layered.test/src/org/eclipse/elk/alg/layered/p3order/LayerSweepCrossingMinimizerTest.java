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
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertThat;
import static org.junit.Assert.assertTrue;

import java.util.Arrays;
import java.util.List;
import java.util.Random;

import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.intermediate.greedyswitch.TestGraphCreator;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.layered.p3order.LayerSweepCrossingMinimizer.CrossMinType;
import org.eclipse.elk.alg.test.PlainJavaInitialization;
import org.eclipse.elk.core.options.PortConstraints;
import org.eclipse.elk.core.options.PortSide;
import org.eclipse.elk.core.util.BasicProgressMonitor;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Ignore;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameters;

import com.google.common.collect.Lists;

// CHECKSTYLEOFF javadoc
// CHECKSTYLEOFF MagicNumber
// CHECKSTYLEOFF MethodName
@SuppressWarnings("restriction")
@RunWith(Parameterized.class)
public class LayerSweepCrossingMinimizerTest extends TestGraphCreator {
    private final CrossMinType crossMinType;
    private LayerSweepCrossingMinimizer crossMin;

    @BeforeClass
    public static void init() {
        PlainJavaInitialization.initializePlainJavaLayout();
    }

    /**
     * Constructor called by Parameterized.
     *
     * @param gT
     *            greedyType
     */
    public LayerSweepCrossingMinimizerTest(final CrossMinType crossMinType) {
        this.crossMinType = crossMinType;
    }

    /**
     * Sets the Parameters to be tested as all elements of the GreedySwitchType enum.
     *
     * @return parameters
     */
    @Parameters(name = "{0}")
    public static Iterable<Object[]> greedyTypes() {
        return Arrays.asList(
                new Object[][] { { CrossMinType.BARYCENTER }, { CrossMinType.ONE_SIDED_GREEDY_SWITCH } });
    }

    @Before
    public void setUp() {
        crossMin = new LayerSweepCrossingMinimizer(crossMinType);
    }

    @Test
    public void givenOneNode_ShouldntCrash() {
        getOneNodeGraph();

        setUpAndMinimizeCrossings();
        // should cause no errors
    }

    @Test
    public void givenNoNode_ShouldntCrash() {

        setUpAndMinimizeCrossings();
        // should cause no errors
    }

    /**
     * <pre>
     * ________
     * |___  *|
     * || |\/ |
     * || |/\ |
     * ||_|  *|
     * |------|
     * </pre>
     *
     * Sweep backward first.
     */
    @Test
    public void givenCompoundGraphWhereOrderIsOnlyCorrectedOnForwardSweep_RemovesCrossing() {
        LNode node = addNodeToLayer(makeLayer(getGraph()));
        LGraph innerGraph = nestedGraph(node);
        LNode leftInnerNode = addNodeToLayer(makeLayer(innerGraph));
        setFixedOrderConstraint(leftInnerNode);
        LNode[] rightInnerNodes = addNodesToLayer(2, makeLayer(innerGraph));
        eastWestEdgeFromTo(leftInnerNode, rightInnerNodes[1]);
        eastWestEdgeFromTo(leftInnerNode, rightInnerNodes[0]);

        Layer rightInnerLayer = innerGraph.getLayers().get(1);
        List<LNode> expectedOrderRightInner = switchOrderOfNodesInLayer(0, 1, rightInnerLayer);
        List<LPort> expectedPortOrderLeft = copyPortsInIndexOrder(leftInnerNode, 0, 1);

        getRandom().setNextBoolean(false);
        setUpAndMinimizeCrossings();

        assertThat(leftInnerNode.getPorts(), is(expectedPortOrderLeft));
        assertThat(rightInnerLayer.getNodes(), is(expectedOrderRightInner));
    }

    /**
     * <pre>
     * ______
     * |*  *|
     * | \/ |
     * | /\ |
     * |*  *|
     * |----|
     * </pre>
     */
    @Test
    public void givenSingleHierarchicalNodeWithCross_RemovesCrossing() {
        LNode node = addNodeToLayer(makeLayer(getGraph()));
        LGraph innerGraph = nestedGraph(node);
        LNode[] innerNodesleft = addNodesToLayer(2, makeLayer(innerGraph));
        LNode[] innerNodesRight = addNodesToLayer(2, makeLayer(innerGraph));
        eastWestEdgeFromTo(innerNodesleft[0], innerNodesRight[1]);
        eastWestEdgeFromTo(innerNodesleft[1], innerNodesRight[0]);

        List<LNode> expectedOrderRightBaryCenter =
                switchOrderOfNodesInLayer(0, 1, innerGraph.getLayers().get(1));

        setUpAndMinimizeCrossings();

        assertThat(innerGraph.getLayers().get(1).getNodes(), is(expectedOrderRightBaryCenter));
    }

    /**
     * <pre>
     * ____  ____
     * |*-+  +-*|
     * |  |\/|  |
     * |*-+/\+-*|
     * |--|  |--|
     * </pre>
     */
    @Test
    public void givenSimpleHierarchicalCross_ShouldResultInNoCrossing() {
        LNode leftOuterNode = addNodeToLayer(makeLayer(getGraph()));
        LNode rightOuterNode = addNodeToLayer(makeLayer(getGraph()));
        LPort[] leftOuterPorts = addPortsOnSide(2, leftOuterNode, PortSide.EAST);
        LPort[] rightOuterPorts = addPortsOnSide(2, rightOuterNode, PortSide.WEST);

        addEdgeBetweenPorts(leftOuterPorts[0], rightOuterPorts[0]);
        addEdgeBetweenPorts(leftOuterPorts[1], rightOuterPorts[1]);

        makeNestedTwoNodeGraphWithEasternPorts(leftOuterNode, leftOuterPorts);

        LGraph rightInnerGraph =
                makeNestedTwoNodeGraphWithWesternPorts(rightOuterNode, rightOuterPorts);

        List<LNode> expectedExternalDummyOrderRight =
                switchOrderOfNodesInLayer(0, 1, rightInnerGraph.getLayers().get(0));
        List<LNode> expectedNormalNodeOrderRight =
                switchOrderOfNodesInLayer(0, 1, rightInnerGraph.getLayers().get(1));

        List<LPort> expectedOrderOfPortsRight =
                Lists.newArrayList(rightOuterPorts[1], rightOuterPorts[0]);

        setOnAllGraphs(LayeredOptions.CROSSING_MINIMIZATION_HIERARCHICAL_SWEEPINESS, 0.1, graph);
        setUpAndMinimizeCrossings();
        List<LNode> actualExternalDummyOrderRight = rightInnerGraph.getLayers().get(0).getNodes();
        assertThat(actualExternalDummyOrderRight, is(expectedExternalDummyOrderRight));
        assertThat(rightOuterNode.getPorts(), is(expectedOrderOfPortsRight));

        List<LNode> actualNormalOrderRight = rightInnerGraph.getLayers().get(1).getNodes();
        assertThat(actualNormalOrderRight, is(expectedNormalNodeOrderRight));
    }

    /**
     * <pre>
     * ____  ____
     * |*-+  +-*|
     * |  |\/|  |
     * |*-+/\+-*|
     * |--|  |--|
     * </pre>
     */
    @Test
    public void givenSimpleHierarchicalCrossSweepingFromRightToLeft_ShouldResultInNoCrossing() {
        LNode leftOuterNode = addNodeToLayer(makeLayer(getGraph()));
        LNode rightOuterNode = addNodeToLayer(makeLayer(getGraph()));
        LPort[] leftOuterPorts = addPortsOnSide(2, leftOuterNode, PortSide.EAST);
        LPort[] rightOuterPorts = addPortsOnSide(2, rightOuterNode, PortSide.WEST);
        addEdgeBetweenPorts(leftOuterPorts[0], rightOuterPorts[0]);
        addEdgeBetweenPorts(leftOuterPorts[1], rightOuterPorts[1]);

        LGraph leftInnerGraph =
                makeNestedTwoNodeGraphWithEasternPorts(leftOuterNode, leftOuterPorts);

        makeNestedTwoNodeGraphWithWesternPorts(rightOuterNode, rightOuterPorts);

        List<LNode> expectedNodeOrderLeft =
                switchOrderOfNodesInLayer(0, 1, leftInnerGraph.getLayers().get(0));

        List<LPort> expectedOrderOfPortsLeft =
                Lists.newArrayList(leftOuterPorts[1], leftOuterPorts[0]);

        getRandom().setNextBoolean(false);
        setUpAndMinimizeCrossings();

        assertThat(leftOuterNode.getPorts(), is(expectedOrderOfPortsLeft));
        assertThat(leftInnerGraph.getLayers().get(0).getNodes(), is(expectedNodeOrderLeft));

    }

    /**
     * <pre>
     * ___        ___
     * | |===*  --| |
     * | |    \/  | |
     * | |    /\  | |
     * |_|===*  --|_|
     *       ___  ___
     *       | |\/| |
     *       |_|/\|_|
     * .           ^
     *           port orders fixed
     * </pre>
     *
     * First Layer and last layer in fixed order.
     *
     * @return graph of the form above.
     */
    @Test
    public void althoughBackwardSweepNotTakenStillCorrectsPortOrder() {
        Layer[] layers = makeLayers(3);
        LNode leftNode = addNodeToLayer(layers[0]);
        LNode[] middleNodes = addNodesToLayer(3, layers[1]);
        LNode[] rightNodes = addNodesToLayer(2, layers[2]);

        eastWestEdgeFromTo(middleNodes[0], rightNodes[0]);
        eastWestEdgeFromTo(middleNodes[1], rightNodes[0]);

        eastWestEdgeFromTo(middleNodes[2], rightNodes[1]);
        eastWestEdgeFromTo(middleNodes[2], rightNodes[1]);
        eastWestEdgeFromTo(leftNode, middleNodes[0]);
        eastWestEdgeFromTo(leftNode, middleNodes[0]);
        eastWestEdgeFromTo(leftNode, middleNodes[1]);
        eastWestEdgeFromTo(leftNode, middleNodes[1]);

        setFixedOrderConstraint(rightNodes[1]);
        setFixedOrderConstraint(rightNodes[0]);
        setFixedOrderConstraint(leftNode);

        setUpIds();
        List<LPort> expectedPortOrderBottomMiddleNode = copyPortsInIndexOrder(middleNodes[2], 1, 0);

        setUpAndMinimizeCrossings();

        Layer middleLayer = getGraph().getLayers().get(1);
        assertTrue(middleLayer.toString(), inOrder(middleNodes[0], middleNodes[1], middleLayer));
        assertThat(middleNodes[2].getPorts(), is(expectedPortOrderBottomMiddleNode));
    }

    /**
     * <pre>
     * ___        ___
     * | |---*  --| |
     * | |    \/  | |
     * | |    /\  | |
     * |_|---*  --|_|
     *       ___  ___
     *       | |\/| |
     *       |_|/\|_|
     *             ^
     *           port orders fixed
     * </pre>
     *
     * First Layer and last layer in fixed order. Should leave this order, since the fixable port
     * order is already corrected before counting.
     *
     * @return graph of the form above.
     */
    @Test
    public void doesNotCountRemovablePortCrossing() {
        Layer[] layers = makeLayers(3);
        LNode leftNode = addNodeToLayer(layers[0]);
        LNode[] middleNodes = addNodesToLayer(3, layers[1]);
        LNode[] rightNodes = addNodesToLayer(2, layers[2]);

        eastWestEdgeFromTo(middleNodes[0], rightNodes[0]);
        eastWestEdgeFromTo(middleNodes[1], rightNodes[0]);

        eastWestEdgeFromTo(middleNodes[2], rightNodes[1]);
        eastWestEdgeFromTo(middleNodes[2], rightNodes[1]);
        eastWestEdgeFromTo(leftNode, middleNodes[0]);
        eastWestEdgeFromTo(leftNode, middleNodes[1]);

        setFixedOrderConstraint(rightNodes[1]);
        setFixedOrderConstraint(rightNodes[0]);
        setFixedOrderConstraint(leftNode);

        setUpIds();

        List<LPort> expectedPortOrderBottomMiddleNode = copyPortsInIndexOrder(middleNodes[2], 1, 0);

        setUpAndMinimizeCrossings();

        Layer middleLayer = getGraph().getLayers().get(1);
        assertTrue(middleLayer.toString(), inOrder(middleNodes[0], middleNodes[1], middleLayer));
        assertThat(middleNodes[2].getPorts(), is(expectedPortOrderBottomMiddleNode));
    }

    /**
     * <pre>
     * *   ___
     *  \  | |
     * *-+-| |
     *   |/|_|
     *   |\
     *   --*
     * </pre>
     */
    @Test
    public void resolvesInLayerPortOrderCrossingsAfterSwitch() {
        LNode[] leftNodes = addNodesToLayer(2, makeLayer(getGraph()));
        LNode[] middleNodes = addNodesToLayer(2, makeLayer(getGraph()));

        eastWestEdgeFromTo(leftNodes[0], middleNodes[1]);
        addInLayerEdge(middleNodes[1], middleNodes[0], PortSide.WEST);
        eastWestEdgeFromTo(leftNodes[1], middleNodes[0]);
        setUpIds();

        List<LNode> expectedNodeOrderMiddleLayer = copyOfSwitchOrderOfNodesInLayer(0, 1, 1);
        List<LPort> expectedPortOrderMiddleTopNode = copyPortsInIndexOrder(middleNodes[0], 1, 0);

        setUpAndMinimizeCrossings();

        if (crossMinType == CrossMinType.BARYCENTER) {
            assertThat(getGraph().getLayers().get(1).getNodes(), is(expectedNodeOrderMiddleLayer));
            assertThat(middleNodes[0].getPorts(), is(expectedPortOrderMiddleTopNode));
        }
    }

    /**
     * Assumes both nodes to actually be in layer.
     *
     * @param firstNode
     * @param secondNode
     * @param layer
     * @return
     */
    private boolean inOrder(final LNode firstNode, final LNode secondNode, final Layer layer) {
        for (LNode node : layer) {
            if (node == firstNode) {
                return true;
            } else if (node == secondNode) {
                return false;
            }
        }
        return false;
    }

    private LGraph makeNestedTwoNodeGraphWithWesternPorts(final LNode rightOuterNode,
            final LPort[] rightOuterPorts) {
        LGraph rightInnerGraph = nestedGraph(rightOuterNode);
        LNode[] rightInnerDummyNodes =
                addExternalPortDummiesToLayer(makeLayer(rightInnerGraph), rightOuterPorts);
        LNode[] rightInnerNodes = addNodesToLayer(2, makeLayer(rightInnerGraph));
        eastWestEdgeFromTo(rightInnerDummyNodes[0], rightInnerNodes[0]);
        eastWestEdgeFromTo(rightInnerDummyNodes[1], rightInnerNodes[1]);
        return rightInnerGraph;
    }

    private LGraph makeNestedTwoNodeGraphWithEasternPorts(final LNode leftOuterNode,
            final LPort[] leftOuterPorts) {
        LGraph leftInnerGraph = nestedGraph(leftOuterNode);
        LNode[] leftInnerNodes = addNodesToLayer(2, makeLayer(leftInnerGraph));
        LNode[] leftInnerDummyNodes =
                addExternalPortDummiesToLayer(makeLayer(leftInnerGraph), leftOuterPorts);
        eastWestEdgeFromTo(leftInnerNodes[0], leftInnerDummyNodes[0]);
        eastWestEdgeFromTo(leftInnerNodes[1], leftInnerDummyNodes[1]);
        return leftInnerGraph;
    }

    /**
     * <pre>
     * ____  ____
     * |*-+  | *|
     * |  |\/|  |
     * |*-+/\| *|
     * |--|  |--|
     * </pre>
     */
    @Test
    public void givenCrossWithNoExternalPortDummiesOnOneNestedGraph_ShouldRemoveCrossing() {
        LNode leftOuterNode = addNodeToLayer(makeLayer());
        LNode rightOuterNode = addNodeToLayer(makeLayer());
        LPort[] leftOuterPorts = addPortsOnSide(2, leftOuterNode, PortSide.EAST);
        LPort[] rightOuterPorts = addPortsOnSide(2, rightOuterNode, PortSide.WEST);
        addEdgeBetweenPorts(leftOuterPorts[0], rightOuterPorts[0]);
        addEdgeBetweenPorts(leftOuterPorts[1], rightOuterPorts[1]);
        makeNestedTwoNodeGraphWithEasternPorts(leftOuterNode, leftOuterPorts);
        LGraph rightInnerGraph = nestedGraph(rightOuterNode);
        addNodesToLayer(2, makeLayer(rightInnerGraph));

        List<LPort> expectedPortOrderRight =
                Lists.newArrayList(rightOuterPorts[1], rightOuterPorts[0]);

        setUpAndMinimizeCrossings();

        assertThat(expectedPortOrderRight, is(rightOuterNode.getPorts()));
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
    public void givenCrossBetweenCompoundAndNonCompoundNodes_ShouldRemoveCrossing() {
        LNode leftOuterNode = addNodeToLayer(makeLayer(getGraph()));
        LNode[] rightNodes = addNodesToLayer(2, makeLayer(getGraph()));
        LPort[] leftOuterPorts = addPortsOnSide(2, leftOuterNode, PortSide.EAST);
        eastWestEdgeFromTo(leftOuterPorts[0], rightNodes[1]);
        eastWestEdgeFromTo(leftOuterPorts[1], rightNodes[0]);
        makeNestedTwoNodeGraphWithEasternPorts(leftOuterNode, leftOuterPorts);

        List<LNode> expectedOrderRight = copyOfSwitchOrderOfNodesInLayer(0, 1, 1);

        setUpAndMinimizeCrossings();

        assertThat(expectedOrderRight, is(getGraph().getLayers().get(1).getNodes()));
    }

    /**
     * <pre>
     * ______
     * |*  *+-*
     * | \/ |
     * | /\ |
     * |*  *+-*
     * |----|
     * </pre>
     */
    @Test
    public void givenCrossInFirstLevelCompoundNode_ShouldRemoveCrossing() {
        // parent graph
        LNode leftOuterNode = addNodeToLayer(makeLayer(getGraph()));
        LNode[] rightNodes = addNodesToLayer(2, makeLayer(getGraph()));
        LPort[] leftOuterPorts = addPortsOnSide(2, leftOuterNode, PortSide.EAST);
        eastWestEdgeFromTo(leftOuterPorts[0], rightNodes[0]);
        eastWestEdgeFromTo(leftOuterPorts[1], rightNodes[1]);

        // child graph
        LGraph leftInnerGraph = nestedGraph(leftOuterNode);
        LNode[] leftInnerNodesleft = addNodesToLayer(2, makeLayer(leftInnerGraph));
        LNode[] leftInnerNodesRight = addNodesToLayer(2, makeLayer(leftInnerGraph));
        LNode[] leftInnerDummyNodes =
                addExternalPortDummiesToLayer(makeLayer(leftInnerGraph), leftOuterPorts);
        eastWestEdgeFromTo(leftInnerNodesRight[0], leftInnerDummyNodes[0]);
        eastWestEdgeFromTo(leftInnerNodesRight[1], leftInnerDummyNodes[1]);
        eastWestEdgeFromTo(leftInnerNodesleft[0], leftInnerNodesRight[1]);
        eastWestEdgeFromTo(leftInnerNodesleft[1], leftInnerNodesRight[0]);

        List<LNode> expectedOrderRight = copyOfSwitchOrderOfNodesInLayer(0, 1, 1);

        setUpAndMinimizeCrossings();

        assertThat(getGraph().getLayers().get(1).getNodes(), is(expectedOrderRight));
    }

    /**
     * <pre>
     *    ___
     * *  | |
     *  \/| |
     *  /\| |
     * *  |_|
     * </pre>
     *
     * port order fixed.
     */
    @Test
    public void givenGraphWithoutNesting_ShouldImproveOnBackwardSweep() {
        Layer leftLayer = makeLayer(getGraph());
        Layer rightLayer = makeLayer(getGraph());
        LNode leftTopNode = addNodeToLayer(leftLayer);
        LNode leftBottomNode = addNodeToLayer(leftLayer);
        LNode rightNode = addNodeToLayer(rightLayer);
        rightNode.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_ORDER);
        eastWestEdgeFromTo(leftTopNode, rightNode);
        eastWestEdgeFromTo(leftBottomNode, rightNode);

        List<LNode> expectedOrderLayerOne = copyOfSwitchOrderOfNodesInLayer(0, 1, 0);

        getRandom().setNextBoolean(false); // sweeps backward
        setUpAndMinimizeCrossings();

        assertThat("Layer one", copyOfNodesInLayer(0), is(expectedOrderLayerOne));
    }

    /**
     * <pre>
     *    *
     *   /
     *  /
     * *  *
     *  \/
     *  /\
     * *  *
     *
     * </pre>
     */
    @Test
    public void givenSimpleGraph_ShouldNotBeReorderedRandomlyOnBackwardSweep() {
        Layer[] layers = makeLayers(2);
        LNode[] leftNodes = addNodesToLayer(2, layers[0]);
        LNode[] rightNodes = addNodesToLayer(3, layers[1]);
        eastWestEdgeFromTo(leftNodes[0], rightNodes[0]);
        eastWestEdgeFromTo(leftNodes[0], rightNodes[2]);
        eastWestEdgeFromTo(leftNodes[1], rightNodes[1]);

        List<LNode> expectedOrderLayerTwo = copyOfSwitchOrderOfNodesInLayer(1, 2, 1);
        List<LNode> expectedOrderLayerOne = copyOfNodesInLayer(0);

        setUpAndMinimizeCrossings();

        assertThat("Layer two", copyOfNodesInLayer(1), is(expectedOrderLayerTwo));
        assertThat("Layer one", copyOfNodesInLayer(0), is(expectedOrderLayerOne));

    }

    /**
     * <pre>
     * ___     ___
     * | |\ /*-| |
     * | |\x/  | |
     * | |/x\  | |
     * |_|/ \*-|_|
     * </pre>
     *
     * port order fixed.
     */
    @Test
    public void givenGraphWhichWorsensOnBackwardSweep_ShouldTakeResultOfForwardSweep() {
        Layer[] layers = makeLayers(3);
        LNode leftNode = addNodeToLayer(layers[0]);
        setFixedOrderConstraint(leftNode);
        LNode[] middleNodes = addNodesToLayer(2, layers[1]);
        LNode rightNode = addNodeToLayer(layers[2]);
        setFixedOrderConstraint(rightNode);
        eastWestEdgeFromTo(leftNode, middleNodes[1]);
        eastWestEdgeFromTo(leftNode, middleNodes[1]);
        eastWestEdgeFromTo(leftNode, middleNodes[0]);
        eastWestEdgeFromTo(leftNode, middleNodes[0]);
        eastWestEdgeFromTo(middleNodes[1], rightNode);
        eastWestEdgeFromTo(middleNodes[0], rightNode);

        List<LNode> expectedOrderLayerTwo = copyOfSwitchOrderOfNodesInLayer(0, 1, 1);

        setUpAndMinimizeCrossings();

        assertThat("Layer two", copyOfNodesInLayer(1), is(expectedOrderLayerTwo));
    }

    /**
     * <pre>
     * ____  ____
     * |*-+  +-*|
     * |  |\/|  |
     * |*-+/\+-*|
     * |--|  |--|
     * </pre>
     *
     * With external origins of port dummies on right node wrong way around.
     */
    @Test
    public void givenNestedGraphWithWronglySortedDummyNodes_ShouldSortAndResolveCrossing() {
        LNode leftOuterNode = addNodeToLayer(makeLayer(getGraph()));
        LNode rightOuterNode = addNodeToLayer(makeLayer(getGraph()));
        LPort[] leftOuterPorts = addPortsOnSide(2, leftOuterNode, PortSide.EAST);
        LPort[] rightOuterPorts = addPortsOnSide(2, rightOuterNode, PortSide.WEST);
        addEdgeBetweenPorts(leftOuterPorts[0], rightOuterPorts[0]);
        addEdgeBetweenPorts(leftOuterPorts[1], rightOuterPorts[1]);

        makeNestedTwoNodeGraphWithEasternPorts(leftOuterNode, leftOuterPorts);

        LGraph rightInnerGraph = nestedGraph(rightOuterNode);
        LPort[] reversedRightPorts = reverse(rightOuterPorts);
        LNode[] rightInnerDummyNodes =
                addExternalPortDummiesToLayer(makeLayer(rightInnerGraph), reversedRightPorts);
        LNode[] rightInnerNodes = addNodesToLayer(2, makeLayer(rightInnerGraph));
        eastWestEdgeFromTo(rightInnerDummyNodes[1], rightInnerNodes[0]);
        eastWestEdgeFromTo(rightInnerDummyNodes[0], rightInnerNodes[1]);

        List<LNode> expectedNormalNodeOrderRight =
                switchOrderOfNodesInLayer(0, 1, rightInnerGraph.getLayers().get(1));

        List<LPort> expectedOrderOfPortsRight =
                Lists.newArrayList(rightOuterPorts[1], rightOuterPorts[0]);

        setUpAndMinimizeCrossings();
        List<LNode> actualExternalDummyOrderRight = rightInnerGraph.getLayers().get(0).getNodes();
        assertThat(actualExternalDummyOrderRight,
                is((List<LNode>) Lists.newArrayList(rightInnerDummyNodes)));

        assertThat(rightOuterNode.getPorts(), is(expectedOrderOfPortsRight));

        List<LNode> actualNormalOrderRight = rightInnerGraph.getLayers().get(1).getNodes();
        assertThat(actualNormalOrderRight, is(expectedNormalNodeOrderRight));
    }

    /**
     * <pre>
     *       /*
     *      / ____
     * *   *==|  |
     * \\ //  |  |
     *  xx    |  |
     * // \\  |  |
     * *   *==|__|
     * </pre>
     *
     * This test was used in development and can be reused if there is a problem with the
     * thoroughness value. It is deactivated because it depends on a random value and is slow.
     */
    @Ignore
    public void givenGraphWhichCanBeRepairedGivenGoodRandomStart_ShouldFindBestSolution() {
        Layer[] layers = makeLayers(3);
        LNode[] leftNodes = addNodesToLayer(2, layers[0]);
        LNode[] middleNodes = addNodesToLayer(2, layers[1]);
        LNode[] rightNodes = addNodesToLayer(3, layers[2]);
        setFixedOrderConstraint(rightNodes[1]);
        eastWestEdgesFromTo(2, leftNodes[0], middleNodes[1]);
        eastWestEdgesFromTo(2, leftNodes[1], middleNodes[0]);
        eastWestEdgeFromTo(middleNodes[0], rightNodes[0]);
        eastWestEdgesFromTo(2, middleNodes[1], rightNodes[1]);
        eastWestEdgesFromTo(2, middleNodes[0], rightNodes[1]);

        List<LNode> expectedWrongOrder = copyOfNodesInLayer(0);
        List<LNode> expectedCorrectOrder = copyOfSwitchOrderOfNodesInLayer(0, 1, 0);

        setUpAndMinimizeCrossings();

        assertThat(getGraph().getLayers().get(0).getNodes(), is(expectedWrongOrder));

        getGraph().setProperty(InternalProperties.RANDOM, new Random()); // AAA! Random in Test!!!
        getGraph().setProperty(LayeredOptions.THOROUGHNESS, 100);
        for (int i = 0; i < 100; i++) {
            setUpAndMinimizeCrossings();
            assertThat(getGraph().getLayers().get(0).getNodes(), is(expectedCorrectOrder));
            getGraph().getLayers().get(0).getNodes().clear();
            getGraph().getLayers().get(0).getNodes().addAll(expectedWrongOrder);
        }
    }

    /**
     * <pre>
     * ____
     * |  |\ /-*
     * |  | x
     * |*-+/ \-*
     * |--|
     * </pre>
     */
    @Test
    public void givenGraphWithNormalEdgeAndHierarchicalEdgeCrossing_ShouldRemoveCrossing() {
        Layer[] layers = makeLayers(2);
        LNode leftOuterNode = addNodeToLayer(layers[0]);
        LNode[] rightNodes = addNodesToLayer(2, layers[1]);
        eastWestEdgeFromTo(leftOuterNode, rightNodes[1]);

        Layer[] innerLayers = makeLayers(2, nestedGraph(leftOuterNode));
        LNode innerNormalNode = addNodeToLayer(innerLayers[0]);
        LPort hierarchPort = addPortOnSide(leftOuterNode, PortSide.EAST);
        LNode innerDummyNode = addExternalPortDummyNodeToLayer(innerLayers[1], hierarchPort);
        eastWestEdgeFromTo(innerNormalNode, innerDummyNode);

        LPort rightTopPort = addPortOnSide(rightNodes[0], PortSide.WEST);
        addEdgeBetweenPorts(hierarchPort, rightTopPort);

        List<LNode> expectedCorrectOrder = copyOfSwitchOrderOfNodesInLayer(0, 1, 1);

        setUpAndMinimizeCrossings();

        assertThat(getGraph().getLayers().get(1).getNodes(), is(expectedCorrectOrder));
    }

    /**
     * <pre>
     * ____
     * |  |\ /-*
     * |  | x
     * |*-+/ \-*
     * |  p
     * |--|
     * </pre>
     *
     * p is a port with no edge.
     */
    @Test
    public void givenGraphWithPortWithoutEdge_ShouldRemoveCrossing() {
        Layer[] layers = makeLayers(2);
        LNode leftOuterNode = addNodeToLayer(layers[0]);
        LNode[] rightNodes = addNodesToLayer(2, layers[1]);
        eastWestEdgeFromTo(leftOuterNode, rightNodes[1]);

        Layer[] innerLayers = makeLayers(2, nestedGraph(leftOuterNode));
        LNode innerNormalNode = addNodeToLayer(innerLayers[0]);
        LPort hierarchPort = addPortOnSide(leftOuterNode, PortSide.EAST);
        LNode innerDummyNode = addExternalPortDummyNodeToLayer(innerLayers[1], hierarchPort);
        eastWestEdgeFromTo(innerNormalNode, innerDummyNode);

        addPortOnSide(leftOuterNode, PortSide.EAST);

        LPort rightTopPort = addPortOnSide(rightNodes[0], PortSide.WEST);
        addEdgeBetweenPorts(hierarchPort, rightTopPort);

        List<LNode> expectedCorrectOrder = copyOfSwitchOrderOfNodesInLayer(0, 1, 1);

        setUpAndMinimizeCrossings();
        assertThat(getGraph().getLayers().get(1).getNodes(), is(expectedCorrectOrder));

    }

    /**
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
    public void showsErrorInCrossingCountingAlgorithmsForFreePortOrder() {
        LNode leftNode = addNodeToLayer(makeLayer(getGraph()));
        LNode[] rightNodes = addNodesToLayer(2, makeLayer(getGraph()));
        setFixedOrderConstraint(rightNodes[0]);
        setFixedOrderConstraint(rightNodes[1]);
        // setFixedOrderConstraint(leftNode);
        LPort leftTopPort = addPortOnSide(leftNode, PortSide.EAST);
        LPort leftLowerPort = addPortOnSide(leftNode, PortSide.EAST);
        LPort rightUpperNodeLowerPort = addPortOnSide(rightNodes[0], PortSide.WEST);
        LPort rightUpperNodeTopPort = addPortOnSide(rightNodes[0], PortSide.WEST);
        LPort rightLowerPort = addPortOnSide(rightNodes[1], PortSide.WEST);
        addEdgeBetweenPorts(leftTopPort, rightUpperNodeTopPort);
        addEdgeBetweenPorts(leftTopPort, rightLowerPort);
        addEdgeBetweenPorts(leftLowerPort, rightUpperNodeLowerPort);

        List<LNode> expectedCorrectOrder = copyOfSwitchOrderOfNodesInLayer(0, 1, 1);

        setUpAndMinimizeCrossings();

        assertThat(getGraph().getLayers().get(1).getNodes(), is(expectedCorrectOrder));
    }

    /**
     * <pre>
     * ___  ____
     * | |\/+-*|
     * |_|/\|  |
     *      |--|
     * </pre>
     */
    @Test
    public void nonRecursive_givenGraphWithCrossingWithFixedPortOrderToHierarchicalGraph_ShouldRemoveCrossing() {
        LNode leftNode = addNodeToLayer(makeLayer(getGraph()));
        setFixedOrderConstraint(leftNode);
        LNode rightOuterNode = addNodeToLayer(makeLayer(getGraph()));
        eastWestEdgeFromTo(leftNode, rightOuterNode);
        LPort rightOuterPort = addPortOnSide(rightOuterNode, PortSide.WEST);
        LPort leftLowerPort = addPortOnSide(leftNode, PortSide.EAST);
        addEdgeBetweenPorts(leftLowerPort, rightOuterPort);

        LGraph rightInnerGraph = nestedGraph(rightOuterNode);
        LNode rightInnerDummyNode =
                addExternalPortDummyNodeToLayer(makeLayer(rightInnerGraph), rightOuterPort);
        LNode rightInnerNode = addNodeToLayer(makeLayer(rightInnerGraph));
        eastWestEdgeFromTo(rightInnerDummyNode, rightInnerNode);

        List<LPort> expectedPortOrderRight = Lists.newArrayList(rightOuterNode.getPorts().get(1),
                rightOuterNode.getPorts().get(0));

        setOnAllGraphs(LayeredOptions.CROSSING_MINIMIZATION_HIERARCHICAL_SWEEPINESS, 1.0, graph);

        setUpAndMinimizeCrossings();
        assertThat(getGraph().getLayers().get(1).getNodes().get(0).getPorts(),
                    is(expectedPortOrderRight));
    }

    /**
     * <pre>
     * ____  _____
     * |  |  +-*-+-*
     * |  |\/|   |
     * |  |/\+-*-+-*
     * |--|  |---|
     * </pre>
     */
    @Test
    public void givenGraphWhichMustTransportSwitchedPortOrderThroughHierarchy_ShouldContainNoCrossing() {
        LNode leftNode = addNodeToLayer(makeLayer(getGraph()));
         setFixedOrderConstraint(leftNode);
        LNode middleOuterNode = addNodeToLayer(makeLayer(getGraph()));
        LNode[] rightNodes = addNodesToLayer(2, makeLayer(getGraph()));

        LPort[] middleOuterRightPorts = addPortsOnSide(2, middleOuterNode, PortSide.EAST);
        eastWestEdgeFromTo(middleOuterRightPorts[0], rightNodes[0]);
        eastWestEdgeFromTo(middleOuterRightPorts[1], rightNodes[1]);

        LPort[] leftPorts = addPortsOnSide(2, leftNode, PortSide.EAST);
        LPort[] middleOuterLeftPorts = addPortsOnSide(2, middleOuterNode, PortSide.WEST);
        addEdgeBetweenPorts(leftPorts[0], middleOuterLeftPorts[0]);
        addEdgeBetweenPorts(leftPorts[1], middleOuterLeftPorts[1]);

        LGraph innerGraph = nestedGraph(middleOuterNode);
        LNode[] leftInnerDummyNodes =
                addExternalPortDummiesToLayer(makeLayer(innerGraph), middleOuterLeftPorts);
        LNode[] innerNodes = addNodesToLayer(2, makeLayer(innerGraph));
        LNode[] rightInnerDummyNodes =
                addExternalPortDummiesToLayer(makeLayer(innerGraph), middleOuterRightPorts);

        eastWestEdgeFromTo(innerNodes[0], rightInnerDummyNodes[0]);
        eastWestEdgeFromTo(innerNodes[1], rightInnerDummyNodes[1]);
        eastWestEdgeFromTo(leftInnerDummyNodes[0], innerNodes[0]);
        eastWestEdgeFromTo(leftInnerDummyNodes[1], innerNodes[1]);

        List<LPort> expextedPorts = copyPortsInIndexOrder(middleOuterNode, 1, 0, 3, 2);

        setUpAndMinimizeCrossings();

        assertThat(getGraph().getLayers().get(1).getNodes().get(0).getPorts(), is(expextedPorts));
    }

    /**
     * <pre>
     *    _____  ____
     * *--+-*-+  |  |
     *    |   |\/|  |
     * *--+-*-+/\|  |
     *    |---|  |--|
     * </pre>
     */
    @Test
    public void givenGraphTransportsSwitchedPortsThroughHierarchyFromEastToWest_ShouldContainNoCrossing() {
        LNode[] leftNodes = addNodesToLayer(2, makeLayer(getGraph()));
        LNode middleOuterNode = addNodeToLayer(makeLayer(getGraph()));
        LNode rightNode = addNodeToLayer(makeLayer(getGraph()));
         setFixedOrderConstraint(rightNode);

        LPort[] leftPortsRightNode = addPortsOnSide(2, rightNode, PortSide.WEST);
        LPort[] middleOuterRightPorts = addPortsOnSide(2, middleOuterNode, PortSide.EAST);
        addEdgeBetweenPorts(middleOuterRightPorts[0], leftPortsRightNode[0]);
        addEdgeBetweenPorts(middleOuterRightPorts[1], leftPortsRightNode[1]);

        LPort[] middleOuterLeftPorts = addPortsOnSide(2, middleOuterNode, PortSide.WEST);
        eastWestEdgeFromTo(leftNodes[0], middleOuterLeftPorts[1]);
        eastWestEdgeFromTo(leftNodes[1], middleOuterLeftPorts[0]);

        LGraph innerGraph = nestedGraph(middleOuterNode);
        LNode[] leftInnerDummyNodes =
                addExternalPortDummiesToLayer(makeLayer(innerGraph), middleOuterLeftPorts);
        LNode[] innerNodes = addNodesToLayer(2, makeLayer(innerGraph));
        LNode[] rightInnerDummyNodes =
                addExternalPortDummiesToLayer(makeLayer(innerGraph), middleOuterRightPorts);

        eastWestEdgeFromTo(innerNodes[0], rightInnerDummyNodes[0]);
        eastWestEdgeFromTo(innerNodes[1], rightInnerDummyNodes[1]);
        eastWestEdgeFromTo(leftInnerDummyNodes[0], innerNodes[0]);
        eastWestEdgeFromTo(leftInnerDummyNodes[1], innerNodes[1]);

        List<LPort> expectedPortOrderMiddleNode =
                copyPortsInIndexOrder(middleOuterNode, 1, 0, 3, 2);

        // sweep backward
        getRandom().setNextBoolean(false);
        setUpAndMinimizeCrossings();

        assertThat(getGraph().getLayers().get(1).getNodes().get(0).getPorts(),
                is(expectedPortOrderMiddleNode));
    }

    /**
     * <pre>
     *  *-------
     * ____  __|_
     * |*-+  +-*|
     * |  |\/|  |
     * |*-+/\+-*|
     * |--|  |--|
     *   |_____*
     *
     * </pre>
     */
    @Test
    public void givenSimpleHierarchicalCrossWithNorthSouthEdge_ShouldResultInNoCrossing() {
        LNode[] leftNodes = addNodesToLayer(3, makeLayer(getGraph()));
        LNode leftOuterNode = leftNodes[1];
        LNode[] rightNodes = addNodesToLayer(3, makeLayer(getGraph()));
        LNode rightOuterNode = rightNodes[1];

        // adds first port, so must be before other edges.
        addNorthSouthEdge(PortSide.NORTH, rightOuterNode, rightNodes[0], leftNodes[0], true);
        LPort[] leftOuterPorts = addPortsOnSide(2, leftOuterNode, PortSide.EAST);
        LPort[] rightOuterPorts = addPortsOnSide(2, rightOuterNode, PortSide.WEST);
        addEdgeBetweenPorts(leftOuterPorts[0], rightOuterPorts[0]);
        addEdgeBetweenPorts(leftOuterPorts[1], rightOuterPorts[1]);

        makeNestedTwoNodeGraphWithEasternPorts(leftOuterNode, leftOuterPorts);

        addNorthSouthEdge(PortSide.SOUTH, leftOuterNode, leftNodes[2], rightNodes[2], false);

        LGraph rightInnerGraph =
                makeNestedTwoNodeGraphWithWesternPorts(rightOuterNode, rightOuterPorts);

        List<LNode> expectedExternalDummyOrderRight =
                switchOrderOfNodesInLayer(0, 1, rightInnerGraph.getLayers().get(0));
        List<LNode> expectedNormalNodeOrderRight =
                switchOrderOfNodesInLayer(0, 1, rightInnerGraph.getLayers().get(1));

        List<LPort> expectedOrderOfPortsRight = Lists.newArrayList(rightOuterNode.getPorts().get(0),
                rightOuterPorts[1], rightOuterPorts[0]);

        setUpAndMinimizeCrossings();
        List<LNode> actualExternalDummyOrderRight = rightInnerGraph.getLayers().get(0).getNodes();
        assertThat(actualExternalDummyOrderRight, is(expectedExternalDummyOrderRight));

        assertThat(rightOuterNode.getPorts(), is(expectedOrderOfPortsRight));

        List<LNode> actualNormalOrderRight = rightInnerGraph.getLayers().get(1).getNodes();
        assertThat(actualNormalOrderRight, is(expectedNormalNodeOrderRight));
    }

    /**
     * <pre>
     * ____  _____
     * |*-+  +-* |
     * |  |\/|   |
     * |*-+/\+-* |
     * |  |  |*-*|
     * |  |  |*-*|
     * |  |  |*-*|
     * |  |  |*-*|
     * |  |  |*-*|
     * |  |  |*-*|
     * |  |  |*-*|
     * |  |  |*-*|
     * |--|  |---|
     * </pre>
     */
    @Test
    public void givenSimpleHierarchicalNodeWithLowConnectivity_DoesNotUseHierarchySweep() {
        if (crossMinType == CrossMinType.BARYCENTER) {

            LNode leftOuterNode = addNodeToLayer(makeLayer(getGraph()));
            LNode rightOuterNode = addNodeToLayer(makeLayer(getGraph()));
            LPort[] leftOuterPorts = addPortsOnSide(2, leftOuterNode, PortSide.EAST);
            LPort[] rightOuterPorts = addPortsOnSide(2, rightOuterNode, PortSide.WEST);

            addEdgeBetweenPorts(leftOuterPorts[0], rightOuterPorts[0]);
            addEdgeBetweenPorts(leftOuterPorts[1], rightOuterPorts[1]);
            LGraph leftInnerGraph =
                    makeNestedTwoNodeGraphWithEasternPorts(leftOuterNode, leftOuterPorts);
            LGraph rightInnerGraph =
                    makeNestedTwoNodeGraphWithWesternPorts(rightOuterNode, rightOuterPorts);
            Layer leftLayer = rightInnerGraph.getLayers().get(1);
            LNode[] leftNodes = addNodesToLayer(7, leftLayer);
            LNode[] rightNodes = addNodesToLayer(7, makeLayer(rightInnerGraph));
            for (int i = 0; i < 7; i++) {
                eastWestEdgeFromTo(leftNodes[i], rightNodes[i]);
            }

            getGraph().setProperty(LayeredOptions.CROSSING_MINIMIZATION_HIERARCHICAL_SWEEPINESS, 0.2);

            List<LNode> expectedExternalDummyOrderRight =
                    Lists.newArrayList(rightInnerGraph.getLayers().get(0));
            List<LPort> expectedOrderOfPortsRight = Lists.newArrayList(rightOuterPorts);
            List<LNode> expectedNodeOrderLeftInner =
                    switchOrderOfNodesInLayer(0, 1, leftInnerGraph.getLayers().get(0));

            setUpAndMinimizeCrossings();
            List<LNode> actualExternalDummyOrderRight =
                    rightInnerGraph.getLayers().get(0).getNodes();
            assertThat(actualExternalDummyOrderRight, is(expectedExternalDummyOrderRight));

            assertThat(rightOuterNode.getPorts(), is(expectedOrderOfPortsRight));

            List<LNode> actualNormalOrderLeft = leftInnerGraph.getLayers().get(0).getNodes();
            assertThat(actualNormalOrderLeft, is(expectedNodeOrderLeftInner));
        }
    }

    /**
     * <pre>
     *
     * __________________
     * |___        ___  |
     * || |===*  --| |  |    ______
     * || |    \/  | |  |    |*  *|
     * || |    /\  | |  |----| \/ |
     * ||_|===*  --|_|  |    | /\ |
     * |      ___  ___  |    |*  *|
     * |      | |\/| |  |    |____|
     * |      |_|/\|_|  |
     * |____________^___|
     * .            |
     *           port orders fixed
     * </pre>
     *
     * First Layer and last layer in fixed order.
     *
     * @return graph of the form above.
     */
    @Test
    public void inRecursiveLayout_althoughBackwardSweepNotTakenStillCorrectsPortOrder() {
        LNode leftOuterNode = addNodeToLayer(makeLayer());
        LNode rightOuterNode = addNodeToLayer(makeLayer());

        eastWestEdgeFromTo(leftOuterNode, rightOuterNode);

        LGraph leftInnerGraph = nestedGraph(leftOuterNode);
        Layer[] layers = makeLayers(3, leftInnerGraph);
        LNode leftNode = addNodeToLayer(layers[0]);
        LNode[] middleNodes = addNodesToLayer(3, layers[1]);
        LNode[] rightNodes = addNodesToLayer(2, layers[2]);

        eastWestEdgeFromTo(middleNodes[0], rightNodes[0]);
        eastWestEdgeFromTo(middleNodes[1], rightNodes[0]);

        eastWestEdgeFromTo(middleNodes[2], rightNodes[1]);
        eastWestEdgeFromTo(middleNodes[2], rightNodes[1]);
        eastWestEdgeFromTo(leftNode, middleNodes[0]);
        eastWestEdgeFromTo(leftNode, middleNodes[0]);
        eastWestEdgeFromTo(leftNode, middleNodes[1]);
        eastWestEdgeFromTo(leftNode, middleNodes[1]);

        setFixedOrderConstraint(rightNodes[1]);
        setFixedOrderConstraint(rightNodes[0]);
        setFixedOrderConstraint(leftNode);

        LGraph rightInnerGraph = nestedGraph(rightOuterNode);
        LNode[] rightGraphLeftNodes = addNodesToLayer(2, makeLayer(rightInnerGraph));
        LNode[] rightGraphRightNodes = addNodesToLayer(2, makeLayer(rightInnerGraph));

        eastWestEdgeFromTo(rightGraphLeftNodes[0], rightGraphRightNodes[1]);
        eastWestEdgeFromTo(rightGraphLeftNodes[1], rightGraphRightNodes[0]);
        setUpIds();

        if (crossMinType == CrossMinType.BARYCENTER) {

            List<LPort> expectedPortOrderBottomMiddleNode =
                    copyPortsInIndexOrder(middleNodes[2], 1, 0);

            setUpAndMinimizeCrossings();

            Layer middleLayer = leftInnerGraph.getLayers().get(1);
            assertTrue(middleLayer.toString(),
                    inOrder(middleNodes[0], middleNodes[1], middleLayer));
            assertThat(middleNodes[2].getPorts(), is(expectedPortOrderBottomMiddleNode));
        }
    }

    /**
     * <pre>
     * _______
     * |*  *-|--*
     * | \/  |
     * | /\  |
     * |*  *-|--*
     * |*--* |
     * |*--* |
     * |*--* |
     * |*--* |
     * |*--* |
     * |*--* |
     * |-----|
     * </pre>
     */
    @Test
    public void recursiveLayout_givenCrossInFirstLevelCompoundNode_sortsPortsAccordingly() {
        // parent graph
        LNode leftOuterNode = addNodeToLayer(makeLayer(getGraph()));
        LNode[] rightNodes = addNodesToLayer(2, makeLayer(getGraph()));
        LPort[] leftOuterPorts = addPortsOnSide(2, leftOuterNode, PortSide.EAST);
        eastWestEdgeFromTo(leftOuterPorts[0], rightNodes[0]);
        eastWestEdgeFromTo(leftOuterPorts[1], rightNodes[1]);

        // child graph
        LGraph leftInnerGraph = nestedGraph(leftOuterNode);
        LNode[] leftInnerNodesLeft = addNodesToLayer(9, makeLayer(leftInnerGraph));
        LNode[] leftInnerNodesRight = addNodesToLayer(9, makeLayer(leftInnerGraph));
        LNode[] leftInnerDummyNodes =
                addExternalPortDummiesToLayer(makeLayer(leftInnerGraph), leftOuterPorts);
        eastWestEdgeFromTo(leftInnerNodesRight[0], leftInnerDummyNodes[0]);
        eastWestEdgeFromTo(leftInnerNodesRight[1], leftInnerDummyNodes[1]);
        eastWestEdgeFromTo(leftInnerNodesLeft[0], leftInnerNodesRight[1]);
        eastWestEdgeFromTo(leftInnerNodesLeft[1], leftInnerNodesRight[0]);

        for (int i = 2; i < leftInnerNodesLeft.length; i++) {
            eastWestEdgeFromTo(leftInnerNodesLeft[i], leftInnerNodesRight[i]);
        }

        List<LNode> expectedOrderRight = copyOfSwitchOrderOfNodesInLayer(0, 1, 1);
        List<LPort> expectedPortOrderLeft = copyPortsInIndexOrder(leftOuterNode, 1, 0);

        setUpAndMinimizeCrossings();

        if (crossMinType == CrossMinType.BARYCENTER) {
            assertThat(leftOuterNode.getPorts(), is(expectedPortOrderLeft));
            assertThat(getGraph().getLayers().get(1).getNodes(), is(expectedOrderRight));
        }
    }

    /**
     * <pre>
     * ___
     * | |==|
     * | |  |
     * | |--+--*
     * |_|--|--*
     *  ^
     *fixed port order
     * </pre>
     */
    @Test
    public void doesNotChangeLayoutDueToHyperedgeCounter() {
        LNode leftNode = addNodeToLayer(makeLayer());
        LNode[] rightNodes = addNodesToLayer(2, makeLayer());
        LPort lowerRightNodePort = addPortOnSide(rightNodes[1], PortSide.WEST);
        LPort lowerRightNodePort2 = addPortOnSide(rightNodes[1], PortSide.WEST);
        LPort[] leftNodePorts = addPortsOnSide(3, leftNode, PortSide.EAST);

        addEdgeBetweenPorts(leftNodePorts[0], lowerRightNodePort2);
        addEdgeBetweenPorts(leftNodePorts[0], lowerRightNodePort2);
        eastWestEdgeFromTo(leftNodePorts[1], rightNodes[0]);
        addEdgeBetweenPorts(leftNodePorts[2], lowerRightNodePort);
        setPortOrderFixed(leftNode);

        List<LNode> expectedOrderRight = Lists.newArrayList(rightNodes);

        getRandom().setNextBoolean(false);
        setUpAndMinimizeCrossings();

        assertThat(getGraph().getLayers().get(1).getNodes(), is(expectedOrderRight));
    }

    /**
     * <pre>
     *   --*
     *   |    ____
     *   //*--|-*|
     * //|    |  |
     *|| --*--|-*|
     *||      |__|
     * \\
     *  \\
     *   \\
     *     *
     * </pre>
     */
    @Test
    public void inLayerCrossingsInLeftMostLayer() {
        LNode[] nodes = addNodesToLayer(4, makeLayer());
        LNode rightNode = addNodeToLayer(makeLayer());
        LPort[] rightOuterPorts = addPortsOnSide(2, rightNode, PortSide.WEST);
        makeNestedTwoNodeGraphWithWesternPorts(rightNode, rightOuterPorts);
        eastWestEdgeFromTo(nodes[1], rightOuterPorts[1]);
        eastWestEdgeFromTo(nodes[2], rightOuterPorts[0]);
        addInLayerEdge(nodes[0], nodes[2], PortSide.WEST);
        final PortSide portSide = PortSide.WEST;
        LPort portOne = addPortOnSide(nodes[1], portSide);
        LPort portTwo = addPortOnSide(nodes[3], portSide);
        addEdgeBetweenPorts(portOne, portTwo);
        addEdgeBetweenPorts(portOne, portTwo);
        setFixedOrderConstraint(rightNode);
        List<LNode> actualOrder = getGraph().getLayers().get(0).getNodes();

        List<LNode> expectedOrder = getListCopyInIndexOrder(actualOrder, 1, 0, 2, 3);

        setUpAndMinimizeCrossings();
        assertThat(actualOrder, is(expectedOrder));
    }

    /**
     * <pre>
     *       *
     *      || 
     *   ----*
     *   |  ____
     *   /-|-*|
     * / | |  |
     *|  --|-*| *
     *|    |__|
     * \
     *  \
     *   \
     *     --*
     *      ||
     *       *
     * </pre>
     */
    @Test
    public void inLayerCrossHierarchyCrossingsInLeftMostLayer() {
        LNode[] nodes = addNodesToLayer(5, makeLayer());
        addPortsOnSide(2, nodes[0], PortSide.WEST);
        addNodeToLayer(makeLayer());
        LPort[] topNodePorts = addPortsOnSide(2, nodes[0], PortSide.WEST);
        LPort[] secondNodePorts = addPortsOnSide(3, nodes[1], PortSide.WEST);
        addEdgeBetweenPorts(topNodePorts[0], secondNodePorts[2]);
        addEdgeBetweenPorts(topNodePorts[1], secondNodePorts[1]);

        LPort[] outerPorts = addPortsOnSide(2, nodes[2], PortSide.WEST);
        LGraph inner = makeNestedTwoNodeGraphWithWesternPorts(nodes[2], outerPorts);

        LPort[] fifthNodePorts = addPortsOnSide(3, nodes[3], PortSide.WEST);
        LPort[] sixthNodePorts = addPortsOnSide(2, nodes[4], PortSide.WEST);
        addEdgeBetweenPorts(fifthNodePorts[0], sixthNodePorts[1]);
        addEdgeBetweenPorts(fifthNodePorts[1], sixthNodePorts[0]);

        addEdgeBetweenPorts(outerPorts[0], fifthNodePorts[2]);

        List<LPort> actualPortOrder = nodes[2].getPorts();
        List<LPort> expectedPortOrder = getListCopyInIndexOrder(actualPortOrder, 1, 0);
        List<LNode> actualNodeOrder = graph.getLayers().get(0).getNodes();
        List<LNode> expectedNodeOrder = Lists.newArrayList(actualNodeOrder);
        List<LNode> actualInnerNodeOrder = inner.getLayers().get(1).getNodes();
        List<LNode> expectedInnerNodeOrder = getListCopyInIndexOrder(actualInnerNodeOrder, 1, 0);

        graph.setProperty(LayeredOptions.CROSSING_MINIMIZATION_HIERARCHICAL_SWEEPINESS, 1.0);
        setUpAndMinimizeCrossings();
        assertThat(actualNodeOrder, is(expectedNodeOrder));
        assertThat(actualPortOrder, is(expectedPortOrder));
        assertThat(actualInnerNodeOrder, is(expectedInnerNodeOrder));

    }

    /**
     * <pre>
     * ______
     * |*  *+-*
     * | \/ |
     * | /\ |
     * |*  *+-*
     * |----|
     * </pre>
     */
    @Test
    public void allRecursive_crossInFirstLevelCompoundNode_ShouldRemoveCrossing() {
        // parent graph
        LNode leftOuterNode = addNodeToLayer(makeLayer(getGraph()));
        LNode[] rightNodes = addNodesToLayer(2, makeLayer(getGraph()));
        LPort[] leftOuterPorts = addPortsOnSide(2, leftOuterNode, PortSide.EAST);
        eastWestEdgeFromTo(leftOuterPorts[0], rightNodes[0]);
        eastWestEdgeFromTo(leftOuterPorts[1], rightNodes[1]);

        // child graph
        LGraph leftInnerGraph = nestedGraph(leftOuterNode);
        LNode[] leftInnerNodesleft = addNodesToLayer(2, makeLayer(leftInnerGraph));
        LNode[] leftInnerNodesRight = addNodesToLayer(2, makeLayer(leftInnerGraph));
        LNode[] leftInnerDummyNodes =
                addExternalPortDummiesToLayer(makeLayer(leftInnerGraph), leftOuterPorts);
        eastWestEdgeFromTo(leftInnerNodesRight[0], leftInnerDummyNodes[0]);
        eastWestEdgeFromTo(leftInnerNodesRight[1], leftInnerDummyNodes[1]);
        eastWestEdgeFromTo(leftInnerNodesleft[0], leftInnerNodesRight[1]);
        eastWestEdgeFromTo(leftInnerNodesleft[1], leftInnerNodesRight[0]);
        getGraph().setProperty(LayeredOptions.CROSSING_MINIMIZATION_HIERARCHICAL_SWEEPINESS, -1.0);

        List<LNode> expectedOrderRight = copyOfSwitchOrderOfNodesInLayer(0, 1, 1);
        List<LNode> expectedInnerOrderRight = Lists.newArrayList(leftInnerNodesRight[1], leftInnerNodesRight[0]);

        setUpAndMinimizeCrossings();

        assertThat(leftInnerGraph.getLayers().get(1).getNodes(), is(expectedInnerOrderRight));
        assertThat(getGraph().getLayers().get(1).getNodes(), is(expectedOrderRight));
    }

    /**
     * <pre>
     * _______  _____
     * |   *-+--+-* |
     * |   *-+--+-* |
     * |_____|  |___|
     * </pre>
     */
    @Test
    public void needsHierarchical_checkWhetherMarked() {
        LNode leftOuterNode = addNodeToLayer(makeLayer());
        LNode rightOuterNode = addNodeToLayer(makeLayer());
        LPort[] leftOuterPorts = addPortsOnSide(2, leftOuterNode, PortSide.EAST);
        LPort[] rightOuterPorts = addPortsOnSide(2, rightOuterNode, PortSide.WEST);
        addEdgeBetweenPorts(leftOuterPorts[0], rightOuterPorts[1]);
        addEdgeBetweenPorts(leftOuterPorts[1], rightOuterPorts[0]);
        makeNestedTwoNodeGraphWithEasternPorts(leftOuterNode, leftOuterPorts);
        makeNestedTwoNodeGraphWithWesternPorts(rightOuterNode, rightOuterPorts);

        graph.setProperty(LayeredOptions.CROSSING_MINIMIZATION_HIERARCHICAL_SWEEPINESS, 0.1);
        
        setUpAndMinimizeCrossings();
        
        List<GraphInfoHolder> graphData = crossMin.getGraphData();
        for (GraphInfoHolder data : graphData) {
            if (data.hasParent()) {
                assertFalse(data.dontSweepInto());
            }
        }
    }
    /**
     * <pre>
     * _______   _____
     * |   *-+-*-+-* |
     * |   *-+-*-+-* |
     * |_____|   |___|
     *         ^
     *   fixed port order
     * </pre>
     */
    @Test
    public void fixedOrderMustTransfer_checkWhetherMarked() {
        LNode leftOuterNode = addNodeToLayer(makeLayer());
        LNode[] middleNodes = addNodesToLayer(2, makeLayer());
        setFixedOrderConstraint(middleNodes[0]);
        setFixedOrderConstraint(middleNodes[1]);
        LNode rightOuterNode = addNodeToLayer(makeLayer());
        LPort[] leftOuterPorts = addPortsOnSide(2, leftOuterNode, PortSide.EAST);
        LPort[] rightOuterPorts = addPortsOnSide(2, rightOuterNode, PortSide.WEST);
        eastWestEdgeFromTo(leftOuterPorts[0], middleNodes[0]);
        eastWestEdgeFromTo(leftOuterPorts[1], middleNodes[1]);
        eastWestEdgeFromTo(middleNodes[0], rightOuterPorts[0]);
        eastWestEdgeFromTo(middleNodes[1], rightOuterPorts[1]);
        makeNestedTwoNodeGraphWithEasternPorts(leftOuterNode, leftOuterPorts);
        makeNestedTwoNodeGraphWithWesternPorts(rightOuterNode, rightOuterPorts);

        setUpAndMinimizeCrossings();

        List<GraphInfoHolder> graphData = crossMin.getGraphData();
        for (GraphInfoHolder data : graphData) {
            if (data.hasParent()) {
                assertFalse(data.dontSweepInto());
            }
        }
    }

    /**
     * <pre>
     * _______ ___ _____
     * |   *-+-|*|-+-* |
     * |     | |_| |   |
     * |     | ___ |   |
     * |   *-+-|*|-+-* |
     * |     | |_| |   |
     * |_____|     |___|
     * </pre>
     */
    @Test
    public void hierarchicalMustTransfer_checkWhetherMarked() {
        LNode leftOuterNode = addNodeToLayer(makeLayer());
        LNode[] middleNodes = addNodesToLayer(2, makeLayer());
        addNodeToLayer(makeLayer(nestedGraph(middleNodes[0])));
        addNodeToLayer(makeLayer(nestedGraph(middleNodes[1])));

        LNode rightOuterNode = addNodeToLayer(makeLayer());
        LPort[] leftOuterPorts = addPortsOnSide(2, leftOuterNode, PortSide.EAST);
        LPort[] rightOuterPorts = addPortsOnSide(2, rightOuterNode, PortSide.WEST);
        eastWestEdgeFromTo(leftOuterPorts[0], middleNodes[0]);
        eastWestEdgeFromTo(leftOuterPorts[1], middleNodes[1]);
        eastWestEdgeFromTo(middleNodes[0], rightOuterPorts[0]);
        eastWestEdgeFromTo(middleNodes[1], rightOuterPorts[1]);
        makeNestedTwoNodeGraphWithEasternPorts(leftOuterNode, leftOuterPorts);
        makeNestedTwoNodeGraphWithWesternPorts(rightOuterNode, rightOuterPorts);
        
        setUpAndMinimizeCrossings();
        
        setOnAllGraphs(LayeredOptions.CROSSING_MINIMIZATION_HIERARCHICAL_SWEEPINESS, -1.0, graph);
        List<GraphInfoHolder> graphData = crossMin.getGraphData();
        assertFalse(graphData.get(leftOuterNode.getNestedGraph().id).dontSweepInto());
        assertFalse(
                graphData.get(rightOuterNode.getNestedGraph().id).dontSweepInto());
    }

    
    /**
     * <pre>
     * _______  _____
     * |   *-+--+   |
     * |   *-+--+   |
     * |_____|  |___|
     *            ^
     *      fixed port order
     * </pre>
     */
    @Test
    public void needsHierarchicalBecauseOfFixedOrder_checkWhetherMarked() {
        LNode leftOuterNode = addNodeToLayer(makeLayer());
        LNode rightOuterNode = addNodeToLayer(makeLayer());
        setFixedOrderConstraint(rightOuterNode);
        LPort[] leftOuterPorts = addPortsOnSide(2, leftOuterNode, PortSide.EAST);
        LPort[] rightOuterPorts = addPortsOnSide(2, rightOuterNode, PortSide.WEST);
        addEdgeBetweenPorts(leftOuterPorts[0], rightOuterPorts[1]);
        addEdgeBetweenPorts(leftOuterPorts[1], rightOuterPorts[0]);
        makeNestedTwoNodeGraphWithEasternPorts(leftOuterNode, leftOuterPorts);
        
        setUpAndMinimizeCrossings();
        
        List<GraphInfoHolder> graphData = crossMin.getGraphData();
        for (GraphInfoHolder data : graphData) {
            if (data.hasParent()) {
                assertFalse(data.dontSweepInto());
            }
        }
    }
    /**
     * <pre>
     * _______  _____
     * |     +--+-* |
     * |     +--+-* |
     * |_____|  |___|
     *    ^
     * fixed port order
     * </pre>
     */
    @Test
    public void needsHierarchicalBecauseOfFixedOrder2_checkWhetherMarked() {
        LNode leftOuterNode = addNodeToLayer(makeLayer());
        LNode rightOuterNode = addNodeToLayer(makeLayer());
        setFixedOrderConstraint(leftOuterNode);
        LPort[] leftOuterPorts = addPortsOnSide(2, leftOuterNode, PortSide.EAST);
        LPort[] rightOuterPorts = addPortsOnSide(2, rightOuterNode, PortSide.WEST);
        addEdgeBetweenPorts(leftOuterPorts[0], rightOuterPorts[1]);
        addEdgeBetweenPorts(leftOuterPorts[1], rightOuterPorts[0]);
        makeNestedTwoNodeGraphWithWesternPorts(rightOuterNode, rightOuterPorts);
        
        setUpAndMinimizeCrossings();
        
        List<GraphInfoHolder> graphData = crossMin.getGraphData();
        for (GraphInfoHolder data : graphData) {
            if (data.hasParent()) {
                assertFalse(data.dontSweepInto());
            }
        }
    }


    /**
     * <pre>
     * _______  _____
     * |   * +--+ * |
     * |   * +--+ * |
     * |_____|  |___|
     * </pre>
     */
    @Ignore // TODO case where bottom up has sortable ports.
    public void needsCrossHierarchy_checkWhetherMarked() {
        LNode leftOuterNode = addNodeToLayer(makeLayer());
        LNode rightOuterNode = addNodeToLayer(makeLayer());
        setFixedOrderConstraint(rightOuterNode);
        LPort[] leftOuterPorts = addPortsOnSide(2, leftOuterNode, PortSide.EAST);
        LPort[] rightOuterPorts = addPortsOnSide(2, rightOuterNode, PortSide.WEST);
        addEdgeBetweenPorts(leftOuterPorts[0], rightOuterPorts[1]);
        addEdgeBetweenPorts(leftOuterPorts[1], rightOuterPorts[0]);

        LGraph leftNestedGraph = nestedGraph(leftOuterNode);
        addNodesToLayer(2, makeLayer(leftNestedGraph));
        addExternalPortDummiesToLayer(makeLayer(leftNestedGraph), leftOuterPorts);
        LGraph rightNestedGraph = nestedGraph(rightOuterNode);
        addNodesToLayer(2, makeLayer(rightNestedGraph));
        addExternalPortDummiesToLayer(makeLayer(rightNestedGraph), rightOuterPorts);

        setUpAndMinimizeCrossings();
        List<GraphInfoHolder> graphData = crossMin.getGraphData();
        for (GraphInfoHolder data : graphData) {
            if (data.hasParent()) {
                assertFalse(data.dontSweepInto());
            }
        }
    }

    /**
     * <pre>
     *    _____
     * *--+-* |
     *    p-* |
     *    |___|
     * </pre>
     * 
     * p is unconnected port.
     */
    @Test
    public void singleCrossHierarchEdge_needsNoCrossHierarchy() {
        LNode leftNode = addNodeToLayer(makeLayer());
        LNode rightOuterNode = addNodeToLayer(makeLayer());
        LPort[] rightOuterPorts = addPortsOnSide(2, rightOuterNode, PortSide.WEST);
        eastWestEdgeFromTo(leftNode, rightOuterPorts[1]);
        setOnAllGraphs(LayeredOptions.CROSSING_MINIMIZATION_HIERARCHICAL_SWEEPINESS, 0.1, graph);
        makeNestedTwoNodeGraphWithWesternPorts(rightOuterNode, rightOuterPorts);
        
        setUpAndMinimizeCrossings();
        List<GraphInfoHolder> graphData = crossMin.getGraphData();
        for (GraphInfoHolder data : graphData) {
            if (data.hasParent()) {
                assertFalse(data.dontSweepInto());
            }
        }
    }

    /**
     * <pre>
     *_____   _______
     *| *-+---p *\  |
     *|   |   | *-* |
     *|   |   |  x  |
     *|   |   | * | |
     *| *-+---+-*-*-+-*
     *|   |   |     |‚
     *|___|   |_____|
     * </pre>
     * 
     * p is unconnected port.
     */
    @Test
    public void moreRandomThanHierarchicalPaths_ShouldNotCrossHierarchy() {
        LNode leftOuterNode = addNodeToLayer(makeLayer());
        LPort[] leftOuterPorts = addPortsOnSide(2, leftOuterNode, PortSide.EAST);
        makeNestedTwoNodeGraphWithEasternPorts(leftOuterNode, leftOuterPorts);
        LNode middleOuterNode = addNodeToLayer(makeLayer());
        LNode rightNode = addNodeToLayer(makeLayer());
        LPort[] rightOuterPorts = addPortsOnSide(1, middleOuterNode, PortSide.EAST);
        LPort[] middleNodeLeftOuterPorts = addPortsOnSide(2, middleOuterNode, PortSide.WEST);
        addEdgeBetweenPorts(leftOuterPorts[0], middleNodeLeftOuterPorts[1]);
        addEdgeBetweenPorts(leftOuterPorts[1], middleNodeLeftOuterPorts[0]);
        eastWestEdgeFromTo(rightOuterPorts[0], rightNode);

        LGraph nestedGraph = nestedGraph(middleOuterNode);
        LNode[] leftDummies = addExternalPortDummiesToLayer(makeLayer(nestedGraph), middleNodeLeftOuterPorts);
        LNode[] leftInnerNodes = addNodesToLayer(4, makeLayer(nestedGraph));
        LNode[] rightInnerNodes = addNodesToLayer(2, makeLayer(nestedGraph));
        LNode[] rightDummies = addExternalPortDummiesToLayer(makeLayer(nestedGraph), rightOuterPorts);
        eastWestEdgeFromTo(leftDummies[0], leftInnerNodes[3]);
        eastWestEdgeFromTo(leftInnerNodes[0], rightInnerNodes[0]);
        eastWestEdgeFromTo(leftInnerNodes[1], rightInnerNodes[0]);
        eastWestEdgeFromTo(leftInnerNodes[1], rightInnerNodes[1]);
        eastWestEdgeFromTo(leftInnerNodes[2], rightInnerNodes[0]);
        eastWestEdgeFromTo(leftInnerNodes[3], rightInnerNodes[1]);
        eastWestEdgeFromTo(rightInnerNodes[1], rightDummies[0]);

        setOnAllGraphs(LayeredOptions.CROSSING_MINIMIZATION_HIERARCHICAL_SWEEPINESS, 0.0, nestedGraph);
        setUpAndMinimizeCrossings();
        List<GraphInfoHolder> graphData = crossMin.getGraphData();
        if (crossMinType == CrossMinType.BARYCENTER) {
            assertTrue(graphData.get(nestedGraph.id).dontSweepInto());
        }
    }

    /**
     * <pre>
     *_____   _______
     *| *-+---p *-| |
     *|   |   | *-* |
     *|   |   |  x  |
     *|   |   | * | |
     *| *-+---+-*-*-+-*
     *|   |   |     |‚
     *|___|   |_____|
     * </pre>
     * 
     * p is unconnected port.
     */
    @Test
    public void moreRandomThanHierarchicalWithNorthSout_ShouldNotCrossHierarchy() {
        LNode leftOuterNode = addNodeToLayer(makeLayer());
        LPort[] leftOuterPorts = addPortsOnSide(2, leftOuterNode, PortSide.EAST);
        makeNestedTwoNodeGraphWithEasternPorts(leftOuterNode, leftOuterPorts);
        LNode middleOuterNode = addNodeToLayer(makeLayer());
        LNode rightNode = addNodeToLayer(makeLayer());
        LPort[] rightOuterPorts = addPortsOnSide(1, middleOuterNode, PortSide.EAST);
        LPort[] middleNodeLeftOuterPorts = addPortsOnSide(2, middleOuterNode, PortSide.WEST);
        addEdgeBetweenPorts(leftOuterPorts[0], middleNodeLeftOuterPorts[1]);
        addEdgeBetweenPorts(leftOuterPorts[1], middleNodeLeftOuterPorts[0]);
        eastWestEdgeFromTo(rightOuterPorts[0], rightNode);

        LGraph nestedGraph = nestedGraph(middleOuterNode);
        LNode[] leftDummies = addExternalPortDummiesToLayer(makeLayer(nestedGraph), middleNodeLeftOuterPorts);
        LNode[] leftInnerNodes = addNodesToLayer(4, makeLayer(nestedGraph));
        LNode[] rightInnerNodes = addNodesToLayer(3, makeLayer(nestedGraph));
        LNode[] rightDummies = addExternalPortDummiesToLayer(makeLayer(nestedGraph), rightOuterPorts);
        eastWestEdgeFromTo(leftDummies[0], leftInnerNodes[3]);
        addNorthSouthEdge(PortSide.NORTH, rightInnerNodes[1], rightInnerNodes[0], leftInnerNodes[0], true);
        eastWestEdgeFromTo(leftInnerNodes[1], rightInnerNodes[1]);
        eastWestEdgeFromTo(leftInnerNodes[1], rightInnerNodes[2]);
        eastWestEdgeFromTo(leftInnerNodes[2], rightInnerNodes[1]);
        eastWestEdgeFromTo(leftInnerNodes[3], rightInnerNodes[2]);
        eastWestEdgeFromTo(rightInnerNodes[1], rightDummies[0]);

        setUpAndMinimizeCrossings();
        List<GraphInfoHolder> graphData = crossMin.getGraphData();
        if (crossMinType == CrossMinType.BARYCENTER) {
            assertTrue(graphData.get(nestedGraph.id).dontSweepInto());
        }
    }

    /**
     * <pre>
     *       -*
     *____  / *
     *|*-0-x /
     *|  |/ x
     *|*-1-/ \
     *|__|    *
     * ^
     * port 0 before 1
     * </pre>
     * 
     * @throws Exception
     */
    // TODO This problem can currently not be solved by our algorithm :-(
    @Ignore
    public void sortOneWrongPort() throws Exception {
        LNode leftNode = addNodeToLayer(makeLayer());
        LPort[] ports = addPortsOnSide(3, leftNode, PortSide.EAST);
        LGraph leftInnerGraph = nestedGraph(leftNode);
        LNode[] leftInnerNodes = addNodesToLayer(2, makeLayer(leftInnerGraph));
        LNode[] dummies = addExternalPortDummiesToLayer(makeLayer(leftInnerGraph), ports);
        eastWestEdgeFromTo(leftInnerNodes[0], dummies[0]);
        eastWestEdgeFromTo(leftInnerNodes[1], dummies[2]);


        LNode[] rightNodes = addNodesToLayer(3, makeLayer());
        eastWestEdgeFromTo(ports[0], rightNodes[2]);
        eastWestEdgeFromTo(ports[1], rightNodes[0]);
        eastWestEdgeFromTo(ports[2], rightNodes[1]);

        List<LPort> expectedPortOrder =
                Lists.newArrayList(leftNode.getPorts().get(1), leftNode.getPorts().get(0), leftNode.getPorts().get(2));

        random.setNextBoolean(false);
        getGraph().setProperty(LayeredOptions.CROSSING_MINIMIZATION_HIERARCHICAL_SWEEPINESS, -1.0);
        new LayerSweepCrossingMinimizer(CrossMinType.BARYCENTER).process(graph, new BasicProgressMonitor());

        assertEquals(leftNode.getPorts(), expectedPortOrder);
    }

    private LPort[] reverse(final LPort[] rightOuterPorts) {
        LPort[] res = new LPort[rightOuterPorts.length];
        for (int i = 0; i < rightOuterPorts.length; i++) {
            LPort lPort = rightOuterPorts[i];
            res[rightOuterPorts.length - 1 - i] = lPort;
        }
        return res;
    }

    private void setUpAndMinimizeCrossings() {
        if (crossMinType == CrossMinType.BARYCENTER) {
            getGraph().setProperty(LayeredOptions.THOROUGHNESS, 1);
        }
        crossMin.process(getGraph(), new BasicProgressMonitor());
    }
}
