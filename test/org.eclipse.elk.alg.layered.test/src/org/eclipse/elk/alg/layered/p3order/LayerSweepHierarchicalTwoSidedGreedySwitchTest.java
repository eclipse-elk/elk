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

import java.util.List;

import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.intermediate.greedyswitch.TestGraphCreator;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.layered.p3order.LayerSweepCrossingMinimizer.CrossMinType;
import org.eclipse.elk.core.data.LayoutMetaDataService;
import org.eclipse.elk.core.options.PortSide;
import org.eclipse.elk.core.util.BasicProgressMonitor;
import org.junit.BeforeClass;
import org.junit.Test;

import com.google.common.collect.Lists;

// CHECKSTYLEOFF javadoc
// CHECKSTYLEOFF MagicNumber
// CHECKSTYLEOFF MethodName
@SuppressWarnings("restriction")
public class LayerSweepHierarchicalTwoSidedGreedySwitchTest extends TestGraphCreator {
    
    @BeforeClass
    public static void initialize() {
        LayoutMetaDataService.getInstance();
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
    public void givencrossInFirstLevelCompoundNodeUsingTwoSidedGreedySwitch_LeavesCrossing() {
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

        List<LNode> expectedSameOrder = copyOfNodesInLayer(1);

        setUpAndMinimizeCrossings();

        assertThat(getGraph().getLayers().get(1).getNodes(), is(expectedSameOrder));
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

        List<LNode> expectedOrderLeft =
                switchOrderOfNodesInLayer(0, 1, innerGraph.getLayers().get(0));

        setUpAndMinimizeCrossings();

        assertThat(innerGraph.getLayers().get(0).getNodes(), is(expectedOrderLeft));
    }

    /**
     * <pre>
     * _________
     * | ___   |  ___
     * | | |--*p _| |
     * | | |   |x | |
     * | |_|--*p -|_|
     * |po fix | po fix
     * |_______|
     * </pre>
     * 
     * @throws Exception
     */
    @Test
    public void onlyTwoSidedGreedySwitchReturnsNoChange() throws Exception {
        LNode leftNode = addNodeToLayer(makeLayer(graph));
        LNode rightNode = addNodeToLayer(makeLayer(graph));
        LPort[] leftPorts = addPortsOnSide(2, leftNode, PortSide.EAST);
        eastWestEdgeFromTo(leftPorts[0], rightNode);
        eastWestEdgeFromTo(leftPorts[1], rightNode);
        setFixedOrderConstraint(rightNode);
        LGraph innerGraph = nestedGraph(leftNode);
        Layer[] innerLayers = makeLayers(2, innerGraph);
        LNode innerLeftNode = addNodeToLayer(innerLayers[0]);
        setFixedOrderConstraint(innerLeftNode);
        LNode[] dummies = addExternalPortDummiesToLayer(innerLayers[1], leftPorts);
        eastWestEdgeFromTo(innerLeftNode, dummies[0]);
        eastWestEdgeFromTo(innerLeftNode, dummies[1]);
        List<LNode> expectedOrderDummies = Lists.newArrayList(innerGraph.getLayers().get(1));
        // sweep backward
        random.setNextBoolean(false);
        setUpAndMinimizeCrossings();

        assertThat(innerGraph.getLayers().get(1).getNodes(), is(expectedOrderDummies));
    }

    /**
     * <pre>
     * _________
     * | ___   |  ___
     * | | |\ *p__| |
     * | | | x |  | |
     * | |_|/ *p--|_|
     * |po fix | po fix
     * |_______|
     * </pre>
     * 
     * @throws Exception
     */
    @Test
    public void crossingBeforeLastLayerCausesCrossingOutside_TwoSidedPrevents() throws Exception {
        LNode leftNode = addNodeToLayer(makeLayer(graph));
        LNode rightNode = addNodeToLayer(makeLayer(graph));
        LPort[] leftPorts = addPortsOnSide(2, leftNode, PortSide.EAST);
        eastWestEdgeFromTo(leftPorts[1], rightNode);
        eastWestEdgeFromTo(leftPorts[0], rightNode);
        setFixedOrderConstraint(rightNode);
        LGraph innerGraph = nestedGraph(leftNode);
        Layer[] innerLayers = makeLayers(2, innerGraph);
        LNode innerLeftNode = addNodeToLayer(innerLayers[0]);
        setFixedOrderConstraint(innerLeftNode);
        LNode[] dummies = addExternalPortDummiesToLayer(innerLayers[1], leftPorts);
        eastWestEdgeFromTo(innerLeftNode, dummies[1]);
        eastWestEdgeFromTo(innerLeftNode, dummies[0]);
        List<LNode> expectedOrderDummies = Lists.newArrayList(innerGraph.getLayers().get(1));

        setUpAndMinimizeCrossings();

        assertThat(innerGraph.getLayers().get(1).getNodes(), is(expectedOrderDummies));
    }

    /**
     * <pre>
     *        _________
     * ___   |   ___ |
     * | |---p* /| | |
     * | |   | x | | |
     * |_|---p* \|_| |
     * po fix| po fix|
     *       |_______|
     * </pre>
     * 
     * @throws Exception
     */
    @Test
    public void crossingBeforeFirstLayerCausesCrossingOutside_TwoSidedPrevents() throws Exception {
        LNode leftNode = addNodeToLayer(makeLayer(graph));
        LNode rightNode = addNodeToLayer(makeLayer(graph));
        LPort[] rightPorts = addPortsOnSide(2, rightNode, PortSide.WEST);
        eastWestEdgeFromTo(leftNode, rightPorts[1]);
        eastWestEdgeFromTo(leftNode, rightPorts[0]);
        setFixedOrderConstraint(leftNode);
        LGraph innerGraph = nestedGraph(rightNode);
        Layer[] innerLayers = makeLayers(2, innerGraph);
        LNode innerRightNode = addNodeToLayer(innerLayers[1]);
        setFixedOrderConstraint(innerRightNode);
        LNode[] dummies = addExternalPortDummiesToLayer(innerLayers[0], rightPorts);
        eastWestEdgeFromTo(dummies[0], innerRightNode);
        eastWestEdgeFromTo(dummies[1], innerRightNode);
        List<LNode> expectedOrderDummies = Lists.newArrayList(innerGraph.getLayers().get(0));
        // backwardSweep
        random.setNextBoolean(false);
        setUpAndMinimizeCrossings();

        assertThat(innerGraph.getLayers().get(0).getNodes(), is(expectedOrderDummies));
    }

    private void setUpAndMinimizeCrossings() {
        LayerSweepCrossingMinimizer crossMin =
                new LayerSweepCrossingMinimizer(CrossMinType.TWO_SIDED_GREEDY_SWITCH);
        crossMin.process(getGraph(), new BasicProgressMonitor());
    }
}
