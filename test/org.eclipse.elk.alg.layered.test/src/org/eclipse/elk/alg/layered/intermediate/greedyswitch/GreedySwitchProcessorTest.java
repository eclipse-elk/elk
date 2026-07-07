/*******************************************************************************
 * Copyright (c) 2016 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.layered.intermediate.greedyswitch;

import static org.hamcrest.CoreMatchers.is;
import static org.junit.Assert.assertThat;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.p3order.LayerSweepCrossingMinimizer;
import org.eclipse.elk.alg.layered.p3order.LayerSweepCrossingMinimizer.CrossMinType;
import org.eclipse.elk.alg.test.PlainJavaInitialization;
import org.eclipse.elk.core.options.PortSide;
import org.eclipse.elk.core.util.BasicProgressMonitor;
import org.eclipse.elk.core.util.IElkProgressMonitor;
import org.junit.BeforeClass;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameters;

/**
 * Tests all variants of greedy switch. All three one-sided methods must behave
 * the same and all three two sided methods must behave the same.
 *
 * @author alan
 */
@RunWith(Parameterized.class)
public class GreedySwitchProcessorTest extends TestGraphCreator {
    private final LayerSweepCrossingMinimizer greedySwitcher;
    private final IElkProgressMonitor monitor;
    private final CrossMinType greedyType;

    /**
     * Constructor called by Parameterized.
     *
     * @param gT
     *            greedyType
     */
    public GreedySwitchProcessorTest(final CrossMinType gT) {
        greedyType = gT;
        greedySwitcher = new LayerSweepCrossingMinimizer(gT);

        monitor = new BasicProgressMonitor();
    }
    
    @BeforeClass
    public static void initPlainJavaLayout() {
        PlainJavaInitialization.initializePlainJavaLayout();
    }

    /**
     * Sets the Parameters to be tested as all elements of the CrossMinType enum.
     *
     * @return parameters
     */
    @Parameters(name = "{0}")
    public static Iterable<Object[]> greedyTypes() {
        return Arrays.asList(
                new Object[][] { { CrossMinType.ONE_SIDED_GREEDY_SWITCH }, { CrossMinType.TWO_SIDED_GREEDY_SWITCH } });
    }

    // CHECKSTYLEOFF javadoc
    // CHECKSTYLEOFF MagicNumber
    @Test
    public void shouldSwitchCross() {
        getCrossFormedGraph();

        List<LNode> expectedOrderLayerOne;
        List<LNode> expectedOrderLayerTwo;
        if (greedyType == CrossMinType.ONE_SIDED_GREEDY_SWITCH) {
            expectedOrderLayerOne = copyOfNodesInLayer(0);
            expectedOrderLayerTwo = copyOfSwitchOrderOfNodesInLayer(0, 1, 1);
        } else {
            expectedOrderLayerOne = copyOfSwitchOrderOfNodesInLayer(0, 1, 0);
            expectedOrderLayerTwo = copyOfNodesInLayer(1);
        }

        startGreedySwitcherWithCurrentType();

        assertThat("Layer one", copyOfNodesInLayer(0), is(expectedOrderLayerOne));
        assertThat("Layer two", copyOfNodesInLayer(1), is(expectedOrderLayerTwo));
    }

    private void startGreedySwitcherWithCurrentType() {
        greedySwitcher.process(getGraph(), monitor);
    }

    @Test
    public void constraintsPreventSwitchInSecondLayer() {
        getCrossFormedGraphWithConstraintsInSecondLayer();

        List<LNode> expectedOrderLayerOne = copyOfSwitchOrderOfNodesInLayer(0, 1, 0);
        List<LNode> expectedOrderLayerTwo = copyOfNodesInLayer(1);

        startGreedySwitcherWithCurrentType();

        assertThat("Layer one", copyOfNodesInLayer(0), is(expectedOrderLayerOne));
        assertThat("Layer two", copyOfNodesInLayer(1), is(expectedOrderLayerTwo));
    }

    @Test
    public void constraintsPreventAnySwitch() {
        getCrossFormedGraphConstraintsPreventAnySwitch();

        List<LNode> expectedOrderLayerOne = copyOfNodesInLayer(0);
        List<LNode> expectedOrderLayerTwo = copyOfNodesInLayer(1);

        startGreedySwitcherWithCurrentType();

        assertThat("Layer one", copyOfNodesInLayer(0), is(expectedOrderLayerOne));
        assertThat("Layer two", copyOfNodesInLayer(1), is(expectedOrderLayerTwo));
    }

    @Test
    public void layoutUnitConstraintPreventsSwitch() {
        getNodesInDifferentLayoutUnitsPreventSwitch();

        List<LNode> expectedOrderLayerTwo = copyOfNodesInLayer(1);

        startGreedySwitcherWithCurrentType();

        assertThat("Layer one", copyOfNodesInLayer(1), is(expectedOrderLayerTwo));
    }

    @Test
    public void oneNode() {
        getOneNodeGraph();
        int layerIndex = 0;
        copyOfSwitchOrderOfNodesInLayer(0, 0, layerIndex);
        // should cause no errors
    }

    @Test
    public void inLayerSwitchable() {
        getInLayerEdgesGraph();

        List<LNode> expectedOrder = copyOfSwitchOrderOfNodesInLayer(0, 1, 1);

        startGreedySwitcherWithCurrentType();

        assertThat(copyOfNodesInLayer(1), is(expectedOrder));
    }

    @Test
    public void multipleEdgesBetweenSameNodes() {
        getMultipleEdgesBetweenSameNodesGraph();

        List<LNode> expectedOrderLayerOne;
        List<LNode> expectedOrderLayerTwo;
        if (greedyType == CrossMinType.ONE_SIDED_GREEDY_SWITCH) {
            expectedOrderLayerOne = copyOfNodesInLayer(0);
            expectedOrderLayerTwo = copyOfSwitchOrderOfNodesInLayer(0, 1, 1);
        } else {
            expectedOrderLayerOne = copyOfSwitchOrderOfNodesInLayer(0, 1, 0);
            expectedOrderLayerTwo = copyOfNodesInLayer(1);
        }

        startGreedySwitcherWithCurrentType();

        assertThat("Layer one", copyOfNodesInLayer(0), is(expectedOrderLayerOne));
        assertThat("Layer two", copyOfNodesInLayer(1), is(expectedOrderLayerTwo));
    }

    @Test
    public void selfLoops() {
        Layer leftLayer = makeLayer(getGraph());
        Layer rightLayer = makeLayer(getGraph());

        LNode topLeft = addNodeToLayer(leftLayer);
        LNode bottomLeft = addNodeToLayer(leftLayer);
        LNode topRight = addNodeToLayer(rightLayer);
        LNode bottomRight = addNodeToLayer(rightLayer);

        LPort topLeftPort = addPortOnSide(topLeft, PortSide.EAST);
        LPort bottomLeftPort = addPortOnSide(bottomLeft, PortSide.EAST);
        setUpIds();
        LGraph selfLoopCrossGraph = getGraph();
        for (Layer layer : selfLoopCrossGraph) {
            for (LNode node : layer) {
                selfLoopOn(node, PortSide.EAST);
                selfLoopOn(node, PortSide.EAST);
                selfLoopOn(node, PortSide.EAST);
                selfLoopOn(node, PortSide.WEST);
                selfLoopOn(node, PortSide.WEST);
                selfLoopOn(node, PortSide.WEST);
            }
        }
        LPort topRightPort = addPortOnSide(topRight, PortSide.WEST);
        LPort bottomRightPort = addPortOnSide(bottomRight, PortSide.WEST);

        addEdgeBetweenPorts(topLeftPort, bottomRightPort);
        addEdgeBetweenPorts(bottomLeftPort, topRightPort);

        List<LNode> expectedOrderLayerOne;
        List<LNode> expectedOrderLayerTwo;
        if (greedyType == CrossMinType.ONE_SIDED_GREEDY_SWITCH) {
            expectedOrderLayerOne = copyOfNodesInLayer(0);
            expectedOrderLayerTwo = copyOfSwitchOrderOfNodesInLayer(0, 1, 1);
        } else {
            expectedOrderLayerOne = copyOfSwitchOrderOfNodesInLayer(0, 1, 0);
            expectedOrderLayerTwo = copyOfNodesInLayer(1);
        }

        startGreedySwitcherWithCurrentType();

        assertThat("Layer one", copyOfNodesInLayer(0), is(expectedOrderLayerOne));
        assertThat("Layer two", copyOfNodesInLayer(1), is(expectedOrderLayerTwo));
    }

    /**
     * <pre>
     *     ______
     *     |____|
     *      |  |
     *      *--+--*
     *         |
     *         *--*
     * </pre>
     *
     * The one-side decider should switch north-south port crossings in the
     * center layer.
     */
    @Test
    public void northSouthPortCrossing() {
        graph = new NorthSouthEdgeTestGraphCreator().getNorthSouthDownwardCrossingGraph();

        int layerIndex = 0;
        List<LNode> expectedOrderTwoSided = new ArrayList<LNode>(copyOfNodesInLayer(layerIndex));
        List<LNode> expectedOrderOneSided = copyOfSwitchOrderOfNodesInLayer(1, 2, layerIndex);

        startGreedySwitcherWithCurrentType();

        if (greedyType == CrossMinType.ONE_SIDED_GREEDY_SWITCH) {
            assertThat(copyOfNodesInLayer(layerIndex), is(expectedOrderOneSided));
        } else {
            assertThat(copyOfNodesInLayer(layerIndex), is(expectedOrderTwoSided));
        }
    }

    /**
     * <pre>
     * *\  --*
     *   \/ /
     * *-*===*
     *  x /
     * * * --*
     * Port order not fixed.
     * </pre>
     */
    @Test
    public void moreComplex() {
        getMoreComplexThreeLayerGraph();
        List<LNode> expectedOrderLayerTwo;
        List<LNode> expectedOrderLayerThree;
        expectedOrderLayerTwo = copyOfNodesInLayer(1);
        expectedOrderLayerThree = copyOfSwitchOrderOfNodesInLayer(0, 1, 2);
        startGreedySwitcherWithCurrentType();
        assertThat("Layer two", copyOfNodesInLayer(1), is(expectedOrderLayerTwo));
        assertThat("Layer three", copyOfNodesInLayer(2), is(expectedOrderLayerThree));
    }

    @Test
    public void switchOnlyForOneSided() {
        getSwitchOnlyOneSided();

        int layerIndex = 1;
        List<LNode> expectedOrderOneSided = copyOfSwitchOrderOfNodesInLayer(0, 1, layerIndex);
        List<LNode> expectedOrderTwoSided = new ArrayList<LNode>(copyOfNodesInLayer(layerIndex));

        startGreedySwitcherWithCurrentType();
        if (greedyType == CrossMinType.ONE_SIDED_GREEDY_SWITCH) {
            assertThat(copyOfNodesInLayer(layerIndex), is(expectedOrderOneSided));
        } else {
            assertThat(copyOfNodesInLayer(layerIndex), is(expectedOrderTwoSided));
        }
    }

    @Test
    public void doesNotWorsenCrossAmount() {
        getGraphWhichCouldBeWorsenedBySwitch();
        List<LNode> expectedOrderFirstLayer = new ArrayList<LNode>(copyOfNodesInLayer(0));
        List<LNode> expectedOrderSecondLayer = new ArrayList<LNode>(copyOfNodesInLayer(1));

        startGreedySwitcherWithCurrentType();

        assertThat("Layer one", copyOfNodesInLayer(0), is(expectedOrderFirstLayer));
        assertThat("Layer two", copyOfNodesInLayer(1), is(expectedOrderSecondLayer));
    }

    @Test
    public void switchMoreThanOnce() {
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

        List<LNode> oneSidedFirstLayer = new ArrayList<LNode>(copyOfNodesInLayer(0));
        List<LNode> oneSidedFirstSwitchSecondLayer = copyOfSwitchOrderOfNodesInLayer(0, 1, 1);
        List<LNode> oneSidedsecondSwitchSecondLayer = getCopyWithSwitchedOrder(2, 3, oneSidedFirstSwitchSecondLayer);
        List<LNode> oneSidedThirdSwitchSecondLayer = getCopyWithSwitchedOrder(1, 2, oneSidedsecondSwitchSecondLayer);

        List<LNode> twoSidedFirstLayer = copyOfSwitchOrderOfNodesInLayer(0, 1, 0);
        List<LNode> twoSidedFirstSwitchSecondLayer = copyOfSwitchOrderOfNodesInLayer(1, 2, 1);
        List<LNode> twoSidedsecondSwitchSecondLayer = getCopyWithSwitchedOrder(0, 1, twoSidedFirstSwitchSecondLayer);
        startGreedySwitcherWithCurrentType();
        if (greedyType == CrossMinType.ONE_SIDED_GREEDY_SWITCH) {
            assertThat("Layer one" + copyOfNodesInLayer(0), copyOfNodesInLayer(0), is(oneSidedFirstLayer));
            assertThat("Layer two " + copyOfNodesInLayer(1), copyOfNodesInLayer(1), is(oneSidedThirdSwitchSecondLayer));
        } else {
            assertThat("Layer one " + copyOfNodesInLayer(0), copyOfNodesInLayer(0), is(twoSidedFirstLayer));
            assertThat("Layer two " + copyOfNodesInLayer(1), copyOfNodesInLayer(1), is(twoSidedsecondSwitchSecondLayer));
        }

    }

}
