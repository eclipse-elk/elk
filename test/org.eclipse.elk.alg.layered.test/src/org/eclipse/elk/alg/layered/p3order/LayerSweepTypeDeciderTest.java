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

import static org.junit.Assert.assertTrue;

import java.util.Arrays;

import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.intermediate.greedyswitch.TestGraphCreator;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.layered.p3order.LayerSweepCrossingMinimizer.CrossMinType;
import org.eclipse.elk.alg.layered.p3order.counting.IInitializable;
import org.eclipse.elk.core.data.LayoutMetaDataService;
import org.eclipse.elk.core.options.PortSide;
import org.junit.BeforeClass;
import org.junit.Test;

@SuppressWarnings("restriction")
public class LayerSweepTypeDeciderTest extends TestGraphCreator {
    
    @BeforeClass
    public static void initialize() {
        LayoutMetaDataService.getInstance();
    }

    /**
     * <pre>
     *   __________
     *   | d--*   |
     *   | |  |   |
     *   | *  d-* |
     * *=+-*--*-*-+=*
     *   |________|
     * </pre>
     */
    @Test
    public void northSouthPortsChangeDecision() {
        Layer[] layers = makeLayers(3);
        LNode leftOuter = addNodeToLayer(layers[0]);
        LNode middleOuter = addNodeToLayer(layers[1]);
        LNode rightOuter = addNodeToLayer(layers[2]);
        LPort leftPortMiddleNode = addPortOnSide(middleOuter, PortSide.WEST);
        LPort rightPortMiddleNode = addPortOnSide(middleOuter, PortSide.EAST);
        eastWestEdgeFromTo(middleOuter, rightOuter);
        eastWestEdgeFromTo(rightPortMiddleNode, rightOuter);
        eastWestEdgeFromTo(leftOuter, leftPortMiddleNode);
        eastWestEdgeFromTo(leftOuter, middleOuter);

        LGraph nestedGraph = nestedGraph(middleOuter);
        Layer[] innerLayers = makeLayers(5, nestedGraph);
        LNode leftDummy = addExternalPortDummyNodeToLayer(innerLayers[0], leftPortMiddleNode);
        LNode[] firstLayer = addNodesToLayer(3, innerLayers[1]);
        LNode[] secondLayer = addNodesToLayer(3, innerLayers[2]);
        LNode[] thirdLayer = addNodesToLayer(2, innerLayers[3]);
        LNode rightDummy = addExternalPortDummyNodeToLayer(innerLayers[4], rightPortMiddleNode);
        eastWestEdgeFromTo(leftDummy, firstLayer[2]);
        eastWestEdgeFromTo(firstLayer[2], secondLayer[2]);
        eastWestEdgeFromTo(secondLayer[2], thirdLayer[1]);
        eastWestEdgeFromTo(thirdLayer[1], rightDummy);

        addNorthSouthEdge(PortSide.NORTH, firstLayer[1], firstLayer[0], secondLayer[0], false);
        addNorthSouthEdge(PortSide.SOUTH, secondLayer[0], secondLayer[1], thirdLayer[0], false);

        graph.id = 0;
        nestedGraph.id = 1;
        setOnAllGraphs(LayeredOptions.CROSSING_MINIMIZATION_HIERARCHICAL_SWEEPINESS, 0.1, graph);
        GraphInfoHolder gd = new GraphInfoHolder(nestedGraph, CrossMinType.BARYCENTER,
                Arrays.asList(new GraphInfoHolder(graph, CrossMinType.BARYCENTER, null)));
        assertTrue(gd.dontSweepInto());
    }

    /**
     * <pre>
     *   ____  
     * *=|=*=|=*
     *   |___|
     * 
     * 
     * </pre>
     */
    @Test
    public void allHierarchicalButMinusOne_StillReturnsBU() throws Exception {
        LNode leftNode = addNodeToLayer(makeLayer());
        LNode middleNode = addNodeToLayer(makeLayer());
        LNode rightNode = addNodeToLayer(makeLayer());
        LPort[] middlePortRight = addPortsOnSide(2, middleNode, PortSide.EAST);
        LPort[] middlePortLeft = addPortsOnSide(2, middleNode, PortSide.WEST);
        eastWestEdgeFromTo(leftNode, middlePortLeft[1]);
        eastWestEdgeFromTo(leftNode, middlePortLeft[0]);
        eastWestEdgeFromTo(middlePortRight[0], rightNode);
        eastWestEdgeFromTo(middlePortRight[1], rightNode);

        LGraph innerGraph = nestedGraph(middleNode);
        Layer[] layers = makeLayers(3, innerGraph);
        LNode[] leftInnerDummyNodes =
                addExternalPortDummiesToLayer(layers[2], middlePortRight);
        LNode node = addNodeToLayer(layers[1]);
        LNode[] rightInnerDummyNodes =
                addExternalPortDummiesToLayer(layers[0], middlePortLeft);
        eastWestEdgeFromTo(leftInnerDummyNodes[0], node);
        eastWestEdgeFromTo(leftInnerDummyNodes[1], node);
        eastWestEdgeFromTo(node, rightInnerDummyNodes[0]);
        eastWestEdgeFromTo(node, rightInnerDummyNodes[1]);

        graph.id = 0;
        innerGraph.id = 1;

        setOnAllGraphs(LayeredOptions.CROSSING_MINIMIZATION_HIERARCHICAL_SWEEPINESS, -1.0, graph);
        GraphInfoHolder gd = new GraphInfoHolder(innerGraph, CrossMinType.BARYCENTER,
                Arrays.asList(new GraphInfoHolder(graph, CrossMinType.BARYCENTER, null)));
        LayerSweepTypeDecider td = new LayerSweepTypeDecider(gd);
        IInitializable.init(Arrays.asList(td), innerGraph.toNodeArray());
        assertTrue(td.useBottomUp());
    }

    /**
     * <pre>
     * _____  
     * |p*=|=*
     * |___|
     * 
     * 
     * </pre>
     */
    @Test
    public void portWithNoEdges_DoesNotCount() throws Exception {
        LNode middleNode = addNodeToLayer(makeLayer());
        LNode rightNode = addNodeToLayer(makeLayer());
        LPort[] middlePortRight = addPortsOnSide(2, middleNode, PortSide.EAST);
        eastWestEdgeFromTo(middlePortRight[0], rightNode);
        eastWestEdgeFromTo(middlePortRight[1], rightNode);

        LGraph innerGraph = nestedGraph(middleNode);
        Layer[] layers = makeLayers(2, innerGraph);
        LNode node = addNodeToLayer(layers[0]);
        LNode[] rightInnerDummyNodes = addExternalPortDummiesToLayer(layers[1], middlePortRight);
        eastWestEdgeFromTo(node, rightInnerDummyNodes[0]);
        eastWestEdgeFromTo(node, rightInnerDummyNodes[1]);
        addPortOnSide(node, PortSide.WEST);

        graph.id = 0;
        innerGraph.id = 1;

        setOnAllGraphs(LayeredOptions.CROSSING_MINIMIZATION_HIERARCHICAL_SWEEPINESS, -0.1, graph);
        GraphInfoHolder gd = new GraphInfoHolder(innerGraph, CrossMinType.BARYCENTER,
                Arrays.asList(new GraphInfoHolder(graph, CrossMinType.BARYCENTER, null)));
        assertTrue(gd.dontSweepInto());
    }


}
