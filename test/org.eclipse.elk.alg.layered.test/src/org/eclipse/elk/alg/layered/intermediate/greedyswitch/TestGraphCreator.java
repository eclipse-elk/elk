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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Deque;
import java.util.List;
import java.util.Random;
import java.util.Set;

import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LNode.NodeType;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.options.GraphProperties;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.core.options.EdgeRouting;
import org.eclipse.elk.core.options.HierarchyHandling;
import org.eclipse.elk.core.options.PortConstraints;
import org.eclipse.elk.core.options.PortSide;
import org.eclipse.elk.graph.properties.IProperty;
import org.eclipse.elk.graph.properties.MapPropertyHolder;

import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;

/**
 * Use to create test graphs. CAUTION: Layout algorithm assumes the ports to be
 * ordered in a clockwise manner. You must think about this yourself when
 * constructing a test graph. This means that the methods for creating edges
 * cannot be used in every case. TODO consider moving all downward into base
 *
 * @author alan
 */
public class TestGraphCreator {
    private int portId = 0;
    private int nodeId = 0;
    /** the graph. */
    protected LGraph graph;

    /**
     * @return the graph
     */
    public LGraph getGraph() {
        setUpGraph(graph);
        return graph;
    }

    private int edgeId = 0;
    protected MockRandom random;

    /**
     * Makes a fancy test graph creator.
     */
    public TestGraphCreator() {
        graph = new LGraph();
        random = new MockRandom();
        setUpGraph(graph);
    }

    protected LGraph setUpGraph(final LGraph g) {
        setUpIds();
        g.setProperty(LayeredOptions.EDGE_ROUTING, EdgeRouting.ORTHOGONAL);
        g.setProperty(InternalProperties.RANDOM, random);
        g.setProperty(LayeredOptions.HIERARCHY_HANDLING, HierarchyHandling.INCLUDE_CHILDREN);
        return g;
    }

    /**
     * Creates empty graph.
     *
     * @return return empty graph.
     */
    public LGraph getEmptyGraph() {
        graph = new LGraph();
        setUpGraph(graph);
        return graph;
    }

    /**
     * Creates two nodes with no connection between them.
     *
     * @return graph with two nodes with no connection between them.
     */
    public LGraph getTwoNodesNoConnectionGraph() {
        Layer layer = makeLayer(graph);
        addNodeToLayer(layer);
        addNodeToLayer(layer);
        setUpIds();
        return graph;
    }

    /**
     * Set layer id globally, node id per layer and port id globally.
     */
    protected void setUpIds() {
        Deque<LGraph> graphs = Lists.newLinkedList();
        graphs.add(graph);
        while (!graphs.isEmpty()) {
            int lId = 0;
            int pId = 0;
            LGraph g = graphs.pop();
            for (Layer l : g) {
                l.id = lId++;
                int i = 0;
                for (LNode n : l) {
                    LGraph nestedGraph = n.getNestedGraph();
                    if (nestedGraph != null) {
                        graphs.add(nestedGraph);
                    }
                    n.id = i++;
                    for (LPort p : n.getPorts()) {
                        p.id = pId++;
                    }
                }
            }
        }
    }

    // CHECKSTYLEOFF MagicNumber
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
    public LGraph getCrossFormedGraph() {

        Layer leftLayer = makeLayer(graph);
        Layer rightLayer = makeLayer(graph);

        LNode topLeft = addNodeToLayer(leftLayer);
        LNode bottomLeft = addNodeToLayer(leftLayer);
        LNode topRight = addNodeToLayer(rightLayer);
        LNode bottomRight = addNodeToLayer(rightLayer);

        eastWestEdgeFromTo(topLeft, bottomRight);
        eastWestEdgeFromTo(bottomLeft, topRight);
        setUpIds();
        return graph;
    }

    /**
     * <pre>
     * *
     *  \\
     *   \\
     * *---*
     * .
     * </pre>
     *
     * @return Graph of the form above.
     */
    public LGraph multipleEdgesAndSingleEdge() {

        Layer leftLayer = makeLayer(graph);
        Layer rightLayer = makeLayer(graph);

        LNode topLeft = addNodeToLayer(leftLayer);
        LNode bottomLeft = addNodeToLayer(leftLayer);
        LNode bottomRight = addNodeToLayer(rightLayer);

        eastWestEdgeFromTo(topLeft, bottomRight);
        eastWestEdgeFromTo(topLeft, bottomRight);
        eastWestEdgeFromTo(bottomLeft, bottomRight);
        setUpIds();
        return graph;
    }

    /**
     * <pre>
     * *  *  <- this node must be ...
     *  \/
     *  /\
     * *  *  <- before this node.
     * .
     * </pre>
     *
     * @return Graph of the form above.
     */
    public LGraph getCrossFormedGraphWithConstraintsInSecondLayer() {
        getCrossFormedGraph();
        Layer layerOne = graph.getLayers().get(1);
        LNode topNode = layerOne.getNodes().get(0);
        LNode secondNode = layerOne.getNodes().get(1);
        setInLayerOrderConstraint(topNode, secondNode);
        setUpIds();
        return graph;
    }

    /**
     * <pre>
     * this node must be.. -> *  *  <- and this node must be ...
     *                         \/
     *                         /\
     *    before this node -> *  *  <- before this node.
     * .
     * </pre>
     *
     * @return Graph of the form above.
     */
    public LGraph getCrossFormedGraphConstraintsPreventAnySwitch() {
        Layer leftLayer = makeLayer(graph);
        Layer rightLayer = makeLayer(graph);

        LNode topLeft = addNodeToLayer(leftLayer);
        LNode bottomLeft = addNodeToLayer(leftLayer);
        LNode topRight = addNodeToLayer(rightLayer);
        LNode bottomRight = addNodeToLayer(rightLayer);

        eastWestEdgeFromTo(topLeft, bottomRight);
        eastWestEdgeFromTo(bottomLeft, topRight);
        setInLayerOrderConstraint(topRight, bottomRight);
        setInLayerOrderConstraint(topLeft, bottomLeft);
        setUpIds();
        return graph;
    }

    /**
     * Creates graph with only one node.
     *
     * @return graph with only one node.
     */
    public LGraph getOneNodeGraph() {
        Layer layer = makeLayer(graph);
        addNodeToLayer(layer);
        setUpIds();
        return graph;
    }

    /**
     * <pre>
     *   --*
     *   |
     * *-+-*-*
     *   |
     *   --*
     * .
     * </pre>
     *
     * @return graph of the form above
     */
    public LGraph getInLayerEdgesGraph() {
        Layer leftLayer = makeLayer(graph);
        Layer middleLayer = makeLayer(graph);
        Layer rightLayer = makeLayer(graph);

        LNode leftNode = addNodeToLayer(leftLayer);
        LNode[] middleNodes = addNodesToLayer(3, middleLayer);
        LNode rightNode = addNodeToLayer(rightLayer);

        // add east side ports first to get expected port ordering
        eastWestEdgeFromTo(middleNodes[1], rightNode);
        eastWestEdgeFromTo(leftNode, middleNodes[1]);
        addInLayerEdge(middleNodes[0], middleNodes[2], PortSide.WEST);
        setUpIds();
        return graph;
    }

    /**
     * <pre>
     *   --*
     *   |
     *   --*
     *
     *  *--*
     * .
     * </pre>
     *
     * @return graph of the form above
     */
    public LGraph getInLayerEdgesGraphWhichResultsInCrossingsWhenSwitched() {
        Layer leftLayer = makeLayer(graph);
        Layer rightLayer = makeLayer(graph);

        LNode leftNode = addNodeToLayer(leftLayer);
        LNode[] rightNodes = addNodesToLayer(3, rightLayer);

        addInLayerEdge(rightNodes[0], rightNodes[1], PortSide.WEST);
        eastWestEdgeFromTo(leftNode, rightNodes[2]);

        setUpIds();
        return graph;
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
     *
     * @return Graph of the form above.
     */
    public LGraph getMultipleEdgesBetweenSameNodesGraph() {
        Layer leftLayer = makeLayer(graph);
        Layer rightLayer = makeLayer(graph);

        LNode topLeft = addNodeToLayer(leftLayer);
        LNode bottomLeft = addNodeToLayer(leftLayer);
        LNode topRight = addNodeToLayer(rightLayer);
        LNode bottomRight = addNodeToLayer(rightLayer);

        eastWestEdgeFromTo(topLeft, bottomRight);
        eastWestEdgeFromTo(topLeft, bottomRight);
        eastWestEdgeFromTo(bottomLeft, topRight);
        eastWestEdgeFromTo(bottomLeft, topRight);

        setUpIds();
        return graph;
    }

    /**
     * <pre>
     * *   *
     *  \ /
     * *-+-*
     *  / \
     * *   *
     * .
     * </pre>
     *
     * @return graph of the form above
     */
    public LGraph getCrossWithExtraEdgeInBetweenGraph() {
        Layer leftLayer = makeLayer(graph);
        Layer rightLayer = makeLayer(graph);

        LNode[] leftNodes = addNodesToLayer(3, leftLayer);
        LNode[] rightNodes = addNodesToLayer(3, rightLayer);

        eastWestEdgeFromTo(leftNodes[0], rightNodes[2]);
        eastWestEdgeFromTo(leftNodes[1], rightNodes[1]);
        eastWestEdgeFromTo(leftNodes[2], rightNodes[0]);

        setUpIds();
        return graph;
    }

    /**
     * Cross formed graph, but each node has three extra self loop edges.
     *
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
    public LGraph getCrossWithManySelfLoopsGraph() {
        Layer leftLayer = makeLayer(graph);
        Layer rightLayer = makeLayer(graph);

        LNode topLeft = addNodeToLayer(leftLayer);
        LNode bottomLeft = addNodeToLayer(leftLayer);
        LNode topRight = addNodeToLayer(rightLayer);
        LNode bottomRight = addNodeToLayer(rightLayer);

        LPort topLeftPort = addPortOnSide(topLeft, PortSide.EAST);
        LPort bottomLeftPort = addPortOnSide(bottomLeft, PortSide.EAST);
        setUpIds();
        LGraph selfLoopCrossGraph = graph;
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
        return selfLoopCrossGraph;
    }

    /**
     * Test for self loop.
     * @param node
     * @param side
     */
    protected void selfLoopOn(final LNode node, final PortSide side) {
        addEdgeBetweenPorts(addPortOnSide(node, side), addPortOnSide(node, side));
    }

    /**
     * <pre>
     * *\  --*
     *   \/ /
     * *-*===*
     *  + /
     * * * --*
     * .
     * </pre>
     *
     * @return Graph of the form above.
     */
    public LGraph getMoreComplexThreeLayerGraph() {
        Layer leftLayer = makeLayer(graph);
        Layer middleLayer = makeLayer(graph);
        Layer rightLayer = makeLayer(graph);

        LNode[] leftNodes = addNodesToLayer(3, leftLayer);
        LNode[] middleNodes = addNodesToLayer(2, middleLayer);
        LNode[] rightNodes = addNodesToLayer(3, rightLayer);
        
        LPort leftMiddleNodePort = addPortOnSide(leftNodes[1], PortSide.EAST);
        LPort middleLowerNodePortEast = addPortOnSide(middleNodes[1], PortSide.EAST);
        LPort middleUpperNodePortEast = addPortOnSide(middleNodes[0], PortSide.EAST);
        LPort rightUpperNodePort = addPortOnSide(rightNodes[0], PortSide.WEST);
        LPort rightMiddleNodePort = addPortOnSide(rightNodes[1], PortSide.WEST);
        setUpIds();

        addEdgeBetweenPorts(middleUpperNodePortEast, rightUpperNodePort);
        addEdgeBetweenPorts(middleUpperNodePortEast, rightMiddleNodePort);
        addEdgeBetweenPorts(middleUpperNodePortEast, rightMiddleNodePort);
        eastWestEdgeFromTo(middleLowerNodePortEast, rightNodes[2]);
        eastWestEdgeFromTo(leftMiddleNodePort, middleNodes[0]);
        eastWestEdgeFromTo(middleNodes[1], rightUpperNodePort);
        eastWestEdgeFromTo(leftMiddleNodePort, middleNodes[1]);
        eastWestEdgeFromTo(leftNodes[2], middleNodes[0]);
        eastWestEdgeFromTo(leftNodes[0], middleNodes[0]);

        setUpIds();
        return graph;
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
     * Port order fixed.
     *
     * @return Graph of the form above.
     */
    public LGraph getFixedPortOrderGraph() {
        Layer leftLayer = makeLayer(graph);
        Layer rightLayer = makeLayer(graph);

        LNode leftNode = addNodeToLayer(leftLayer);
        leftNode.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_ORDER);

        LNode rightTopNode = addNodeToLayer(rightLayer);
        LNode rightBottomNode = addNodeToLayer(rightLayer);

        eastWestEdgeFromTo(leftNode, rightBottomNode);
        eastWestEdgeFromTo(leftNode, rightTopNode);

        setUpIds();
        return graph;
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
     * Port order not fixed
     *
     * @return Graph of the form above.
     */
    public LGraph getGraphNoCrossingsDueToPortOrderNotFixed() {
        Layer leftLayer = makeLayer(graph);
        Layer rightLayer = makeLayer(graph);

        LNode leftNode = addNodeToLayer(leftLayer);

        LNode rightTopNode = addNodeToLayer(rightLayer);
        LNode rightBottomNode = addNodeToLayer(rightLayer);

        eastWestEdgeFromTo(leftNode, rightBottomNode);
        eastWestEdgeFromTo(leftNode, rightTopNode);

        setUpIds();
        return graph;
    }

    /**
     * <pre>
     * *  *---*
     *  \/
     *  /\
     * *  *---*
     * .
     * </pre>
     *
     * @return Graph of the form above.
     */
    public LGraph getSwitchOnlyOneSided() {
        Layer[] layers = makeLayers(3);

        LNode[] leftNodes = addNodesToLayer(2, layers[0]);
        LNode[] middleNodes = addNodesToLayer(2, layers[1]);
        LNode[] rightNodes = addNodesToLayer(2, layers[2]);

        eastWestEdgeFromTo(middleNodes[0], rightNodes[0]);
        eastWestEdgeFromTo(middleNodes[1], rightNodes[1]);
        eastWestEdgeFromTo(leftNodes[0], middleNodes[1]);
        eastWestEdgeFromTo(leftNodes[1], middleNodes[0]);

        setUpIds();
        return graph;
    }

    /**
     * <pre>
     * *--*  *
     *     \/
     *     /\
     * *--*  *
     * .
     * </pre>
     *
     * @return Graph of the form above.
     */
    public LGraph getSwitchOnlyEastOneSided() {
        Layer[] layers = makeLayers(3);

        LNode[] leftNodes = addNodesToLayer(2, layers[0]);
        LNode[] middleNodes = addNodesToLayer(2, layers[1]);
        LNode[] rightNodes = addNodesToLayer(2, layers[2]);

        eastWestEdgeFromTo(leftNodes[0], middleNodes[0]);
        eastWestEdgeFromTo(leftNodes[1], middleNodes[1]);
        eastWestEdgeFromTo(middleNodes[0], rightNodes[1]);
        eastWestEdgeFromTo(middleNodes[1], rightNodes[0]);

        setUpIds();
        return graph;
    }

    /**
     * <pre>
     * ____
     * |  |----
     * |__|\  |
     * ____ | |
     * |  |/  |
     * |__|---|
     * .
     * </pre>
     *
     * Port order fixed.
     *
     * @return Graph of the form above.
     */
    public LGraph getFixedPortOrderInLayerEdgesDontCrossEachOther() {
        Layer layer = makeLayer(graph);
        LNode[] nodes = addNodesToLayer(2, layer);
        setFixedOrderConstraint(nodes[0]);
        setFixedOrderConstraint(nodes[1]);
        // must add ports and edges manually, due to clockwise port ordering
        LPort upperPortUpperNode = addPortOnSide(nodes[0], PortSide.EAST);
        LPort lowerPortUpperNode = addPortOnSide(nodes[0], PortSide.EAST);
        LPort upperPortLowerNode = addPortOnSide(nodes[1], PortSide.EAST);
        LPort lowerPortLowerNode = addPortOnSide(nodes[1], PortSide.EAST);
        addEdgeBetweenPorts(upperPortUpperNode, lowerPortLowerNode);
        addEdgeBetweenPorts(lowerPortUpperNode, upperPortLowerNode);
        setUpIds();
        return graph;
    }

    /**
     * <pre>
     * ____
     * |  |----
     * |__|\  |
     * ____ | |
     * |  |-+--
     * |__|-|
     * .
     * </pre>
     *
     * Port order fixed.
     *
     * @return Graph of the form above.
     */
    public LGraph getFixedPortOrderInLayerEdgesWithCrossings() {
        Layer layer = makeLayer(graph);
        LNode[] nodes = addNodesToLayer(2, layer);
        setFixedOrderConstraint(nodes[0]);
        setFixedOrderConstraint(nodes[1]);
        addInLayerEdge(nodes[0], nodes[1], PortSide.EAST);
        addInLayerEdge(nodes[0], nodes[1], PortSide.EAST);

        setUpIds();
        return graph;
    }

    /**
     * <pre>
     *      ____
     *    / |  |
     * *-+--|  |\
     *   | /|  |-+-*
     * *-++-|__| |
     *   ||      |
     *   || ___  |
     *   | \| | /
     * *-+--| |/
     * *-+--|_|
     *    \
     *     \
     *      *
     * .
     * </pre>
     *
     * Port order fixed.
     *
     * @return Graph of the form above.
     */
    public LGraph getMoreComplexInLayerGraph() {
        Layer[] layers = makeLayers(3);
        LNode[] leftNodes = addNodesToLayer(4, layers[0]);
        LNode[] middleNodes = addNodesToLayer(3, layers[1]);
        LNode rightNode = addNodeToLayer(layers[2]);
        setFixedOrderConstraint(middleNodes[0]);
        setFixedOrderConstraint(middleNodes[1]);

        eastWestEdgeFromTo(leftNodes[1], middleNodes[0]);

        eastWestEdgeFromTo(leftNodes[3], middleNodes[1]);
        eastWestEdgeFromTo(leftNodes[2], middleNodes[1]);
        addInLayerEdge(middleNodes[0], middleNodes[1], PortSide.WEST);
        eastWestEdgeFromTo(leftNodes[0], middleNodes[0]);
        addInLayerEdge(middleNodes[0], middleNodes[2], PortSide.WEST);

        addInLayerEdge(middleNodes[0], middleNodes[1], PortSide.EAST);
        eastWestEdgeFromTo(middleNodes[0], rightNode);

        setUpIds();
        return graph;

    }

    /**
     * <pre>
     * *==*  *
     *     \/
     *     /\
     * *==*  *
     * .
     * </pre>
     *
     * First Layer and last layer in fixed order.
     *
     * @return graph of the form above.
     */
    public LGraph getGraphWhichCouldBeWorsenedBySwitch() {
        Layer[] layers = makeLayers(3);
        LNode[] leftNodes = addNodesToLayer(2, layers[0]);
        LNode[] middleNodes = addNodesToLayer(2, layers[1]);
        LNode[] rightNodes = addNodesToLayer(2, layers[2]);

        setInLayerOrderConstraint(leftNodes[0], leftNodes[1]);
        setInLayerOrderConstraint(rightNodes[0], rightNodes[1]);

        eastWestEdgeFromTo(middleNodes[0], rightNodes[1]);
        eastWestEdgeFromTo(middleNodes[1], rightNodes[0]);
        eastWestEdgeFromTo(leftNodes[0], middleNodes[0]);
        eastWestEdgeFromTo(leftNodes[0], middleNodes[0]);
        eastWestEdgeFromTo(leftNodes[1], middleNodes[1]);
        eastWestEdgeFromTo(leftNodes[1], middleNodes[1]);

        setUpIds();
        return graph;
    }

    /**
     * <pre>
     *     * <-- this ...
     *    /
     * *-+-* <-- cannot switch with this
     *  / _|__
     * *  |  |
     *    |__|
     *
     * .
     * </pre>
     *
     * @return graph of the form above.
     */
    public LGraph getNodesInDifferentLayoutUnitsPreventSwitch() {
        Layer[] layers = makeLayers(2);
        LNode[] leftNodes = addNodesToLayer(2, layers[0]);
        LNode[] rightNodes = addNodesToLayer(3, layers[1]);

        eastWestEdgeFromTo(leftNodes[1], rightNodes[0]);

        addNorthSouthEdge(PortSide.EAST, rightNodes[2], rightNodes[1], leftNodes[0], true);

        rightNodes[1].setProperty(InternalProperties.IN_LAYER_LAYOUT_UNIT, rightNodes[2]);
        rightNodes[2].setProperty(InternalProperties.IN_LAYER_LAYOUT_UNIT, rightNodes[2]);

        setUpIds();
        return graph;
    }

    /**
     * <pre>
     *  ---*
     *  |
     *  | ____
     * *+-|  |
     * *+-|  |
     *   \|__|
     * Port order not fixed.
     * .
     * </pre>
     *
     * @return Graph of the form above.
     */
    public LGraph multipleInBetweenLayerEdgesIntoNodeWithNoFixedPortOrder() {
        Layer leftLayer = makeLayer(graph);
        LNode[] leftNodes = addNodesToLayer(2, leftLayer);
        Layer rightLayer = makeLayer(graph);
        LNode[] rightNodes = addNodesToLayer(2, rightLayer);

        addInLayerEdge(rightNodes[0], rightNodes[1], PortSide.WEST);
        eastWestEdgeFromTo(leftNodes[0], rightNodes[1]);
        eastWestEdgeFromTo(leftNodes[0], rightNodes[1]);
        setUpIds();
        return graph;
    }

    /**
     * <pre>
     *  ---*
     *  |
     *  | ____
     * *+-|  |
     * *+-|  |
     *  | |__|
     *   \
     *    *
     * Port order not fixed.
     * .
     * </pre>
     *
     * @return Graph of the form above.
     */
    public LGraph multipleInBetweenLayerEdgesIntoNodeWithNoFixedPortOrderCauseCrossings() {
        Layer leftLayer = makeLayer(graph);
        LNode[] leftNodes = addNodesToLayer(2, leftLayer);
        Layer rightLayer = makeLayer(graph);
        LNode[] rightNodes = addNodesToLayer(3, rightLayer);

        addInLayerEdge(rightNodes[0], rightNodes[2], PortSide.WEST);
        eastWestEdgeFromTo(leftNodes[0], rightNodes[1]);
        eastWestEdgeFromTo(leftNodes[0], rightNodes[1]);
        setUpIds();
        return graph;
    }

    /**
     * <pre>
     * *----*
     *  \\
     *   \--*
     *    --*
     * *---/
     *  \---*
     * </pre>
     *
     * .
     *
     * @return Graph of the form above.
     */
    public LGraph getSwitchedProblemGraph() {
        LNode[] leftNodes = addNodesToLayer(2, makeLayer(graph));
        LNode[] rightNodes = addNodesToLayer(4, makeLayer(graph));

        eastWestEdgeFromTo(leftNodes[1], rightNodes[2]);
        eastWestEdgeFromTo(leftNodes[1], rightNodes[3]);

        eastWestEdgeFromTo(leftNodes[0], rightNodes[0]);
        eastWestEdgeFromTo(leftNodes[0], rightNodes[1]);
        eastWestEdgeFromTo(leftNodes[0], rightNodes[2]);

        setUpIds();
        return graph;
    }

    /**
     * <pre>
     * *   *<- Into same port
     *  \//
     *  //\
     * *   *
     * </pre>
     *
     * .
     *
     * @return Graph of the form above.
     */
    public LGraph twoEdgesIntoSamePort() {
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
        LPort topRightSecondPort = addPortOnSide(topRight, PortSide.WEST);

        addEdgeBetweenPorts(bottomLeftFirstPort, topRightFirstPort);
        addEdgeBetweenPorts(bottomLeftSecondPort, topRightSecondPort);
        setUpIds();
        return graph;
    }

    /**
     * <pre>
     * *---* <- Into same port
     *   /
     *  /
     * *---*
     * </pre>
     *
     * .
     *
     * @return Graph of the form above.
     */
    public LGraph twoEdgesIntoSamePortCrossesWhenSwitched() {
        Layer leftLayer = makeLayer(graph);
        Layer rightLayer = makeLayer(graph);

        LNode topLeft = addNodeToLayer(leftLayer);
        LNode bottomLeft = addNodeToLayer(leftLayer);
        LNode topRight = addNodeToLayer(rightLayer);
        LNode bottomRight = addNodeToLayer(rightLayer);

        LPort topRightPort = addPortOnSide(topRight, PortSide.WEST);
        LPort bottomLeftPort = addPortOnSide(bottomLeft, PortSide.EAST);
        addEdgeBetweenPorts(bottomLeftPort, topRightPort);

        LPort topLeftPort = addPortOnSide(topLeft, PortSide.EAST);
        addEdgeBetweenPorts(topLeftPort, topRightPort);

        eastWestEdgeFromTo(bottomLeft, bottomRight);

        setUpIds();
        return graph;
    }

    /**
     * <pre>
     * *  *
     *  \/
     *  /\
     * *--*<- Into same port
     * </pre>
     *
     * .
     *
     * @return Graph of the form above.
     */
    public LGraph twoEdgesIntoSamePortResolvesCrossingWhenSwitched() {
        Layer leftLayer = makeLayer(graph);
        Layer rightLayer = makeLayer(graph);

        LNode topLeft = addNodeToLayer(leftLayer);
        LNode bottomLeft = addNodeToLayer(leftLayer);
        LNode topRight = addNodeToLayer(rightLayer);
        LNode bottomRight = addNodeToLayer(rightLayer);

        LPort topLeftPort = addPortOnSide(topLeft, PortSide.EAST);
        LPort bottomLeftPort = addPortOnSide(bottomLeft, PortSide.EAST);
        LPort bottomRightPort = addPortOnSide(bottomRight, PortSide.WEST);

        addEdgeBetweenPorts(topLeftPort, bottomRightPort);
        addEdgeBetweenPorts(bottomLeftPort, bottomRightPort);

        eastWestEdgeFromTo(bottomLeft, topRight);

        setUpIds();
        return graph;
    }

    /**
     * <pre>
     * *---* <- into same port
     *   /
     *  /
     * *---*
     * ^
     * |
     * Two edges into same port.
     * </pre>
     *
     * @return Graph of the form above.
     */
    public LGraph twoEdgesIntoSamePortFromEastWithFixedPortOrder() {
        Layer leftLayer = makeLayer(graph);
        Layer rightLayer = makeLayer(graph);

        LNode topLeft = addNodeToLayer(leftLayer);
        LNode bottomLeft = addNodeToLayer(leftLayer);
        LNode topRight = addNodeToLayer(rightLayer);
        LNode bottomRight = addNodeToLayer(rightLayer);

        setFixedOrderConstraint(bottomLeft);
        setFixedOrderConstraint(topRight);

        LPort topLeftPort = addPortOnSide(topLeft, PortSide.EAST);
        LPort bottomLeftPort = addPortOnSide(bottomLeft, PortSide.EAST);
        LPort topRightPort = addPortOnSide(topRight, PortSide.WEST);
        LPort bottomRightPort = addPortOnSide(bottomRight, PortSide.WEST);

        addEdgeBetweenPorts(bottomLeftPort, bottomRightPort);
        addEdgeBetweenPorts(bottomLeftPort, topRightPort);
        addEdgeBetweenPorts(topLeftPort, topRightPort);

        setUpIds();
        return graph;
    }

    /**
     * Return nodes as 2d array.
     * @param g
     * @return
     */
    public LNode[][] getCurrentOrder(final LGraph g) {
        LNode[][] nodeOrder = new LNode[g.getLayers().size()][];
        List<Layer> layers = g.getLayers();
        for (int i = 0; i < layers.size(); i++) {
            Layer layer = layers.get(i);
            List<LNode> nodes = layer.getNodes();
            nodeOrder[i] = new LNode[nodes.size()];
            for (int j = 0; j < nodes.size(); j++) {
                nodeOrder[i][j] = nodes.get(j);
            }
        }
        return nodeOrder;
    }

    // CHECKSTYLEOFF Javadoc

    public void setAsNorthSouthNode(final LNode node) {
        node.setType(NodeType.NORTH_SOUTH_PORT);
    }

    public void addNorthSouthEdge(final PortSide side, final LNode nodeWithNSPorts, final LNode northSouthDummy,
            final LNode nodeWithEastWestPorts, final boolean nodeWithEastWestPortsIsOrigin) {
        boolean normalNodeEastOfNsPortNode = nodeWithEastWestPorts.getLayer().getIndex() < nodeWithNSPorts.getLayer()
                .getIndex();
        PortSide direction = normalNodeEastOfNsPortNode ? PortSide.WEST : PortSide.EAST;

        PortSide targetNodePortSide = direction.opposed();
        LPort normalNodePort = addPortOnSide(nodeWithEastWestPorts, targetNodePortSide);

        LPort dummyNodePort = addPortOnSide(northSouthDummy, direction);

        if (nodeWithEastWestPortsIsOrigin) {
            addEdgeBetweenPorts(normalNodePort, dummyNodePort);
        } else {
            addEdgeBetweenPorts(dummyNodePort, normalNodePort);
        }

        northSouthDummy.setProperty(InternalProperties.IN_LAYER_LAYOUT_UNIT, nodeWithNSPorts);
        northSouthDummy.setProperty(InternalProperties.ORIGIN, nodeWithNSPorts);

        setAsNorthSouthNode(northSouthDummy);

        LPort originPort = addPortOnSide(nodeWithNSPorts, side);
        dummyNodePort.setProperty(InternalProperties.ORIGIN, originPort);
        originPort.setProperty(InternalProperties.PORT_DUMMY, northSouthDummy);

        List<LNode> baryAssoc = Lists.newArrayList(northSouthDummy);

        List<LNode> otherBaryAssocs = nodeWithNSPorts.getProperty(InternalProperties.BARYCENTER_ASSOCIATES);
        if (otherBaryAssocs == null) {
            nodeWithNSPorts.setProperty(InternalProperties.BARYCENTER_ASSOCIATES, baryAssoc);
        } else {
            otherBaryAssocs.addAll(baryAssoc);
        }

        if (side == PortSide.NORTH) {
            northSouthDummy.getProperty(InternalProperties.IN_LAYER_SUCCESSOR_CONSTRAINTS).add(nodeWithNSPorts);
        } else {
            nodeWithNSPorts.getProperty(InternalProperties.IN_LAYER_SUCCESSOR_CONSTRAINTS).add(northSouthDummy);
        }
    }

    public void setInLayerOrderConstraint(final LNode thisNode, final LNode beforeThisNode) {
        List<LNode> scndNodeAsList = Lists.newArrayList(beforeThisNode);
        thisNode.setProperty(InternalProperties.IN_LAYER_SUCCESSOR_CONSTRAINTS, scndNodeAsList);
    }

    public void setAsLongEdgeDummy(final LNode node) {
        node.setType(NodeType.LONG_EDGE);
        node.setProperty(InternalProperties.IN_LAYER_LAYOUT_UNIT, null);
    }

    public void setPortOrderFixed(final LNode node) {
        node.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_ORDER);
        node.getGraph().getProperty(InternalProperties.GRAPH_PROPERTIES).add(GraphProperties.NON_FREE_PORTS);
    }

    public Layer[] makeLayers(final int amount) {
        return makeLayers(amount, graph);
    }

    public Layer[] makeLayers(final int amount, final LGraph g) {
        Layer[] layers = new Layer[amount];
        for (int i = 0; i < layers.length; i++) {
            layers[i] = makeLayer(g);
        }
        return layers;
    }

    public MapPropertyHolder setFixedOrderConstraint(final LNode node) {
        return node.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_ORDER);
    }

    public void setFixedOrderConstraint(final LNode[] nodes) {
        for (LNode n : nodes) {
            setFixedOrderConstraint(n);
        }
    }

    public void addInLayerEdge(final LNode nodeOne, final LNode nodeTwo, final PortSide portSide) {
        LPort portOne = addPortOnSide(nodeOne, portSide);
        LPort portTwo = addPortOnSide(nodeTwo, portSide);
        addEdgeBetweenPorts(portOne, portTwo);
    }

    public void addInLayerEdge(final LNode nodeOne, final LPort portTwo, final PortSide portSide) {
        LPort portOne = addPortOnSide(nodeOne, portSide);
        addEdgeBetweenPorts(portOne, portTwo);
    }

    public void addInLayerEdge(final LPort portOne, final LNode nodeTwo) {
        PortSide portSide = portOne.getSide();
        LPort portTwo = addPortOnSide(nodeTwo, portSide);
        addEdgeBetweenPorts(portOne, portTwo);
    }

    public LNode[] addNodesToLayer(final int amountOfNodes, final Layer leftLayer) {
        LNode[] nodes = new LNode[amountOfNodes];
        for (int j = 0; j < amountOfNodes; j++) {
            nodes[j] = addNodeToLayer(leftLayer);
        }
        return nodes;
    }

    protected Layer makeLayer() {
        return makeLayer(graph);
    }

    public Layer makeLayer(final LGraph g) {
        List<Layer> layers = g.getLayers();
        Layer layer = new Layer(g);
        layers.add(layer);
        return layer;
    }

    public LNode addNodeToLayer(final Layer layer) {
        LNode node = new LNode(layer.getGraph());
        node.setType(NodeType.NORMAL);
        node.setProperty(InternalProperties.IN_LAYER_LAYOUT_UNIT, node);
        node.setLayer(layer);
        node.id = nodeId++;
        return node;
    }

    public void eastWestEdgeFromTo(final LNode left, final LNode right) {
        LPort leftPort = addPortOnSide(left, PortSide.EAST);
        LPort rightPort = addPortOnSide(right, PortSide.WEST);
        addEdgeBetweenPorts(leftPort, rightPort);
    }

    public void addEdgeBetweenPorts(final LPort from, final LPort to) {
        LEdge edge = new LEdge();
        edge.setSource(from);
        edge.setTarget(to);
        edge.id = edgeId++;
    }

    /**
     * Sets port constraints to fixed!
     *
     * @param node
     * @param portSide
     * @return
     */
    public LPort addPortOnSide(final LNode node, final PortSide portSide) {
        LPort port = addPortTo(node);
        port.setSide(portSide);
        if (!node.getProperty(LayeredOptions.PORT_CONSTRAINTS).isSideFixed()) {
            node.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        }
        return port;
    }

    public LPort[] addPortsOnSide(final int n, final LNode node, final PortSide portSide) {
        LPort[] ps = new LPort[n];
        for (int i = 0; i < ps.length; i++) {
            ps[i] = addPortOnSide(node, portSide);
        }
        return ps;
    }

    private LPort addPortTo(final LNode node) {
        LPort port = new LPort();
        port.setNode(node);
        port.id = portId++;
        return port;
    }

    /**
     * <pre>
     *     ____
     *     |  |
     * ----0  | <--- two edges into one port
     * | / |  |
     * | +-0  |
     * |/| |__|
     * ||\
     * ||  *
     * | \
     * ----*
     * </pre>
     *
     * Port order not fixed.
     *
     * @return Graph of the form above.
     */
    public LGraph multipleEdgesIntoOnePortAndFreePortOrder() {
        LNode[] nodes = addNodesToLayer(3, makeLayer(graph));
        addInLayerEdge(nodes[0], nodes[2], PortSide.WEST);
        LPort portUpperNode = addPortOnSide(nodes[0], PortSide.WEST);
        addEdgeBetweenPorts(portUpperNode, addPortOnSide(nodes[2], PortSide.WEST));
        addEdgeBetweenPorts(portUpperNode, addPortOnSide(nodes[1], PortSide.WEST));
        setUpIds();
        return graph;
    }

    /**
     * <pre>
     * ___
     * | |\    *
     * | |\\  /
     * |_|-++-+*
     *     || \
     *     ====*
     * </pre>
     *
     * left Node has fixed Port Order.
     *
     * @return Graph of the form above.
     */
    public LGraph getOnlyCorrectlyImprovedByBestOfForwardAndBackwardSweepsInSingleLayer() {
        LNode leftNode = addNodeToLayer(makeLayer(graph));
        setFixedOrderConstraint(leftNode);
        LNode[] rightNodes = addNodesToLayer(3, makeLayer(graph));
        eastWestEdgeFromTo(leftNode, rightNodes[2]);
        eastWestEdgeFromTo(leftNode, rightNodes[2]);
        eastWestEdgeFromTo(leftNode, rightNodes[1]);
        addInLayerEdge(rightNodes[0], rightNodes[2], PortSide.WEST);
        setUpIds();
        return graph;
    }

    protected <T> List<T> getListCopyInIndexOrder(final List<T> li, final int... is) {
        List<T> list = new ArrayList<>();
        for (int i : is) {
            list.add(li.get(i));
        }
        return list;
    }

    protected <T> T[] getArrayInIndexOrder(final T[] arr, final int... is) {
        T[] copy = Arrays.copyOf(arr, arr.length);
        int j = 0;
        for (int i : is) {
            copy[j++] = arr[i];
        }
        return copy;
    }

    protected <T> List<T> copyOfListSwitchingOrder(final int i, final int j, final List<T> list) {
        List<T> listCopy = new ArrayList<T>(list);
        T first = listCopy.get(i);
        T second = listCopy.get(j);
        listCopy.set(i, second);
        listCopy.set(j, first);
        return listCopy;
    }

    protected <T> T[] switchOrderInArray(final int i, final int j, final T[] arr) {
        T[] copy = Arrays.copyOf(arr, arr.length);
        T first = arr[i];
        T snd = arr[j];
        copy[j] = first;
        copy[i] = snd;
        return copy;
    }

    protected LNode addExternalPortDummyNodeToLayer(final Layer layer, final LPort port) {
        LNode node = addNodeToLayer(layer);
        node.setProperty(InternalProperties.ORIGIN, port);
        node.setType(NodeType.EXTERNAL_PORT);
        node.setProperty(InternalProperties.EXT_PORT_SIDE, port.getSide());
        port.setProperty(InternalProperties.PORT_DUMMY, node);
        port.setProperty(InternalProperties.INSIDE_CONNECTIONS, true);
        Set<GraphProperties> ps = node.getGraph().getProperty(InternalProperties.GRAPH_PROPERTIES);
        ps.add(GraphProperties.EXTERNAL_PORTS);
        return node;
    }

    protected LNode[] addExternalPortDummiesToLayer(final Layer layer, final LPort[] ports) {
        LNode[] nodes = new LNode[ports.length];
        PortSide side = ports[0].getSide();
        for (int i = 0; i < ports.length; i++) {
            int portIndex = side == PortSide.EAST ? i : ports.length - 1 - i;
            nodes[i] = addExternalPortDummyNodeToLayer(layer, ports[portIndex]);
        }
        return nodes;
    }

    protected LGraph nestedGraph(final LNode node) {
        node.setProperty(InternalProperties.COMPOUND_NODE, true);
        LGraph nestedGraph = node.getNestedGraph();
        if (nestedGraph == null) {
            nestedGraph = new LGraph();
            setUpGraph(nestedGraph);
            node.setNestedGraph(nestedGraph);
            nestedGraph.setParentNode(node);
        }
        return nestedGraph;
    }

    protected List<LNode> switchOrderOfNodesInLayer(final int nodeOne, final int nodeTwo, final Layer layer) {
        List<LNode> nodes = layer.getNodes();
        LNode firstNode = nodes.get(nodeOne);
        LNode secondNode = nodes.get(nodeTwo);
        List<LNode> switchedList = new ArrayList<LNode>(nodes);
        switchedList.set(nodeOne, secondNode);
        switchedList.set(nodeTwo, firstNode);
        return switchedList;
    }

    protected List<LNode> copyOfNodesInLayer(final int layerIndex) {
        return new ArrayList<LNode>(graph.getLayers().get(layerIndex).getNodes());
    }

    protected List<LNode> copyOfSwitchOrderOfNodesInLayer(final int nodeOne, final int nodeTwo, final int layerIndex) {
        List<LNode> layer = copyOfNodesInLayer(layerIndex);
        return getCopyWithSwitchedOrder(nodeOne, nodeTwo, layer);
    }

    protected List<LNode> getCopyWithSwitchedOrder(final int nodeOne, final int nodeTwo, final List<LNode> layer) {
        LNode firstNode = layer.get(nodeOne);
        LNode secondNode = layer.get(nodeTwo);
        List<LNode> switchedList = new ArrayList<LNode>(layer);
        switchedList.set(nodeOne, secondNode);
        switchedList.set(nodeTwo, firstNode);
        return switchedList;
    }

    protected void eastWestEdgesFromTo(final int numberOfNodes, final LNode left, final LNode right) {
        for (int i = 0; i < numberOfNodes; i++) {
            eastWestEdgeFromTo(left, right);
        }
    }

    protected void eastWestEdgeFromTo(final LPort leftPort, final LNode rightNode) {
        addEdgeBetweenPorts(leftPort, addPortOnSide(rightNode, PortSide.WEST));
    }

    protected List<LPort> copyPortsInIndexOrder(final LNode node, final int... indices) {
        List<LPort> res = Lists.newArrayList();
        for (int i : indices) {
            res.add(node.getPorts().get(i));
        }
        return res;
    }

    protected void eastWestEdgeFromTo(final LNode left, final LPort right) {
        addEdgeBetweenPorts(addPortOnSide(left, PortSide.EAST), right);
    }

    public MockRandom getRandom() {
        return random;
    }

    public void setRandom(final MockRandom random) {
        this.random = random;
    }

    protected static class MockRandom extends Random {
        public MockRandom() {
            super();
        }

        private static final long serialVersionUID = 1L;
        private float current = 0;
        private double changeBy = 0.0001;
        private boolean nextBoolean = true;

        @Override
        public boolean nextBoolean() {
            return nextBoolean;
        }

        @Override
        public float nextFloat() {
            current += changeBy;
            return current;
        }

        @Override
        public double nextDouble() {
            return nextFloat();
        }

        public void setNextBoolean(final boolean nextBoolean) {
            this.nextBoolean = nextBoolean;
        }

        public void setChangeBy(final double changeBy) {
            this.changeBy = changeBy;
        }

        public void setCurrent(final float current) {
            this.current = current;
        }

    }

    protected <T> void setOnAllGraphs(final IProperty<T> prop, final T val, final LGraph graph) {
        graph.setProperty(prop, val);
        for (LNode node : Iterables.concat(graph)) {
            LGraph nestedGraph = node.getNestedGraph();
            if (nestedGraph != null) {
                setOnAllGraphs(prop, val, nestedGraph);
            }
        }
    }

}
