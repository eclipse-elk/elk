/*******************************************************************************
 * Copyright (c) 2013, 2020 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.layered;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LNode.NodeType;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.core.options.EdgeRouting;
import org.eclipse.elk.core.options.PortConstraints;
import org.eclipse.elk.core.options.PortSide;
import org.eclipse.elk.graph.properties.MapPropertyHolder;

import com.google.common.collect.Lists;

/**
 * Use to create test graphs. <b>Caution:</b> Layout algorithm assumes the ports to be ordered in a clockwise manner.
 * You must think about this yourself when constructing a test graph. This means that the methods for creating edges
 * cannot be used in every case.
 */
public class TestGraphCreator {
    private int portId = 0;
    private int nodeId = 0;
    /** the graph. */
    protected final LGraph graph;

    /**
     * @return the graph
     */
    protected LGraph getGraph() {
        setUpIds();
        graph.setProperty(LayeredOptions.EDGE_ROUTING, EdgeRouting.ORTHOGONAL);
        return graph;
    }

    private int edgeId = 0;
    protected MockRandom random;

    /**
     * Makes a fancy test graph creator.
     */
    public TestGraphCreator() {
        graph = new LGraph();
        graph.setProperty(LayeredOptions.EDGE_ROUTING, EdgeRouting.ORTHOGONAL);
        random = new MockRandom();
        graph.setProperty(InternalProperties.RANDOM, random);
    }

    /**
     * Creates empty graph.
     * 
     * @return return empty graph.
     */
    public LGraph getEmptyGraph() {
        setUpIds();
        graph.setProperty(LayeredOptions.EDGE_ROUTING, EdgeRouting.ORTHOGONAL);
        return graph;
    }

    /**
     * Creates two nodes with no connection between them.
     * 
     * @return graph with two nodes with no connection between them.
     */
    public LGraph getTwoNodesNoConnectionGraph() {
        Layer layer = makeLayer();
        addNodeToLayer(layer);
        addNodeToLayer(layer);
        setUpIds();
        return graph;
    }

    protected void setUpIds() {
        int lId = 0;
        for (Layer l : graph) {
            l.id = lId++;
            int i = 0;
            for (LNode n : l) {
                n.id = i++;
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
        Layer leftLayer = makeLayer();
        Layer rightLayer = makeLayer();

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
        Layer leftLayer = makeLayer();
        Layer rightLayer = makeLayer();

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
        Layer leftLayer = makeLayer();
        Layer rightLayer = makeLayer();

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
        Layer layer = makeLayer();
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
        Layer leftLayer = makeLayer();
        Layer middleLayer = makeLayer();
        Layer rightLayer = makeLayer();

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
        Layer leftLayer = makeLayer();
        Layer rightLayer = makeLayer();

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
        Layer leftLayer = makeLayer();
        Layer rightLayer = makeLayer();

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
        Layer leftLayer = makeLayer();
        Layer rightLayer = makeLayer();

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
        LGraph selfLoopCrossGraph = getCrossFormedGraph();
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
        return selfLoopCrossGraph;
    }

    private void selfLoopOn(final LNode node, final PortSide side) {
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
     * Port order not fixed.
     * 
     * @return Graph of the form above.
     */
    public LGraph getMoreComplexThreeLayerGraph() {
        Layer leftLayer = makeLayer();
        Layer middleLayer = makeLayer();
        Layer rightLayer = makeLayer();

        LNode[] leftNodes = addNodesToLayer(3, leftLayer);
        LNode[] middleNodes = addNodesToLayer(2, middleLayer);
        LNode[] rightNodes = addNodesToLayer(3, rightLayer);

        eastWestEdgeFromTo(leftNodes[0], middleNodes[0]);
        eastWestEdgeFromTo(leftNodes[1], middleNodes[0]);
        eastWestEdgeFromTo(leftNodes[1], middleNodes[1]);
        eastWestEdgeFromTo(leftNodes[2], middleNodes[0]);
        eastWestEdgeFromTo(middleNodes[0], rightNodes[0]);
        eastWestEdgeFromTo(middleNodes[0], rightNodes[1]);
        eastWestEdgeFromTo(middleNodes[0], rightNodes[1]);
        eastWestEdgeFromTo(middleNodes[1], rightNodes[0]);
        eastWestEdgeFromTo(middleNodes[1], rightNodes[2]);

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
        Layer leftLayer = makeLayer();
        Layer rightLayer = makeLayer();

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
        Layer leftLayer = makeLayer();
        Layer rightLayer = makeLayer();

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

        eastWestEdgeFromTo(leftNodes[0], middleNodes[1]);
        eastWestEdgeFromTo(leftNodes[1], middleNodes[0]);
        eastWestEdgeFromTo(middleNodes[0], rightNodes[0]);
        eastWestEdgeFromTo(middleNodes[1], rightNodes[1]);

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
        Layer layer = makeLayer();
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
        Layer layer = makeLayer();
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

        eastWestEdgeFromTo(leftNodes[0], middleNodes[0]);
        eastWestEdgeFromTo(leftNodes[0], middleNodes[0]);
        eastWestEdgeFromTo(leftNodes[1], middleNodes[1]);
        eastWestEdgeFromTo(leftNodes[1], middleNodes[1]);
        eastWestEdgeFromTo(middleNodes[0], rightNodes[1]);
        eastWestEdgeFromTo(middleNodes[1], rightNodes[0]);

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

        eastWestEdgeFromTo(leftNodes[0], rightNodes[1]);
        eastWestEdgeFromTo(leftNodes[1], rightNodes[0]);

        LPort southPort = addPortOnSide(rightNodes[1], PortSide.SOUTH);
        LPort northPort = addPortOnSide(rightNodes[2], PortSide.NORTH);

        addEdgeBetweenPorts(southPort, northPort);

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
        Layer leftLayer = makeLayer();
        LNode[] leftNodes = addNodesToLayer(2, leftLayer);
        Layer rightLayer = makeLayer();
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
        Layer leftLayer = makeLayer();
        LNode[] leftNodes = addNodesToLayer(2, leftLayer);
        Layer rightLayer = makeLayer();
        LNode[] rightNodes = addNodesToLayer(3, rightLayer);

        addInLayerEdge(rightNodes[0], rightNodes[2], PortSide.WEST);
        eastWestEdgeFromTo(leftNodes[0], rightNodes[1]);
        eastWestEdgeFromTo(leftNodes[0], rightNodes[1]);
        setUpIds();
        return graph;
    }

    /**
     * <pre>
     * ----*
     * |  /
     * *++-*
     *  /\ 
     * *--+* 
     *     *
     * </pre>
     * 
     * Should have four crossings.
     * 
     * @return Graph of the form above.
     */
    public LGraph shouldSwitchThreeTimesGraph() {
        LNode[] leftNodes = addNodesToLayer(2, makeLayer());
        LNode[] rightNodes = addNodesToLayer(4, makeLayer());

        eastWestEdgeFromTo(leftNodes[1], rightNodes[0]);
        eastWestEdgeFromTo(leftNodes[1], rightNodes[2]);

        eastWestEdgeFromTo(leftNodes[0], rightNodes[0]);
        eastWestEdgeFromTo(leftNodes[0], rightNodes[1]);
        eastWestEdgeFromTo(leftNodes[0], rightNodes[3]);

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
        LNode[] leftNodes = addNodesToLayer(2, makeLayer());
        LNode[] rightNodes = addNodesToLayer(4, makeLayer());

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
        Layer leftLayer = makeLayer();
        Layer rightLayer = makeLayer();

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
        Layer leftLayer = makeLayer();
        Layer rightLayer = makeLayer();

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
        Layer leftLayer = makeLayer();
        Layer rightLayer = makeLayer();

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
        Layer leftLayer = makeLayer();
        Layer rightLayer = makeLayer();

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
     * Returns the nodes in the graph as two-dimensional array of LNodes.
     * 
     * @return graph as LNode[][].
     */
    public LNode[][] getCurrentOrder() {
        LNode[][] nodeOrder = new LNode[graph.getLayers().size()][];
        List<Layer> layers = graph.getLayers();
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

    protected void setAsNorthSouthNode(final LNode node) {
        node.setType(NodeType.NORTH_SOUTH_PORT);
    }

    protected void addNorthSouthEdge(final PortSide side, final LNode nodeWithNSPorts, final LNode northSouthDummy,
            final LNode nodeWithEastWestPorts, final boolean nodeWithEastWestPortsIsOrigin) {
        boolean normalNodeEastOfNsPortNode =
                nodeWithEastWestPorts.getLayer().getIndex() < nodeWithNSPorts.getLayer().getIndex();
        PortSide direction = normalNodeEastOfNsPortNode ? PortSide.WEST : PortSide.EAST;

        PortSide targetNodePortSide = direction == PortSide.WEST ? PortSide.EAST : PortSide.WEST;
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

    protected void setInLayerOrderConstraint(final LNode thisNode, final LNode beforeThisNode) {
        List<LNode> scndNodeAsList = Lists.newArrayList(beforeThisNode);
        thisNode.setProperty(InternalProperties.IN_LAYER_SUCCESSOR_CONSTRAINTS, scndNodeAsList);
    }

    protected void setAsLongEdgeDummy(final LNode node) {
        node.setType(NodeType.LONG_EDGE);
    }

    protected void setPortOrderFixed(final LNode node) {
        node.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_ORDER);
    }

    protected Layer[] makeLayers(final int amount) {
        Layer[] layers = new Layer[amount];
        for (int i = 0; i < layers.length; i++) {
            layers[i] = makeLayer();
        }
        return layers;
    }

    protected MapPropertyHolder setFixedOrderConstraint(final LNode node) {
        return node.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_ORDER);
    }

    protected void addInLayerEdge(final LNode nodeOne, final LNode nodeTwo, final PortSide portSide) {
        LPort portOne = addPortOnSide(nodeOne, portSide);
        LPort portTwo = addPortOnSide(nodeTwo, portSide);
        addEdgeBetweenPorts(portOne, portTwo);
    }

    protected LNode[] addNodesToLayer(final int amountOfNodes, final Layer leftLayer) {
        LNode[] nodes = new LNode[amountOfNodes];
        for (int j = 0; j < amountOfNodes; j++) {
            nodes[j] = addNodeToLayer(leftLayer);
        }
        return nodes;
    }

    protected Layer makeLayer() {
        return makeLayer(graph);
    }

    protected Layer makeLayer(final LGraph graph) {
        List<Layer> layers = graph.getLayers();
        Layer layer = new Layer(graph);
        layers.add(layer);
        return layer;
    }

    protected LNode addNodeToLayer(final Layer layer) {
        LNode node = new LNode(graph);
        node.setType(NodeType.NORMAL);
        node.setProperty(InternalProperties.IN_LAYER_LAYOUT_UNIT, node);
        node.setLayer(layer);
        node.id = nodeId++;
        return node;
    }

    protected void eastWestEdgeFromTo(final LNode left, final LNode right) {
        LPort leftPort = addPortOnSide(left, PortSide.EAST);
        LPort rightPort = addPortOnSide(right, PortSide.WEST);
        addEdgeBetweenPorts(leftPort, rightPort);
    }

    protected void eastWestEdgeFromTo(final LPort leftPort, final LNode right) {
        LPort rightPort = addPortOnSide(right, PortSide.WEST);
        addEdgeBetweenPorts(leftPort, rightPort);
    }

    protected void addEdgeBetweenPorts(final LPort from, final LPort to) {
        LEdge edge = new LEdge();
        edge.setSource(from);
        edge.setTarget(to);
        edge.id = edgeId++;
    }

    protected LPort addPortOnSide(final LNode node, final PortSide portSide) {
        LPort port = addPortTo(node);
        port.setSide(portSide);
        return port;
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
        LNode[] nodes = addNodesToLayer(3, makeLayer());
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
        LNode leftNode = addNodeToLayer(makeLayer());
        setFixedOrderConstraint(leftNode);
        LNode[] rightNodes = addNodesToLayer(3, makeLayer());
        eastWestEdgeFromTo(leftNode, rightNodes[2]);
        eastWestEdgeFromTo(leftNode, rightNodes[2]);
        eastWestEdgeFromTo(leftNode, rightNodes[1]);
        addInLayerEdge(rightNodes[0], rightNodes[2], PortSide.WEST);
        setUpIds();
        return graph;
    }

    protected final class MockRandom extends Random {
        private static final long serialVersionUID = 1L;
        private boolean nextBoolean = true;
        private double changeBy = 0.01;
        private double currentDouble = 0.01;

        public MockRandom() {
        }

        @Override
        public boolean nextBoolean() {
            return nextBoolean;
        }

        @Override
        public double nextDouble() {
            return currentDouble += changeBy;
        }

        @Override
        public float nextFloat() {
            return (float) nextDouble();
        }

        public void setNextBoolean(final boolean b) {
            nextBoolean = b;
        }

        public void setChangeBy(final double d) {
            changeBy = d;
        }
    }

    protected <T> List<T> getListInIndexOrder(final List<T> li, final int... is) {
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

    protected <T> List<T> switchOrderInList(final int i, final int j, final List<T> list) {
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

}
