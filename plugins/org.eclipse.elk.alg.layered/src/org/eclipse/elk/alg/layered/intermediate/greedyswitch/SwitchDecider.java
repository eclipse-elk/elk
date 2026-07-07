/*******************************************************************************
 * Copyright (c) 2014, 2015 Kiel University and others.
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
package org.eclipse.elk.alg.layered.intermediate.greedyswitch;

import java.util.List;

import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LNode.NodeType;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.p3order.GraphInfoHolder;
import org.eclipse.elk.alg.layered.p3order.counting.CrossingsCounter;
import org.eclipse.elk.core.options.PortSide;
import org.eclipse.elk.core.util.Pair;

/**
 * This class decides whether two neighboring nodes should be switched. There are two variants:
 * <ul>
 * <li>OneSided – The traditional risky way: The decider checks if a switch would reduce crossings
 * on the given side of the layer whose nodes are to be switched.
 * <li>TwoSided – The faithless way: The decider checks if a switch would reduce crossings on both
 * sides of the layer whose nodes are to be switched.
 * </ul>
 */
public final class SwitchDecider {
    private final LNode[] freeLayer;
    private final CrossingsCounter leftInLayerCounter;
    private final CrossingsCounter rightInLayerCounter;
    private final NorthSouthEdgeNeighbouringNodeCrossingsCounter northSouthCounter;
    private final CrossingMatrixFiller crossingMatrixFiller;
    private GraphInfoHolder graphData;
    private CrossingsCounter parentCrossCounter;
    private boolean countCrossingsCausedByPortSwitch;

    /**
     * Creates SwitchDecider.
     * 
     * @param freeLayerIndex
     *            The freeLayer to switch in.
     * @param graph
     *            The graph as LNode[][]
     * @param crossingMatrixFiller
     *            the crossing matrix filler
     * @param portPositions
     *            port position array
     * @param oneSided
     *            whether greedy mode is one sided or two-sided (not one-sided)
     * @param graphData
     *            collected graphData
     */
    public SwitchDecider(final int freeLayerIndex, final LNode[][] graph,
            final CrossingMatrixFiller crossingMatrixFiller, final int[] portPositions,
            final GraphInfoHolder graphData, final boolean oneSided) {
        this.crossingMatrixFiller = crossingMatrixFiller;
        this.graphData = graphData;
        if (freeLayerIndex >= graph.length) {
            throw new IndexOutOfBoundsException("Greedy SwitchDecider: Free layer not in graph.");
        }
        freeLayer = graph[freeLayerIndex];

        leftInLayerCounter = new CrossingsCounter(portPositions);
        leftInLayerCounter.initPortPositionsForInLayerCrossings(freeLayer, PortSide.WEST);
        rightInLayerCounter = new CrossingsCounter(portPositions);
        rightInLayerCounter.initPortPositionsForInLayerCrossings(freeLayer, PortSide.EAST);
        northSouthCounter = new NorthSouthEdgeNeighbouringNodeCrossingsCounter(freeLayer);
        countCrossingsCausedByPortSwitch = !oneSided && graphData.hasParent() && !graphData.dontSweepInto()
                && freeLayer[0].getType() == NodeType.EXTERNAL_PORT;
        if (countCrossingsCausedByPortSwitch) {
            initParentCrossingsCounters(freeLayerIndex, graph.length);
        }
    }

    private void initParentCrossingsCounters(final int freeLayerIndex, final int length) {
        GraphInfoHolder parentGraphData = graphData.parentGraphData();
        LNode[][] parentNodeOrder = parentGraphData.currentNodeOrder();
        int[] portPos = parentGraphData.portPositions();
        parentCrossCounter = new CrossingsCounter(portPos);
        int parentNodeLayerPos = graphData.parent().getLayer().id;
        LNode[] leftLayer = parentNodeLayerPos > 0 ? parentNodeOrder[parentNodeLayerPos - 1] : new LNode[0];
        LNode[] middleLayer = parentNodeOrder[parentNodeLayerPos];
        LNode[] rightLayer = parentNodeLayerPos < parentNodeOrder.length - 1 ? parentNodeOrder[parentNodeLayerPos + 1]
                : new LNode[0];
        boolean rightMostLayer = freeLayerIndex == length - 1;
        if (rightMostLayer) {
            parentCrossCounter.initForCountingBetween(middleLayer, rightLayer);
        } else {
            parentCrossCounter.initForCountingBetween(leftLayer, middleLayer);
        }
    }

    /**
     * Notifies in-layer counter of node switch for efficiency reasons.
     * 
     * @param upperNode
     *            a node
     * @param lowerNode
     *            a node
     */
    public void notifyOfSwitch(final LNode upperNode, final LNode lowerNode) {
        leftInLayerCounter.switchNodes(upperNode, lowerNode, PortSide.WEST);
        rightInLayerCounter.switchNodes(upperNode, lowerNode, PortSide.EAST);
        if (countCrossingsCausedByPortSwitch) {
            LPort upperPort = (LPort) upperNode.getProperty(InternalProperties.ORIGIN);
            LPort lowerPort = (LPort) lowerNode.getProperty(InternalProperties.ORIGIN);
            parentCrossCounter.switchPorts(upperPort, lowerPort);
        }
    }

    /**
     * @return whether switching the nodes represented by the indices would reduce the number of crossings.
     */
    public boolean doesSwitchReduceCrossings(final int upperNodeIndex, final int lowerNodeIndex) {
        
        if (constraintsPreventSwitch(upperNodeIndex, lowerNodeIndex)) {
            return false;
        }

        LNode upperNode = freeLayer[upperNodeIndex];
        LNode lowerNode = freeLayer[lowerNodeIndex];
        
        Pair<Integer, Integer> leftInlayer =
                leftInLayerCounter.countInLayerCrossingsBetweenNodesInBothOrders(upperNode, lowerNode, PortSide.WEST);
        Pair<Integer, Integer> rightInlayer =
                rightInLayerCounter.countInLayerCrossingsBetweenNodesInBothOrders(upperNode, lowerNode, PortSide.EAST);
        northSouthCounter.countCrossings(upperNode, lowerNode);
        int upperLowerCrossings =
                crossingMatrixFiller.getCrossingMatrixEntry(upperNode, lowerNode)
                        + leftInlayer.getFirst() + rightInlayer.getFirst()
                        + northSouthCounter.getUpperLowerCrossings();
        int lowerUpperCrossings =
                crossingMatrixFiller.getCrossingMatrixEntry(lowerNode, upperNode)
                        + leftInlayer.getSecond() + rightInlayer.getSecond()
                        + northSouthCounter.getLowerUpperCrossings();

        if (countCrossingsCausedByPortSwitch) {
            LPort upperPort = (LPort) upperNode.getProperty(InternalProperties.ORIGIN);
            LPort lowerPort = (LPort) lowerNode.getProperty(InternalProperties.ORIGIN);
            Pair<Integer, Integer> crossingNumbers =
                    parentCrossCounter.countCrossingsBetweenPortsInBothOrders(upperPort, lowerPort);
            upperLowerCrossings += crossingNumbers.getFirst();
            lowerUpperCrossings += crossingNumbers.getSecond();
        }

        return upperLowerCrossings > lowerUpperCrossings;
    }


    /**
     * Check if in layer {@link InternalProperties.IN_LAYER_SUCCESSOR_CONSTRAINTS} or
     * {@link InternalProperties.IN_LAYER_LAYOUT_UNIT} constraints prevent a possible switch or if
     * the nodes are a normal node and a north south port dummy.
     * 
     * @param nodeIndex
     *            the index of the upper node, assuming a left-right order.
     * @param lowerNodeIndex
     *            the index of the lower node, assuming a left-right order.
     * @return true if constraints should prevent switching.
     */
    private boolean constraintsPreventSwitch(final int nodeIndex, final int lowerNodeIndex) {
        LNode upperNode = freeLayer[nodeIndex];
        LNode lowerNode = freeLayer[lowerNodeIndex];

        return haveSuccessorConstraints(upperNode, lowerNode)
                || haveLayoutUnitConstraints(upperNode, lowerNode)
                || areNormalAndNorthSouthPortDummy(upperNode, lowerNode);
    }

    private boolean haveSuccessorConstraints(final LNode upperNode, final LNode lowerNode) {
        List<LNode> constraints =
                upperNode.getProperty(InternalProperties.IN_LAYER_SUCCESSOR_CONSTRAINTS);
        boolean hasSuccessorConstraint =
                constraints != null && constraints.size() != 0 && constraints.contains(lowerNode);
        return hasSuccessorConstraint;
    }

    private boolean haveLayoutUnitConstraints(final LNode upperNode, final LNode lowerNode) {
        boolean neitherNodeIsLongEdgeDummy =
                upperNode.getType() != NodeType.LONG_EDGE
                        && lowerNode.getType() != NodeType.LONG_EDGE;

        // If upperNode and lowerNode are part of a layout unit not only containing themselves,
        // then the layout units must be equal for a switch to be allowed.
        LNode upperLayoutUnit = upperNode.getProperty(InternalProperties.IN_LAYER_LAYOUT_UNIT);
        LNode lowerLayoutUnit = lowerNode.getProperty(InternalProperties.IN_LAYER_LAYOUT_UNIT);
        
        boolean areInDifferentLayoutUnits = upperLayoutUnit != lowerLayoutUnit;

        // FIXME the following predicate is problematic, layout units are represented by a regular node,
        // thus 'upperNode' can be 'upperLayoutUnit' and still have more nodes in the layout unit
        boolean nodesHaveLayoutUnits =
                partOfMultiNodeLayoutUnit(upperNode, upperLayoutUnit)
                        || partOfMultiNodeLayoutUnit(lowerNode, lowerLayoutUnit);

        boolean upperNodeHasNorthernEdges = hasEdgesOnSide(upperNode, PortSide.NORTH);
        boolean lowerNodeHasSouthernEdges = hasEdgesOnSide(lowerNode, PortSide.SOUTH);

        // hotfix for #162, if north or south edges are present, there must be a layout unit
        nodesHaveLayoutUnits |= hasEdgesOnSide(upperNode, PortSide.SOUTH) || hasEdgesOnSide(lowerNode, PortSide.NORTH); 

        boolean hasLayoutUnitConstraint = (nodesHaveLayoutUnits && areInDifferentLayoutUnits)
                || (upperNodeHasNorthernEdges || lowerNodeHasSouthernEdges);

        return neitherNodeIsLongEdgeDummy && hasLayoutUnitConstraint;
    }

    private boolean hasEdgesOnSide(final LNode node, final PortSide side) {
        Iterable<LPort> ports = node.getPortSideView(side);
        for (LPort port : ports) {
            if (port.getProperty(InternalProperties.PORT_DUMMY) != null
                    || port.getConnectedEdges().iterator().hasNext()) {
                return true;
            }
        }
        return false;
    }

    private boolean partOfMultiNodeLayoutUnit(final LNode node, final LNode layoutUnit) {
        return layoutUnit != null && layoutUnit != node;
    }

    private boolean areNormalAndNorthSouthPortDummy(final LNode upperNode, final LNode lowerNode) {
        return isNorthSouthPortNode(upperNode) && isNormalNode(lowerNode)
                || isNorthSouthPortNode(lowerNode) && isNormalNode(upperNode);
    }

    private boolean isNormalNode(final LNode node) {
        return node.getType() == NodeType.NORMAL;
    }

    private boolean isNorthSouthPortNode(final LNode node) {
        return node.getType() == NodeType.NORTH_SOUTH_PORT;
    }


    /**
     * The side on which to count crossings for the one-sided SwitchDecider.
     * 
     * @author alan
     *
     */
    public enum CrossingCountSide {
        /** Consider crossings to the west of the free layer. */
        WEST,
        /** Consider crossings to the east of the free layer. */
        EAST
    }

}