/*******************************************************************************
 * Copyright (c) 2015 Kiel University and others.
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

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LNode.NodeType;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.p3order.counting.CrossMinUtil;
import org.eclipse.elk.core.options.PortSide;

/**
 * Counts the crossings caused by the order of north south port dummies when their respective normal
 * node in the same layer has a fixed port order. Also counts crossings between north south edges
 * and long edge dummies.
 * 
 * @author alan
 */
public class NorthSouthEdgeNeighbouringNodeCrossingsCounter {

    private int upperLowerCrossings;
    private int lowerUpperCrossings;
    private final Map<LPort, Integer> portPositions;
    private final LNode[] layer;

    /**
     * Creates a counter for north south port crossings.
     * 
     * @param nodes
     *            the order of nodes in the layer in question.
     */
    public NorthSouthEdgeNeighbouringNodeCrossingsCounter(final LNode[] nodes) {
        layer = nodes;
        portPositions = new HashMap<LPort, Integer>();
        initializePortPositions();
    }

    /**
     * Since accessing the index of a port is linear, the positions are saved. To prevent dependency
     * problems, they are saved in a field of this Object. Ports are numbered as they are in the
     * list returned by getPorts().
     */
    private void initializePortPositions() {
        for (LNode node : layer) {
            setPortIdsOn(node, PortSide.SOUTH);
            setPortIdsOn(node, PortSide.NORTH);
        }
    }

    private void setPortIdsOn(final LNode node, final PortSide side) {
        Iterable<LPort> ports = CrossMinUtil.inNorthSouthEastWestOrder(node, side);
        int portId = 0;
        for (LPort port : ports) {
            portPositions.put(port, portId++);
        }
    }

    /**
     * Counts north south port crossings and crossings between north south ports and dummy nodes,
     * for uppperNode and lowerNode.
     * 
     * @param upperNode
     *            The first node.
     * @param lowerNode
     *            The second node.
     */
    public void countCrossings(final LNode upperNode, final LNode lowerNode) {
        upperLowerCrossings = 0;
        lowerUpperCrossings = 0;

        processIfTwoNorthSouthNodes(upperNode, lowerNode);

        processIfNorthSouthLongEdgeDummyCrossing(upperNode, lowerNode);

        processIfNormalNodeWithNSPortsAndLongEdgeDummy(upperNode, lowerNode);
    }

    private void processIfTwoNorthSouthNodes(final LNode upperNode, final LNode lowerNode) {
        if (isNorthSouth(upperNode) && isNorthSouth(lowerNode) && !haveDifferentOrigins(upperNode, lowerNode)) {
            if (isNorthOfNormalNode(upperNode)) {
                countCrossingsOfTwoNorthSouthDummies(upperNode, lowerNode);
            } else {
                countCrossingsOfTwoNorthSouthDummies(lowerNode, upperNode);
            }
        }
    }

    private void countCrossingsOfTwoNorthSouthDummies(final LNode furtherFromNormalNode,
            final LNode closerToNormalNode) {
        if (originPortPositionOf(furtherFromNormalNode) > originPortPositionOf(closerToNormalNode)) {
            List<LPort> closerEastPorts = closerToNormalNode.getPortSideView(PortSide.EAST);
            upperLowerCrossings = closerEastPorts.isEmpty() ? 0 : closerEastPorts.get(0).getDegree();
            List<LPort> furtherWestPorts = furtherFromNormalNode.getPortSideView(PortSide.WEST);
            lowerUpperCrossings = furtherWestPorts.isEmpty() ? 0 : furtherWestPorts.get(0).getDegree();
        } else {
            List<LPort> closerWestPorts = closerToNormalNode.getPortSideView(PortSide.WEST);
            upperLowerCrossings = closerWestPorts.isEmpty() ? 0 : closerWestPorts.get(0).getDegree();
            List<LPort> furtherEastPorts = furtherFromNormalNode.getPortSideView(PortSide.EAST);
            lowerUpperCrossings = furtherEastPorts.isEmpty() ? 0 : furtherEastPorts.get(0).getDegree();
        }
    }

    private void processIfNorthSouthLongEdgeDummyCrossing(final LNode upperNode,
            final LNode lowerNode) {
        if (isNorthSouth(upperNode) && isLongEdgeDummy(lowerNode)) {
            if (isNorthOfNormalNode(upperNode)) {
                upperLowerCrossings = 1;
            } else {
                lowerUpperCrossings = 1;
            }
        } else if (isNorthSouth(lowerNode) && isLongEdgeDummy(upperNode)) {
            if (isNorthOfNormalNode(lowerNode)) {
                lowerUpperCrossings = 1;
            } else {
                upperLowerCrossings = 1;
            }
        }
    }

    private void processIfNormalNodeWithNSPortsAndLongEdgeDummy(final LNode upperNode,
            final LNode lowerNode) {

        if (isNormal(upperNode) && isLongEdgeDummy(lowerNode)) {
            upperLowerCrossings = numberOfNorthSouthEdges(upperNode, PortSide.SOUTH);
            lowerUpperCrossings = numberOfNorthSouthEdges(upperNode, PortSide.NORTH);
        }
        if (isNormal(lowerNode) && isLongEdgeDummy(upperNode)) {
            upperLowerCrossings = numberOfNorthSouthEdges(lowerNode, PortSide.NORTH);
            lowerUpperCrossings = numberOfNorthSouthEdges(lowerNode, PortSide.SOUTH);
        }
    }

    private int numberOfNorthSouthEdges(final LNode node, final PortSide side) {
        int numberOfEdges = 0;
        for (LPort port : node.getPortSideView(side)) {
            numberOfEdges += hasConnectedNorthSouthEdge(port) ? 1 : 0;
        }
        return numberOfEdges;
    }

    private boolean hasConnectedNorthSouthEdge(final LPort port) {
        return port.getProperty(InternalProperties.PORT_DUMMY) != null;
    }

    private boolean haveDifferentOrigins(final LNode upperNode, final LNode lowerNode) {
        return originOf(upperNode) != originOf(lowerNode);
    }

    private int originPortPositionOf(final LNode node) {
        LPort origin = originPortOf(node);
        final LPort port = origin;
        return portPositions.get(port);
    }

    private LPort originPortOf(final LNode node) {
        LPort port = node.getPorts().get(0);
        LPort origin = (LPort) port.getProperty(InternalProperties.ORIGIN);
        return origin;
    }

    private boolean isNorthOfNormalNode(final LNode upperNode) {
        return originPortOf(upperNode).getSide() == PortSide.NORTH;
    }

    private LNode originOf(final LNode node) {
        return (LNode) node.getProperty(InternalProperties.ORIGIN);
    }

    private boolean isLongEdgeDummy(final LNode node) {
        return node.getType() == NodeType.LONG_EDGE;
    }

    private boolean isNorthSouth(final LNode node) {
        return node.getType() == NodeType.NORTH_SOUTH_PORT;
    }

    private boolean isNormal(final LNode node) {
        return node.getType() == NodeType.NORMAL;
    }

    /**
     * Get crossing count.
     * 
     * @return the crossings between when ordered upper - lower.
     */
    public int getUpperLowerCrossings() {
        return upperLowerCrossings;
    }

    /**
     * 
     * Get crossing count.
     * 
     * @return the crossings between when ordered lower - upper.
     */
    public int getLowerUpperCrossings() {
        return lowerUpperCrossings;
    }

}
