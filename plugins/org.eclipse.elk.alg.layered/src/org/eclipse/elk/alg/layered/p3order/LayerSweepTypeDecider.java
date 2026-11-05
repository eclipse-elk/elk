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

import java.util.ArrayList;
import java.util.List;

import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LNode.NodeType;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.p3order.counting.IInitializable;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.core.options.PortSide;

import org.eclipse.elk.alg.layered.graph.LGraphUtil;

/**
 * In order to decide whether to sweep into the graph or not, we compare the number of paths to nodes whose position is
 * decided on by random decisions to the number of paths to nodes whose position depends on cross-hierarchy edges. By
 * calculating (pathsToRandom - pathsToHierarchical) / allPaths, this value will always be between -1 (many cross
 * hierarchical paths) and +1 (many random paths). By setting the boundary
 * CROSSING_MINIMIZATION_HIERARCHICAL_SWEEPINESS, we can choose how likely it is to be hierarchical or more bottom up.
 * 
 * <p>
 * Must be initialized using {@link IInitializable#init(java.util.List)}!
 * </p>
 */
public class LayerSweepTypeDecider implements IInitializable {

    private NodeInfo[][] nodeInfo;
    private GraphInfoHolder graphData;

    /**
     * Creates LayerSweepTypeDecider for deciding whether to sweep into graphs or not.
     * <p>
     * Must be initialized using {@link IInitializable#init(java.util.List)}!
     * </p>
     * 
     * @param graphData
     *            the graph holder object for auxiliary graph information needed in crossing minimization
     */
    public LayerSweepTypeDecider(final GraphInfoHolder graphData) {
        this.graphData = graphData;
        nodeInfo = new NodeInfo[graphData.currentNodeOrder().length][];
    }

    /**
     * Decide whether to use bottom up or cross-hierarchical sweep method.
     * 
     * @return decision
     */
    public boolean useBottomUp() {
        double boundary = graphData.lGraph().getProperty(LayeredOptions.CROSSING_MINIMIZATION_HIERARCHICAL_SWEEPINESS);
        if (bottomUpForced(boundary) || rootNode() || fixedPortOrder() || fewerThanTwoInOutEdges()) {
            return true;
        }

        if (graphData.crossMinDeterministic()) {
            return false;
        }

        int pathsToRandom = 0;
        int pathsToHierarchical = 0;

        List<LNode> nsPortDummies = new ArrayList<>();
        for (LNode[] layer : graphData.currentNodeOrder()) {
            for (LNode node : layer) {
                // We must visit all sources of edges first, so we collect north south dummies for later.
                if (isNorthSouthDummy(node)) {
                    nsPortDummies.add(node);
                    continue;
                }

                NodeInfo currentNode = nodeInfoFor(node);
                // Check for hierarchical port dummies or random influence.
                if (isExternalPortDummy(node)) {
                    currentNode.hierarchicalInfluence = 1;
                    if (isEasternDummy(node)) {
                        pathsToHierarchical += currentNode.connectedEdges;
                    }
                } else if (hasNoWesternPorts(node)) {
                    currentNode.randomInfluence = 1;
                } else if (hasNoEasternPorts(node)) {
                    pathsToRandom += currentNode.connectedEdges;
                }

                // Increase counts of paths by the number outgoing edges times the influence
                // and transfer information to targets.
                for (LEdge edge : node.getOutgoingEdges()) {
                    pathsToRandom += currentNode.randomInfluence;
                    pathsToHierarchical += currentNode.hierarchicalInfluence;
                    transferInfoToTarget(currentNode, edge);
                }

                // Do the same for north/south dummies: Increase counts of paths by the number outgoing edges times the
                // influence and transfer information to dummies.
                Iterable<LPort> northSouthPorts =
                        LGraphUtil.concat(node.getPortSideView(PortSide.NORTH), node.getPortSideView(PortSide.SOUTH));
                for (LPort port : northSouthPorts) {
                    LNode nsDummy = port.getProperty(InternalProperties.PORT_DUMMY);
                    if (nsDummy != null) {
                        pathsToRandom += currentNode.randomInfluence;
                        pathsToHierarchical += currentNode.hierarchicalInfluence;
                        transferInfoTo(currentNode, nsDummy);
                    }
                }
            }

            // Now process nsPortDummies
            for (LNode node : nsPortDummies) {
                NodeInfo currentNode = nodeInfoFor(node);
                for (LEdge edge : node.getOutgoingEdges()) {
                    pathsToRandom += currentNode.randomInfluence;
                    pathsToHierarchical += currentNode.hierarchicalInfluence;
                    transferInfoToTarget(currentNode, edge);
                }
            }
            nsPortDummies.clear();
        }

        double allPaths = pathsToRandom + pathsToHierarchical;
        double normalized = allPaths == 0 ? Double.POSITIVE_INFINITY : (pathsToRandom - pathsToHierarchical) / allPaths;
        return normalized >= boundary;
    }

    private boolean fixedPortOrder() {
        return graphData.parent().getProperty(LayeredOptions.PORT_CONSTRAINTS).isOrderFixed();
    }

    private void transferInfoToTarget(final NodeInfo currentNode, final LEdge edge) {
        LNode target = targetNode(edge);
        transferInfoTo(currentNode, target);
    }

    private void transferInfoTo(final NodeInfo currentNode, final LNode target) {
        NodeInfo targetNodeInfo = nodeInfoFor(target);
        targetNodeInfo.transfer(currentNode);
        targetNodeInfo.connectedEdges++;
    }

    /**
     * For a single node, collects number of paths to nodes with random or hierarchical nodes.
     * 
     */
    static class NodeInfo {
        private int connectedEdges;
        private int hierarchicalInfluence;
        private int randomInfluence;

        public void transfer(final NodeInfo nodeInfo) {
            hierarchicalInfluence += nodeInfo.hierarchicalInfluence;
            randomInfluence += nodeInfo.randomInfluence;
            connectedEdges += nodeInfo.connectedEdges;
        }

        @Override
        public String toString() {
            return "NodeInfo [connectedEdges=" + connectedEdges + ", hierarchicalInfluence=" + hierarchicalInfluence
                    + ", randomInfluence=" + randomInfluence + "]";
        }
    }

    private boolean fewerThanTwoInOutEdges() {
        return graphData.parent().getPortSideView(PortSide.EAST).size() < 2
                && graphData.parent().getPortSideView(PortSide.WEST).size() < 2;
    }

    private boolean rootNode() {
        return !graphData.hasParent();
    }

    private boolean bottomUpForced(final double boundary) {
        return boundary < -1;
    }

    private LNode targetNode(final LEdge edge) {
        return edge.getTarget().getNode();
    }

    private boolean hasNoEasternPorts(final LNode node) {
        List<LPort> eastPorts = node.getPortSideView(PortSide.EAST);
        return eastPorts.isEmpty() || !java.util.stream.StreamSupport.stream(eastPorts.spliterator(), false).anyMatch(p -> p.getConnectedEdges().iterator().hasNext());
    }

    private boolean hasNoWesternPorts(final LNode node) {
        List<LPort> westPorts = node.getPortSideView(PortSide.WEST);
        return westPorts.isEmpty() || !java.util.stream.StreamSupport.stream(westPorts.spliterator(), false).anyMatch(p -> p.getConnectedEdges().iterator().hasNext());
    }

    private boolean isExternalPortDummy(final LNode node) {
        return node.getType() == NodeType.EXTERNAL_PORT;
    }

    private boolean isNorthSouthDummy(final LNode node) {
        return node.getType() == NodeType.NORTH_SOUTH_PORT;
    }

    private boolean isEasternDummy(final LNode node) {
        return originPort(node).getSide() == PortSide.EAST;
    }

    private LPort originPort(final LNode node) {
        return (LPort) node.getProperty(InternalProperties.ORIGIN);
    }

    private NodeInfo nodeInfoFor(final LNode n) {
        return nodeInfo[n.getLayer().id][n.id];
    }

    @Override
    public void initAtLayerLevel(final int l, final LNode[][] nodeOrder) {
        nodeOrder[l][0].getLayer().id = l;
        nodeInfo[l] = new NodeInfo[nodeOrder[l].length];
    }

    @Override
    public void initAtNodeLevel(final int l, final int n, final LNode[][] nodeOrder) {
        LNode node = nodeOrder[l][n];
        node.id = n;
        nodeInfo[l][n] = new NodeInfo();
    }
}
