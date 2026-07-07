/*******************************************************************************
 * Copyright (c) 2010, 2015 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.layered.p3order;

import java.util.Collections;
import java.util.List;

import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.PortType;
import org.eclipse.elk.alg.layered.p3order.counting.IInitializable;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.core.options.PortSide;

import com.google.common.collect.Lists;

/**
 * Calculates port ranks and distributes ports.
 * The <em>rank</em> of a port is a floating point number that represents its position
 * inside the containing layer. This depends on the node order of that layer and on the
 * port constraints of the nodes. Port ranks are used by {@link ICrossingMinimizationHeuristic}s
 * for calculating barycenter or median values for nodes. Furthermore, they are used in this
 * class for distributing the ports of nodes where the order of ports is not fixed,
 * which has to be done as the last step of each crossing minimization processor.
 * There are different ways to determine port ranks, therefore that is done in concrete subclasses.
 * <p>
 * Must be initialized using {@link IInitializable#init(List, LNode[][])}.
 *
 * @author cds
 * @author ima
 * @author msp
 */
public abstract class AbstractBarycenterPortDistributor implements ISweepPortDistributor, IInitializable {

    /** port ranks array in which the results of ranks calculation are stored. */
    private float[] portRanks;
    private float minBarycenter;
    private float maxBarycenter;
    private int[][] nodePositions;
    private float[] portBarycenter;
    private List<LPort> inLayerPorts;

    /**
     * Constructs a port distributor for the given array of port ranks. 
     * All ports are required to be assigned ids in the range of the given array.
     *
     * @param numLayers
     *            the number of layers in the graph.
     */
    public AbstractBarycenterPortDistributor(final int numLayers) {
        inLayerPorts = Lists.newArrayList();
        nodePositions = new int[numLayers][];
    }

    /**
     * Returns the array of port ranks.
     *
     * @return the array of port ranks
     */
    public float[] getPortRanks() {
        return portRanks;
    }


    // /////////////////////////////////////////////////////////////////////////////
    // Port Rank Assignment

    @Override
    public boolean distributePortsWhileSweeping(final LNode[][] nodeOrder, 
            final int currentIndex,
            final boolean isForwardSweep) {
        
        updateNodePositions(nodeOrder, currentIndex);
        LNode[] freeLayer = nodeOrder[currentIndex];
        PortSide side = isForwardSweep ? PortSide.WEST : PortSide.EAST;
        if (isNotFirstLayer(nodeOrder.length, currentIndex, isForwardSweep)) {
    
            LNode[] fixedLayer = nodeOrder[isForwardSweep ? currentIndex - 1 : currentIndex + 1];
            calculatePortRanks(fixedLayer, portTypeFor(isForwardSweep));
            for (LNode node : freeLayer) {
                distributePorts(node, side);
            }
    
            calculatePortRanks(freeLayer, portTypeFor(!isForwardSweep));
            for (LNode node : fixedLayer) {
                if (!hasNestedGraph(node)) {
                    distributePorts(node, side.opposed());
                }
            }
        } else {
            for (LNode node : freeLayer) {
                distributePorts(node, side);
            }
        }
        // Barycenter port distributor can not be used with always improving crossing minimization heuristics
        // which do not need to count.
        return false;
    }

    /**
     * Determine ranks for all ports of specific type in the given layer.
     * The ranks are written to the {@link #getPortRanks()} array.
     *
     * @param layer
     *            a layer as node array
     * @param portType
     *            the port type to consider
     */
    public final void calculatePortRanks(final LNode[] layer, final PortType portType) {
        float consumedRank = 0;
        for (int nodeIx = 0; nodeIx < layer.length; nodeIx++) {
            consumedRank += calculatePortRanks(layer[nodeIx], consumedRank, portType);
        }
    }

    /**
     * Assign port ranks for the input or output ports of the given node. If the node's port
     * constraints imply a fixed order, the ports are assumed to be pre-ordered in the usual way,
     * i.e. in clockwise order north - east - south - west.
     * The ranks are written to the {@link #getPortRanks()} array.
     *
     * @param node
     *            a node
     * @param rankSum
     *            the sum of ranks of preceding nodes in the same layer
     * @param type
     *            the port type to consider
     * @return the rank consumed by the given node; the following node's ranks start at
     *         {@code rankSum + consumedRank}
     * @see {@link org.eclipse.alg.layered.intermediate.PortListSorter} 
     */
    protected abstract float calculatePortRanks(LNode node, float rankSum, PortType type);


    // /////////////////////////////////////////////////////////////////////////////
    // Port Distribution

    private void distributePorts(final LNode node, final PortSide side) {
        if (!node.getProperty(LayeredOptions.PORT_CONSTRAINTS).isOrderFixed()) {
            // distribute ports in sweep direction and on north south side of node.
            distributePorts(node, node.getPorts(side));
            distributePorts(node, node.getPorts(PortSide.SOUTH));
            distributePorts(node, node.getPorts(PortSide.NORTH));
            // sort the ports by considering the side, type, and barycenter values
            sortPorts(node);
        }
    }

    /**
     * Distribute the ports of the given node by their sides, connected ports, and input or output
     * type.
     *
     * @param node
     *            node whose ports shall be sorted
     */
    private void distributePorts(final LNode node, final Iterable<LPort> ports) {
        inLayerPorts.clear();
        iteratePortsAndCollectInLayerPorts(node, ports);

        if (!inLayerPorts.isEmpty()) {
            calculateInLayerPortsBarycenterValues(node);
        }
    }

    private void iteratePortsAndCollectInLayerPorts(final LNode node, final Iterable<LPort> ports) {
        minBarycenter = 0.0f;
        maxBarycenter = 0.0f;

        // a float value large enough to ensure that barycenters of south ports work fine
        final float absurdlyLargeFloat = 2 * node.getLayer().getNodes().size() + 1;
        // calculate barycenter values for the ports of the node
        PortIteration: 
        for (LPort port : ports) {
            
            boolean northSouthPort = port.getSide() == PortSide.NORTH || port.getSide() == PortSide.SOUTH;
            float sum = 0;

            if (northSouthPort) {
                // Find the dummy node created for the port
                LNode portDummy = port.getProperty(InternalProperties.PORT_DUMMY);
                if (portDummy == null) {
                    continue;
                }

                sum += dealWithNorthSouthPorts(absurdlyLargeFloat, port, portDummy);

            } else {
                // add up all ranks of connected ports
                for (LEdge outgoingEdge : port.getOutgoingEdges()) {
                    LPort connectedPort = outgoingEdge.getTarget();
                    if (connectedPort.getNode().getLayer() == node.getLayer()) {
                        inLayerPorts.add(port);
                        continue PortIteration;
                    } else {
                        // outgoing edges go to the subsequent layer and are seen clockwise
                        sum += portRanks[connectedPort.id];
                    }
                }
                for (LEdge incomingEdge : port.getIncomingEdges()) {
                    LPort connectedPort = incomingEdge.getSource();
                    if (connectedPort.getNode().getLayer() == node.getLayer()) {
                        inLayerPorts.add(port);
                        continue PortIteration;
                    } else {
                        // incoming edges go to the preceding layer and are seen
                        // counter-clockwise
                        sum -= portRanks[connectedPort.id];
                    }
                }
            }

            if (port.getDegree() > 0) {
                portBarycenter[port.id] = sum / port.getDegree();
                minBarycenter = Math.min(minBarycenter, portBarycenter[port.id]);
                maxBarycenter = Math.max(maxBarycenter, portBarycenter[port.id]);
            } else if (northSouthPort) {
                // For northern and southern ports, the sum directly corresponds to the
                // barycenter value to be used.
                portBarycenter[port.id] = sum;
            }
        }
    }

    private void calculateInLayerPortsBarycenterValues(final LNode node) {
        // go through the list of in-layer ports and calculate their barycenter values
        int nodeIndexInLayer = positionOf(node) + 1;
        int layerSize = node.getLayer().getNodes().size() + 1;
        for (LPort inLayerPort : inLayerPorts) {
            // add the indices of all connected in-layer ports
            int sum = 0;
            int inLayerConnections = 0;

            for (LPort connectedPort : inLayerPort.getConnectedPorts()) {
                if (connectedPort.getNode().getLayer() == node.getLayer()) {
                    sum += positionOf(connectedPort.getNode()) + 1;
                    inLayerConnections++;
                }
            }
            // The port's barycenter value is the mean index of connected nodes. If that
            // value is lower than the node's index, most in-layer edges point upwards, so we want
            // the port to be placed at the top of the side. If the value is higher than the
            // nodes's index, most in-layer edges point downwards, so we want the port to be
            // placed at the bottom of the side.
            float barycenter = (float) sum / inLayerConnections;

            PortSide portSide = inLayerPort.getSide();

            if (portSide == PortSide.EAST) {
                if (barycenter < nodeIndexInLayer) {
                    // take a low value in order to have the port above
                    portBarycenter[inLayerPort.id] = minBarycenter - barycenter;
                } else {
                    // take a high value in order to have the port below
                    portBarycenter[inLayerPort.id] = maxBarycenter + (layerSize - barycenter);
                }
            } else if (portSide == PortSide.WEST) {
                if (barycenter < nodeIndexInLayer) {
                    // take a high value in order to have the port above
                    portBarycenter[inLayerPort.id] = maxBarycenter + barycenter;
                } else {
                    // take a low value in order to have the port below
                    portBarycenter[inLayerPort.id] = minBarycenter - (layerSize - barycenter);
                }
            }
        }
    }

    private float dealWithNorthSouthPorts(final float absurdlyLargeFloat, final LPort port,
            final LNode portDummy) {
        // Find out if it's an input port, an output port, or both
        boolean input = false;
        boolean output = false;
        for (LPort portDummyPort : portDummy.getPorts()) {
            if (portDummyPort.getProperty(InternalProperties.ORIGIN) == port) {
                if (!portDummyPort.getOutgoingEdges().isEmpty()) {
                    output = true;
                } else if (!portDummyPort.getIncomingEdges().isEmpty()) {
                    input = true;
                }
            }
        }
        float sum = 0;
        if (input && input ^ output) {
            // It's an input port; the index of its dummy node is its inverted sortkey
            // (for southern input ports, the key must be larger than the ones
            // assigned to output ports or input&&output ports)
            sum = port.getSide() == PortSide.NORTH ? -positionOf(portDummy)
                    : absurdlyLargeFloat - positionOf(portDummy);
        } else if (output && input ^ output) {
            // It's an output port; the index of its dummy node is its sort key
            // (for northern output ports, the key must be larger than the ones assigned
            // to input ports or input&&output ports, which are negative and 0,
            // respectively)
            sum = positionOf(portDummy) + 1.0f;
        } else if (input && output) {
            // It's both, an input and an output port; it must sit between input and
            // output ports
            // North: input ports < 0.0, output ports > 0.0
            // South: input ports > FLOAT_MAX / 2, output ports near zero
            sum = port.getSide() == PortSide.NORTH ? 0.0f : absurdlyLargeFloat / 2f;
        }
        return sum;
    }

    private int positionOf(final LNode node) {
            return nodePositions[node.getLayer().id][node.id];
    }

    private void updateNodePositions(final LNode[][] nodeOrder, final int currentIndex) {
        LNode[] layer = nodeOrder[currentIndex];
        for (int i = 0; i < layer.length; i++) {
            LNode node = layer[i];
            nodePositions[node.getLayer().id][node.id] = i;
        }
    }

    private boolean hasNestedGraph(final LNode node) {
        return node.getNestedGraph() != null;
    }

    private boolean isNotFirstLayer(final int length, final int currentIndex,
            final boolean isForwardSweep) {
        return isForwardSweep ? currentIndex != 0 : currentIndex != length - 1;
    }

    private PortType portTypeFor(final boolean isForwardSweep) {
        return isForwardSweep ? PortType.OUTPUT : PortType.INPUT;
    }

    /**
     * Sort the ports of a node using the given relative position values. These values are
     * interpreted as a hint for the clockwise order of ports.
     *
     * @param node
     *            a node
     */
    private void sortPorts(final LNode node) {
        Collections.sort(node.getPorts(), (port1, port2) -> {
            PortSide side1 = port1.getSide();
            PortSide side2 = port2.getSide();

            if (side1 != side2) {
                // sort according to the node side
                return side1.ordinal() - side2.ordinal();
            } else {
                float port1Bary = portBarycenter[port1.id];
                float port2Bary = portBarycenter[port2.id];
                if (port1Bary == 0 && port2Bary == 0) {
                    return 0;
                } else if (port1Bary == 0) {
                    return -1;
                } else if (port2Bary == 0) {
                    return 1;
                } else {
                    // sort according to the position value
                    return Float.compare(port1Bary, port2Bary);
                }
            }
        });
    }

    private int nPorts;

    @Override
    public void initAtLayerLevel(final int l, final LNode[][] nodeOrder) {
        nodePositions[l] = new int[nodeOrder[l].length];
    }

    @Override
    public void initAtNodeLevel(final int l, final int n, final LNode[][] nodeOrder) {
        LNode node = nodeOrder[l][n];
        node.id = n;
        nodePositions[l][n] = n;
    }

    @Override
    public void initAtPortLevel(final int l, final int n, final int p, final LNode[][] nodeOrder) {
        nodeOrder[l][n].getPorts().get(p).id = nPorts++;
    }

    @Override
    public void initAfterTraversal() {
        portRanks = new float[nPorts];
        portBarycenter = new float[nPorts];
    }
}
