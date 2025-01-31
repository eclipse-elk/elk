/*******************************************************************************
 * Copyright (c) 2020 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 * 
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.alg.layered.intermediate.preserveorder;

import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.stream.Collectors;

import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.graph.LNode.NodeType;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.OrderingStrategy;
import org.eclipse.elk.core.options.PortSide;

/**
 * Orders {@link LPort}s in the same layer by their edges {@link InternalProperties#MODEL_ORDER}
 * or the order of the nodes they connect to in the previous layer.
 * This takes into account that ports that connect to the same (long edge) target should be ordered next to each other.
 * Outgoing ports are ordered before incoming ports.
 * Incoming ports choose a position that does not create conflicts with the previous layer.
 * FIXME optimize this by setting ids on all already handled nodes in the previous and current layer (e.g. after the 
 * node order was set) to avoid searching for the correct order and using the ids instead to compare.
 */
public class ModelOrderPortComparator implements Comparator<LPort> {

    /**
     * The model order associated to a target node of this node.
     */
    private final Map<LNode, Integer> targetNodeModelOrder;
    
    /**
     * Whether ports should have a model order.
     */
    private boolean portModelOrder;

    /**
     * The previous layer.
     */
    private LNode[] previousLayer;
    
    private OrderingStrategy strategy;
    
    /**
     * Each port has an entry of ports for which it is bigger.
     */
    private HashMap<LPort, HashSet<LPort>> biggerThan = new HashMap<>();
    /**
     * Each port has an entry of ports for which it is smaller.
     */
    private HashMap<LPort, HashSet<LPort>> smallerThan = new HashMap<>();

    /**
     * Creates a comparator to compare {@link LPort}s in the same layer.
     * 
     * @param previousLayer The previous layer
     * @param targetNodeModelOrder The minimal model order connecting to a target node.
     */
    public ModelOrderPortComparator(final Layer previousLayer, final OrderingStrategy strategy,
            final Map<LNode, Integer> targetNodeModelOrder, final boolean portModelOrder) {
        this.previousLayer = new LNode[previousLayer.getNodes().size()];
        this.strategy = strategy;
        previousLayer.getNodes().toArray(this.previousLayer);
        this.targetNodeModelOrder = targetNodeModelOrder;
        this.portModelOrder = portModelOrder;
    }

    /**
     * Creates a comparator to compare {@link LPort}s in the same layer.
     * 
     * @param previousLayer The previous layer
     * @param targetNodeModelOrder The minimal model order connecting to a target node.
     */
    public ModelOrderPortComparator(final LNode[] previousLayer, final OrderingStrategy strategy,
            final Map<LNode, Integer> targetNodeModelOrder, final boolean portModelOrder) {
        this.previousLayer = previousLayer;
        this.strategy = strategy;
        this.targetNodeModelOrder = targetNodeModelOrder;
        this.portModelOrder = portModelOrder;
    }

    @Override
    public int compare(final LPort originalP1, final LPort originalP2) {
        LPort p1 = originalP1;
        LPort p2 = originalP2;

        if (!biggerThan.containsKey(p1)) {
            biggerThan.put(p1, new HashSet<>());
        } else if (biggerThan.get(p1).contains(p2)) {
            return 1;
        }
        if (!biggerThan.containsKey(p2)) {
            biggerThan.put(p2, new HashSet<>());
        } else if (biggerThan.get(p2).contains(p1)) {
            return -1;
        }
        if (!smallerThan.containsKey(p1)) {
            smallerThan.put(p1, new HashSet<>());
        } else if (smallerThan.get(p1).contains(p2)) {
            return -1;
        }
        if (!smallerThan.containsKey(p2)) {
            smallerThan.put(p2, new HashSet<>());
        } else if (biggerThan.get(p2).contains(p1)) {
            return 1;
        }
        // Sort nodes by their ports sides NORTH < EAST < SOUTH < WEST.
        if (p1.getSide() != p2.getSide()) {
            int result =  new PortSideComparator().compare(p1.getSide(), p2.getSide());
            if (result > 0) {
                updateBiggerAndSmallerAssociations(p1, p2, 1);
            } else {
                updateBiggerAndSmallerAssociations(p2, p1, 1);
            }
            return result;
        }
        int reverseOrder = 1;
        
        // Sort incoming edges by sorting their ports by the order of the nodes they connect to.
        if (!p1.getIncomingEdges().isEmpty() && !p2.getIncomingEdges().isEmpty()) {
            if (p1.getSide() == PortSide.WEST && p2.getSide() == PortSide.WEST
                    || p1.getSide() == PortSide.NORTH && p2.getSide() == PortSide.NORTH
                    || p1.getSide() == PortSide.SOUTH && p2.getSide() == PortSide.SOUTH) {
                // Some ports are ordered in the way around.
                // Previously this did not matter, since the north south processor did override the ordering.
                reverseOrder = -reverseOrder;
            }

            LPort p1SourcePort = p1.getIncomingEdges().get(0).getSource();
            LPort p2SourcePort = p2.getIncomingEdges().get(0).getSource();
            LNode p1Node = p1SourcePort.getNode();
            LNode p2Node = p2SourcePort.getNode();
            if (p1Node.equals(p2Node)) {
                // If both connect to the same node, check their occurrence order since these ports have to be handled
                // first and should hence be sorted.
                for (LPort port : p1Node.getPorts()) {
                    if (p1SourcePort.equals(port)) {
                        updateBiggerAndSmallerAssociations(p2, p1, reverseOrder);
                        return -reverseOrder;
                    } else if (p2SourcePort.equals(port)) {
                        updateBiggerAndSmallerAssociations(p1, p2, reverseOrder);
                        return reverseOrder;
                    }
                }
            }
            // If both ports connect to long edges in the same layer, reverse the order.
            if (p1SourcePort.getNode().getType() == NodeType.LONG_EDGE
                && p2SourcePort.getNode().getType() == NodeType.LONG_EDGE
                && p1Node.getLayer().id == p2Node.getLayer().id && p1Node.getLayer().id == p1.getNode().getLayer().id) {
                //       _
                // n1  n2_|
                // ||    ||
                // |__2__||
                // ___1___|
                //
                Layer previousLayer = p1Node.getLayer();
                // if it is a SOUTH or EAST port reverse the order.
                // checkReferenceLayer already updates the transitive ordering association.
                int inPreviousLayer =  checkReferenceLayer(previousLayer, p1Node, p2Node, p1, p2);
                if (inPreviousLayer != 0) {
                    if (p1.getSide() == PortSide.EAST && p2.getSide() == PortSide.EAST) {
                        // Some ports are ordered in the way around.
                        // Previously this did not matter, since the north south processor did override the ordering.
                        reverseOrder = -reverseOrder;
                    }
                    if (inPreviousLayer > 0) {
                        updateBiggerAndSmallerAssociations(p1, p2, reverseOrder);
                        return reverseOrder;
                    } else {
                        updateBiggerAndSmallerAssociations(p2, p1, reverseOrder);
                        return -reverseOrder;
                    }
                }
            }
            
            // Check which of the nodes connects first to the previous layer.
            // checkReferenceLayer already updates the transitive ordering association.
            int inPreviousLayer =  checkReferenceLayer(Arrays.stream(previousLayer).collect(Collectors.toList()), p1Node, p2Node, p1, p2);
            if (inPreviousLayer != 0) {
                if (inPreviousLayer > 0) {
                    updateBiggerAndSmallerAssociations(p1, p2, reverseOrder);
                    return reverseOrder;
                } else {
                    updateBiggerAndSmallerAssociations(p2, p1, reverseOrder);
                    return -reverseOrder;
                }
            }
            // If both ports do not connect to the previous layer, use the port model order.
            // If port order is used instead of edge order, consult it to make decisions.
            // If not both ports have a model order fall back to the edge order.
            if (portModelOrder) {
                int result =  checkPortModelOrder(p1, p2);
                if (result != 0) {
                    if (result > 0) {
                        updateBiggerAndSmallerAssociations(p1, p2, reverseOrder);
                        return reverseOrder;
                    } else {
                        updateBiggerAndSmallerAssociations(p2, p1, reverseOrder);
                        return -reverseOrder;
                    }
                }
            }
        }

        // Sort outgoing edges by sorting their ports based on the model order of their edges.
        if (!p1.getOutgoingEdges().isEmpty() && !p2.getOutgoingEdges().isEmpty()) {
            if (p1.getSide() == PortSide.WEST && p2.getSide() == PortSide.WEST
                    || p1.getSide() == PortSide.SOUTH && p2.getSide() == PortSide.SOUTH) {
                // Some ports are ordered in the way around.
                // Previously this did not matter, since the north south processor did override the ordering.
                reverseOrder = -reverseOrder;
            }
            LNode p1TargetNode = p1.getProperty(InternalProperties.LONG_EDGE_TARGET_NODE);
            LNode p2TargetNode = p2.getProperty(InternalProperties.LONG_EDGE_TARGET_NODE);

            // If node order is preferred and both edges have a target node use the model order of the long edge target.
            if (this.strategy == OrderingStrategy.PREFER_NODES && p1TargetNode != null && p2TargetNode != null
                    && p1TargetNode.hasProperty(InternalProperties.MODEL_ORDER)
                    && p2TargetNode.hasProperty(InternalProperties.MODEL_ORDER)) {
                int p1MO = p1TargetNode.getProperty(InternalProperties.MODEL_ORDER);
                int p2MO = p2TargetNode.getProperty(InternalProperties.MODEL_ORDER);
                if (p1MO > p2MO) {
                    updateBiggerAndSmallerAssociations(p1, p2, reverseOrder);
                    return reverseOrder;
                } else {
                    updateBiggerAndSmallerAssociations(p2, p1, reverseOrder);
                    return -reverseOrder;
                }
            }

            // If port order is used instead of edge order, consult it to make decisions.
            // If not both ports have a model order fall back to the edge order.
            if (portModelOrder) {
                int result =  checkPortModelOrder(p1, p2);
                if (result != 0) {
                    if (result > 0) {
                        updateBiggerAndSmallerAssociations(p1, p2, reverseOrder);
                        return reverseOrder;
                    } else {
                        updateBiggerAndSmallerAssociations(p2, p1, reverseOrder);
                        return -reverseOrder;
                    }
                }
                // If one of the ports has no model order, find something else to compare them which is the edge order.
            }
            
            int p1Order = 0;
            int p2Order = 0;
            if (p1.getOutgoingEdges().get(0).hasProperty(InternalProperties.MODEL_ORDER)) {
                p1Order = p1.getOutgoingEdges().get(0).getProperty(InternalProperties.MODEL_ORDER);
            }
            if (p2.getOutgoingEdges().get(0).hasProperty(InternalProperties.MODEL_ORDER)) {
                p2Order = p2.getOutgoingEdges().get(0).getProperty(InternalProperties.MODEL_ORDER);
            }
            
            // If both ports have the same target nodes, make sure that the backward edge is below the normal edge.
            if (p1TargetNode != null && p1TargetNode.equals(p2TargetNode)) {
                // If both edges are reversed or not reversed, just use their model order.
                if (p1Order > p2Order) {
                    updateBiggerAndSmallerAssociations(p1, p2, reverseOrder);
                    return reverseOrder;
                } else {
                    updateBiggerAndSmallerAssociations(p2, p1, reverseOrder);
                    return -reverseOrder;
                }
            }
            // Use precomputed ordering value if possible to utilize order inheritence of edges connected to a node.
            // This allows to bundle edges leading to the same node, disregarding their model order.
            if (targetNodeModelOrder != null) {
                if (targetNodeModelOrder.containsKey(p1TargetNode)) {
                    p1Order = targetNodeModelOrder.get(p1TargetNode);
                }
                if (targetNodeModelOrder.containsKey(p2TargetNode)) {
                    p2Order = targetNodeModelOrder.get(p2TargetNode);
                }
            }
            // If the nodes have different targets just use their order.
            if (p1Order > p2Order) {
                updateBiggerAndSmallerAssociations(p1, p2, reverseOrder);
                return reverseOrder;
            } else {
                updateBiggerAndSmallerAssociations(p2, p1, reverseOrder);
                return -reverseOrder;
            }
            
        }
        // Sort outgoing ports before incoming ports.
        if (!p1.getIncomingEdges().isEmpty() && !p2.getOutgoingEdges().isEmpty()) {
            updateBiggerAndSmallerAssociations(p1, p2, reverseOrder);
            return 1;
        } else if (!p1.getOutgoingEdges().isEmpty() && !p2.getIncomingEdges().isEmpty()) {
            updateBiggerAndSmallerAssociations(p2, p1, reverseOrder);
            return -1;
        } else if (p1.hasProperty(InternalProperties.MODEL_ORDER) && p2.hasProperty(InternalProperties.MODEL_ORDER)) {
            // The ports have no edges.
            // Use the port model order to compare them. This can always be very bad since these unconnected ports
            // can transitively order nodes that should be ordered differently.
            // This can only be prevented, if one handles the sorting such that unconnected ports are handled last.
            int p1MO = p1.getProperty(InternalProperties.MODEL_ORDER);
            int p2MO = p2.getProperty(InternalProperties.MODEL_ORDER);
            // Still check the side, since WEST and SOUTH must be the other way around.
            if (p1.getSide() == PortSide.WEST && p2.getSide() == PortSide.WEST
                    || p1.getSide() == PortSide.SOUTH && p2.getSide() == PortSide.SOUTH) {
                // Some ports are ordered in the way around.
                // Previously this did not matter, since the north south processor did override the ordering.
                reverseOrder = -reverseOrder;
            }
            
            if (p1MO > p2MO) {
                updateBiggerAndSmallerAssociations(p1, p2, reverseOrder);
                return reverseOrder;
            } else {
                updateBiggerAndSmallerAssociations(p2, p1, reverseOrder);
                return -reverseOrder;
            }
        } else {
            updateBiggerAndSmallerAssociations(p2, p1, reverseOrder);
            return -reverseOrder;
        }
    }
    
    /**
     * Compares two ports by considering the model order of them.
     * 
     * @param p1 The first port.
     * @param p2 The second port.
     * @return a negative integer, zero, or a positive integer as the
     *         first argument is less than, equal to, or greater than the
     *         second.
     */
    public int checkPortModelOrder(final LPort p1, final LPort p2) {
        if (p1.hasProperty(InternalProperties.MODEL_ORDER)
                && p2.hasProperty(InternalProperties.MODEL_ORDER)) {
            return Integer.compare(p1.getProperty(InternalProperties.MODEL_ORDER),
                    p2.getProperty(InternalProperties.MODEL_ORDER));
        }
        return 0;
    }
    
    private void updateBiggerAndSmallerAssociations(final LPort biggerOri, final LPort smallerOri, int reverseOrder) {
        LPort bigger = biggerOri;
        LPort smaller = smallerOri;
        if (reverseOrder < 0) {
            bigger = smallerOri;
            smaller = biggerOri;
        }
        HashSet<LPort> biggerPortBiggerThan = biggerThan.get(bigger);
        HashSet<LPort> smallerPortBiggerThan = biggerThan.get(smaller);
        HashSet<LPort> biggerPortSmallerThan = smallerThan.get(bigger);
        HashSet<LPort> smallerPortSmallerThan = smallerThan.get(smaller);
        biggerPortBiggerThan.add(smaller);
        smallerPortSmallerThan.add(bigger);
        for (LPort verySmall : smallerPortBiggerThan) {
            biggerPortBiggerThan.add(verySmall);
            smallerThan.get(verySmall).add(bigger);
            smallerThan.get(verySmall).addAll(biggerPortSmallerThan);
        }
        

        for (LPort veryBig : biggerPortSmallerThan) {
            smallerPortSmallerThan.add(veryBig);
            biggerThan.get(veryBig).add(smaller);
            biggerThan.get(veryBig).addAll(smallerPortBiggerThan);
        }
    }
    
    /**
     * Given a previous layer, check which of the two reference nodes of a port is the first in it.
     * 
     * @param layer The layer to check
     * @param p1Node The reference node of port p1
     * @param p2Node The reference node of port p2
     * @param p1 The first port
     * @param p2 The second port
     * @return A comparator value showing which port should be first.
     */
    private int checkReferenceLayer(Iterable<LNode> layer, LNode p1Node, LNode p2Node, LPort p1, LPort p2) {
        for (LNode node : layer) {
            if (node.equals(p1Node)) {
                // If the first node is found first, it has to be above the second one and does hence have a smaller
                // Model order.
                return -1;
            } else if (node.equals(p2Node)) {
                return 1;
            }
        }
        return 0; // Should never happen, the previous layer needs these nodes in it.
    }
    
    private class PortSideComparator implements Comparator<PortSide> {

        /* (non-Javadoc)
         * @see java.util.Comparator#compare(java.lang.Object, java.lang.Object)
         */
        @Override
        public int compare(final PortSide ps1, final PortSide ps2) {
            return Integer.compare(ps1.ordinal(), ps2.ordinal());
        }
        
    }
    
    /**
     * Clears the transitive ordering.
     */
    public void clearTransitiveOrdering() {
        this.biggerThan = new HashMap<>();
        this.smallerThan = new HashMap<>();
    }
}
