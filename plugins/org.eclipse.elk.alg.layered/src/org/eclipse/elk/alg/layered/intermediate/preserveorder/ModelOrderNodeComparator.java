/*******************************************************************************
 * Copyright (c) 2020, 2023 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 * 
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.alg.layered.intermediate.preserveorder;

import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;

import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LNode.NodeType;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.LongEdgeOrderingStrategy;
import org.eclipse.elk.alg.layered.options.OrderingStrategy;

/**
 * Orders {@link LNode}s in the same layer by {@link InternalProperties#MODEL_ORDER}
 * or the order of the nodes they connect to in the previous layer.
 */
public class ModelOrderNodeComparator implements Comparator<LNode> {
    
    /**
     * The previous layer.
     */
    private LNode[] previousLayer;
    
    /**
     * The ordering strategy.
     */
    private final OrderingStrategy orderingStrategy;
    
    /**
     * Each node has an entry of nodes for which it is bigger.
     */
    private HashMap<LNode, HashSet<LNode>> biggerThan = new HashMap<>();
    /**
     * Each node has an entry of nodes for which it is smaller.
     */
    private HashMap<LNode, HashSet<LNode>> smallerThan = new HashMap<>();
    
    /**
     * Dummy node sorting strategy when compared to nodes with no connection to the previous layer.
     */
    private LongEdgeOrderingStrategy longEdgeNodeOrder = LongEdgeOrderingStrategy.EQUAL;
    
    /**
     * Whether the node comparator was called before ports.
     */
    private boolean beforePorts;
    
    /**
     * Creates a comparator to compare {@link LNode}s in the same layer.
     * 
     * @param thePreviousLayer The previous layer
     * @param orderingStrategy The ordering strategy
     * @param longEdgeOrderingStrategy The strategy to order dummy nodes and nodes with no connection the previous layer
     */
    public ModelOrderNodeComparator(final Layer thePreviousLayer, final OrderingStrategy orderingStrategy,
            final LongEdgeOrderingStrategy longEdgeOrderingStrategy, boolean beforePorts) {
        this(orderingStrategy, longEdgeOrderingStrategy, beforePorts);
        this.previousLayer = new LNode[thePreviousLayer.getNodes().size()];
        thePreviousLayer.getNodes().toArray(this.previousLayer);
    }

    /**
     * Creates a comparator to compare {@link LNode}s in the same layer.
     * 
     * @param previousLayer The previous layer
     * @param orderingStrategy The ordering strategy
     * @param longEdgeOrderingStrategy The strategy to order dummy nodes and nodes with no connection the previous layer
     */
    public ModelOrderNodeComparator(final LNode[] previousLayer, final OrderingStrategy orderingStrategy,
            final LongEdgeOrderingStrategy longEdgeOrderingStrategy, boolean beforePorts) {
        this(orderingStrategy, longEdgeOrderingStrategy, beforePorts);
        this.previousLayer = previousLayer;
    }
    
    private ModelOrderNodeComparator(final OrderingStrategy orderingStrategy,
            final LongEdgeOrderingStrategy longEdgeOrderingStrategy, boolean beforePorts) {
        this.orderingStrategy = orderingStrategy;
        this.longEdgeNodeOrder = longEdgeOrderingStrategy;
        this.beforePorts = beforePorts;
    }

    @Override
    public int compare(final LNode n1, final LNode n2) {
        if (!biggerThan.containsKey(n1)) {
            biggerThan.put(n1, new HashSet<>());
        } else if (biggerThan.get(n1).contains(n2)) {
            return 1;
        }
        if (!biggerThan.containsKey(n2)) {
            biggerThan.put(n2, new HashSet<>());
        } else if (biggerThan.get(n2).contains(n1)) {
            return -1;
        }
        if (!smallerThan.containsKey(n1)) {
            smallerThan.put(n1, new HashSet<>());
        } else if (smallerThan.get(n1).contains(n2)) {
            return -1;
        }
        if (!smallerThan.containsKey(n2)) {
            smallerThan.put(n2, new HashSet<>());
        } else if (biggerThan.get(n2).contains(n1)) {
            return 1;
        }
        // If no model order is set, the one node is a dummy node and the nodes should be ordered
        // by the connected edges.
        // This kind of ordering should be preferred, if the order of the edges has priority.
        if (orderingStrategy == OrderingStrategy.PREFER_EDGES || !n1.hasProperty(InternalProperties.MODEL_ORDER)
                || !n2.hasProperty(InternalProperties.MODEL_ORDER)) {
            // In this case the order of the connected nodes in the previous layer should be respected
            LPort p1SourcePort = null;
            for (LPort p : n1.getPorts()) {
                // Get the first port that actually connects to a previous layer.
                if (!p.getIncomingEdges().isEmpty()) {
                    if (p.getIncomingEdges().get(0).getSource().getNode().getLayer().id == (n1.getLayer().id - 1)) {
                        p1SourcePort = p.getIncomingEdges().get(0).getSource();
                    }
                }
            }
            
            LPort p2SourcePort = null;
            for (LPort p : n2.getPorts()) {
                // Get the first port that actually connects to a previous layer.
                if (!p.getIncomingEdges().isEmpty()) {
                    if (p.getIncomingEdges().get(0).getSource().getNode().getLayer().id == (n2.getLayer().id - 1)) {
                        p2SourcePort = p.getIncomingEdges().get(0).getSource();
                    }
                }
            }
            
            // Case both nodes have connections to the previous layer.
            if (p1SourcePort != null && p2SourcePort != null) {
                LNode p1Node = p1SourcePort.getNode();
                LNode p2Node = p2SourcePort.getNode();
                
                // If both nodes connect to the same node the order of their corresponding ports in the previous
                // layer should be used to order them.
                if (p1Node != null && p1Node.equals(p2Node)) {
                    // We are not allowed to look at the model order of the edges but we have to look at the actual
                    // port ordering since the edge order might be broken here.
                    // If we have long edges the edge model order might not respect the ordering.
                    for (LPort port : p1Node.getPorts()) {
                        if (port.equals(p1SourcePort)) {
                            // Case the port is the one connecting to n1, therefore, n1 has a smaller model order
                            updateBiggerAndSmallerAssociations(n2, n1);
                            return -1;
                        } else if (port.equals(p2SourcePort)) {
                            // Case the port is the one connecting to n2, therefore, n1 has a bigger model order
                            updateBiggerAndSmallerAssociations(n1, n2);
                            return 1;
                        }
                    }
                    assert (false);
                    // Cannot happen, since both nodes have a connection to the previous layer.
                    // Nevertheless, if it does happen provide well defined behavior and use the model order of the edges.
                    int n1EdgeOrder = getModelOrderFromConnectedEdges(n1);
                    int n2EdgeOrder = getModelOrderFromConnectedEdges(n2);
                    if (n1EdgeOrder > n2EdgeOrder) {
                        updateBiggerAndSmallerAssociations(n1, n2);
                        return 1;
                    } else {
                        // I assume that equal is not an alternative here.
                        updateBiggerAndSmallerAssociations(n2, n1);
                        return -1;
                    }
                }
                
                // If the nodes do not connect to the same node, the nodes are ordered by the nodes they connect to.
                // One can disregard the model order here 
                // since the ordering in the previous layer does already reflect it.
                for (LNode previousNode : previousLayer) {
                    if (previousNode.equals(p1Node)) {
                        updateBiggerAndSmallerAssociations(n2, n1);
                        return -1;
                    } else if (previousNode.equals(p2Node)) {
                        updateBiggerAndSmallerAssociations(n1, n2);
                        return 1;
                    }
                }
                // Again, I assume that such a node must exist, and I should never get here.
            }
            
            // One node has no source port
            if (p1SourcePort != null && p2SourcePort == null || p1SourcePort == null && p2SourcePort != null) {
                // Check whether one of them is a helper dummy node.
                int comparedWithLongEdgeFeedback = handleHelperDummyNodes(n1, n2);
                if (comparedWithLongEdgeFeedback != 0) {
                    if (comparedWithLongEdgeFeedback > 0) {
                        updateBiggerAndSmallerAssociations(n1, n2);
                    } else {
                        updateBiggerAndSmallerAssociations(n2, n1);
                    }
                    return comparedWithLongEdgeFeedback;
                }
                
                // Otherwise use the model order of the connected edges if one of them is no dummy node.
                if (!n1.hasProperty(InternalProperties.MODEL_ORDER) || !n2.hasProperty(InternalProperties.MODEL_ORDER)) {
                    int n1ModelOrder = getModelOrderFromConnectedEdges(n1);
                    int n2ModelOrder = getModelOrderFromConnectedEdges(n2);
                    if (n1ModelOrder > n2ModelOrder) {
                        updateBiggerAndSmallerAssociations(n1, n2);
                        return 1;
                    } else {
                        updateBiggerAndSmallerAssociations(n2, n1);
                        return -1;
                    }
                }
            }

            // Both nodes are not connected to the previous layer
            if (p1SourcePort == null && p2SourcePort == null) {
                // Check whether one of them is a helper dummy node.
                int comparedWithLongEdgeFeedback = handleHelperDummyNodes(n1, n2);
                if (comparedWithLongEdgeFeedback != 0) {
                    if (comparedWithLongEdgeFeedback > 0) {
                        updateBiggerAndSmallerAssociations(n1, n2);
                    } else {
                        updateBiggerAndSmallerAssociations(n2, n1);
                    }
                    return comparedWithLongEdgeFeedback;
                }
            }
            
            // Fall through case.
            // Both nodes are not connected to the previous layer and are normal nodes.
            // The node model order shall be used to order them.
        }
        // Order nodes by their order in the model.
        // This is also the fallback case if one of the nodes is not connected to the previous layer.
        if (n1.hasProperty(InternalProperties.MODEL_ORDER) && n2.hasProperty(InternalProperties.MODEL_ORDER)) {
            int n1ModelOrder = n1.getProperty(InternalProperties.MODEL_ORDER);
            int n2ModelOrder = n2.getProperty(InternalProperties.MODEL_ORDER);
            if (n1ModelOrder > n2ModelOrder) {
                updateBiggerAndSmallerAssociations(n1, n2);
                return 1;
            } else {
                updateBiggerAndSmallerAssociations(n2, n1);
                return -1;
            }
        } else {
            // If they have no model order, I should still make an ordering decision that I save somehow.
            // This decision is somehow random I just sort the first one before the second one.
            updateBiggerAndSmallerAssociations(n2, n1);
            return -1;
        }
    }
    
    /**
     * The {@link InternalProperties#MODEL_ORDER} of the first incoming edge of a node.
     * 
     * @param n The node
     * @return The model order of the first incoming edge of the given node.
     * Returns Integer.MAX_VALUE if no such edge exists.
     */
    private int getModelOrderFromConnectedEdges(final LNode n) {
        LPort sourcePort = n.getPorts().stream().filter(p -> !p.getIncomingEdges().isEmpty()).findFirst().orElse(null);
        if (sourcePort != null) {
            LEdge edge = sourcePort.getIncomingEdges().get(0);
            if (edge != null) {
                return edge.getProperty(InternalProperties.MODEL_ORDER);
            }
        }
        // Set to MAX_INT to sort dummy nodes over nodes without a connection to the previous layer.
        // Set to -MAX_INT to sort dummy nodes under nodes without a connection to the previous layer.
        // Set to 0 if you do not care about their order.
        // One of this has to be chosen, since dummy nodes are not comparable with nodes
        // that do not have a connection to the previous layer.
        return longEdgeNodeOrder.returnValue();
    }
    
    private void updateBiggerAndSmallerAssociations(final LNode bigger, final LNode smaller) {
        HashSet<LNode> biggerNodeBiggerThan = biggerThan.get(bigger);
        HashSet<LNode> smallerNodeBiggerThan = biggerThan.get(smaller);
        HashSet<LNode> biggerNodeSmallerThan = smallerThan.get(bigger);
        HashSet<LNode> smallerNodeSmallerThan = smallerThan.get(smaller);
        biggerNodeBiggerThan.add(smaller);
        smallerNodeSmallerThan.add(bigger);
        for (LNode verySmall : smallerNodeBiggerThan) {
            biggerNodeBiggerThan.add(verySmall);
            smallerThan.get(verySmall).add(bigger);
            smallerThan.get(verySmall).addAll(biggerNodeSmallerThan);
        }
        

        for (LNode veryBig : biggerNodeSmallerThan) {
            smallerNodeSmallerThan.add(veryBig);
            biggerThan.get(veryBig).add(smaller);
            biggerThan.get(veryBig).addAll(smallerNodeBiggerThan);
        }
    }
    
    private int handleHelperDummyNodes(LNode n1, LNode n2) {
        if (n1.getType() == NodeType.LONG_EDGE && n2.getType() == NodeType.NORMAL) {
            // n1 could be a long edge node feedback node.
            
            LPort dummyNodeSourcePort = getFirstIncomingSourcePortOfNode(n1);
            LNode dummyNodeSourceNode = dummyNodeSourcePort.getNode();
            LPort dummyNodeTargetPort = getFirstOutgoingTargetPortOfNode(n1);
            LNode dummyNodeTargetNode = dummyNodeTargetPort.getNode();
            int dummyLayerId = n1.getLayer().id;
            
            // Check whether the dummy node is feedback source or feedback target if not return.
            if (dummyNodeSourceNode.getLayer().id != dummyLayerId && dummyNodeTargetNode.getLayer().id != dummyLayerId) {
                return 0;
            }
            // Case the source of the dummy is the same node as n2, than the dummy node is routed below.
            if (dummyNodeSourceNode.equals(n2)) {
                updateBiggerAndSmallerAssociations(n1, n2);
                return 1;
            } else {
                // Calculate whether the dummy node leads to the target node.
                if (dummyNodeTargetNode.equals(n2)) {
                    updateBiggerAndSmallerAssociations(n1, n2);
                    return 1;
                }
                
                // If the two nodes are not the same, order them based on the source model order.
                return this.compare(dummyNodeSourceNode, n2);
            }
        } else if (n1.getType() == NodeType.NORMAL && n2.getType() == NodeType.LONG_EDGE) {
            // n2 could be a long edge node feedback node.
            LPort dummyNodeSourcePort = getFirstIncomingSourcePortOfNode(n2);
            LNode dummyNodeSourceNode = dummyNodeSourcePort.getNode();
            LPort dummyNodeTargetPort = getFirstOutgoingTargetPortOfNode(n2);
            LNode dummyNodeTargetNode = dummyNodeTargetPort.getNode();
            int dummyLayerId = n1.getLayer().id;
            
            // Check whether the dummy node is feedback source or feedback target if not return.
            if (dummyNodeSourceNode.getLayer().id != dummyLayerId && dummyNodeTargetNode.getLayer().id != dummyLayerId) {
                return 0;
            }
            // Case the source of the dummy is the same node as n2, than the dummy node is routed below.
            if (dummyNodeSourceNode.equals(n1)) {
                updateBiggerAndSmallerAssociations(n2, n1);
                return -1;
            } else {
                // Calculate whether the dummy node leads to the target node.
                if (dummyNodeTargetNode.equals(n1)) {
                    updateBiggerAndSmallerAssociations(n2, n1);
                    return -1;
                }
                
                // If the two nodes are not the same, order them based on the source model order.
                return this.compare(n1, dummyNodeSourceNode);
            }
        } else if (n1.getType() == NodeType.LONG_EDGE && n2.getType() == NodeType.LONG_EDGE) {
            // One of these edges is a feedback edge. This has to be the case since at least one of them is not connected
            // to the previous layer.
            // If only one is a long edge feedback node, I have a problem since these are not comparable.
            // I must find the reference node in the current layer of each edge and use it instead.
            LPort n1dummyNodeSourcePort = getFirstIncomingSourcePortOfNode(n1);
            LPort n1dummyNodeTargetPort = getFirstOutgoingTargetPortOfNode(n1);
            LNode n1dummySourceNode = n1dummyNodeSourcePort.getNode();
            LNode n1dummyTargetNode = n1dummyNodeTargetPort.getNode();
            int n1LayerId = n1.getLayer().id;
            boolean n1SourceFeedbackNode = false;
            boolean n1TargetFeedbackNode = false;
            
            LPort n2dummyNodeSourcePort = getFirstIncomingSourcePortOfNode(n2);
            LPort n2dummyNodeTargetPort = getFirstOutgoingTargetPortOfNode(n2);
            LNode n2dummySourceNode = n2dummyNodeSourcePort.getNode();
            LNode n2dummyTargetNode = n2dummyNodeTargetPort.getNode();
            int n2LayerId = n2.getLayer().id;
            boolean n2SourceFeedbackNode = false;
            boolean n2TargetFeedbackNode = false;
            
            LNode n1ReferenceNode = n1;
            LNode n2ReferenceNode = n2;
            if (n1dummySourceNode.getLayer().id == n1LayerId) {
                // This means that n1dummySourceNode is the reference node that we need to consider for ordering n1;
                n1SourceFeedbackNode = true;
                n1ReferenceNode = n1dummySourceNode;
            } else if (n1dummyTargetNode.getLayer().id == n1LayerId) {
                // This means that n1dummyNodeTargetPort is the reference node that we need to consider for ordering n1;
                n1TargetFeedbackNode = true;
                n1ReferenceNode = n1dummyTargetNode;
            }
            if (n2dummySourceNode.getLayer().id == n2LayerId) {
                // This means that n2dummySourceNode is the reference node that we need to consider for ordering n2;
                n2SourceFeedbackNode = true;
                n2ReferenceNode = n2dummySourceNode;
            } else if (n2dummyTargetNode.getLayer().id == n2LayerId) {
                // This means that n2dummyNodeTargetPort is the reference node that we need to consider for ordering n2;
                n2TargetFeedbackNode = true;
                n2ReferenceNode = n2dummyTargetNode;
            }
            
            // After this each reference node should be a real node in this layer that I can use to compare n1 and n2.
            
            // Case both are on the same node.
            if (n1ReferenceNode.equals(n2ReferenceNode)) {
                // Find the first port that occurs on the node.
                // Since we have a feedback node, we need to reverse the decision.
                // If the side we connect to is WEST. we also have to reverse the decision.
                // This may be a problem since here the node comparator depends on the port order in the same layer.
                if (this.beforePorts) {
                    // The order on the reference node ports may just be wrong since the ports are not yet sorted.
                    // If both reference nodes are source feedback nodes, then the model order (considering reversing) does the trick.
                    // If one is source one is target, I will just order them but it will always create a crossing.
                    // If both are target, I should not have this problem and can never be here, since these nodes
                    // should have a previous layer node.
                    if (n1SourceFeedbackNode && n2SourceFeedbackNode) {
                        int returnValue = new ModelOrderPortComparator(previousLayer, orderingStrategy, null, n2TargetFeedbackNode)
                            .compare(n1dummyNodeSourcePort, n2dummyNodeSourcePort);
                        if (returnValue > 0) {
                            updateBiggerAndSmallerAssociations(n2, n1);
                            return 1;
                        } else {
                            updateBiggerAndSmallerAssociations(n1, n2);
                            return -1;
                        }
                    } else if (n1SourceFeedbackNode && n2TargetFeedbackNode) {
                        updateBiggerAndSmallerAssociations(n2, n1);
                        return 1;
                    } else if (n1TargetFeedbackNode && n2SourceFeedbackNode) {
                        updateBiggerAndSmallerAssociations(n1, n2);
                        return -1;
                    } else if (n1TargetFeedbackNode && n2TargetFeedbackNode) {
                        // In this case, there must be incoming edges that can be used for ordering.
                        return 0;
                    }
                } else {
                    // In this case, the order of the ports can just be used.
                    // Since the order of WEST ports is reversed and the order for feedback edges connecting to WEST
                    // ports is reversed, I may just do nothing here.
                    for (LPort port : n1ReferenceNode.getPorts()) {
                        if (n1dummyNodeSourcePort.equals(port)) {
                            updateBiggerAndSmallerAssociations(n2, n1);
                            return -1;
                        } else if (n2dummyNodeSourcePort.equals(port)) {
                            updateBiggerAndSmallerAssociations(n1, n2);
                            return 1;
                        }
                    }
                }
            }
            
            // If the nodes are different, just compare them since one should be a normal non-feedback dummy now.
            // Worst case would be that one is a dummy and one a normal dangling node, where this would just create a 
            // static ordering.
            return this.compare(n1ReferenceNode, n2ReferenceNode);
        } else {
            // These nodes are just two normal nodes and need to be handled by node model order.
            return 0;
        }
    }
    
    /**
     * Helper method to get the first incoming port of a node.
     * 
     * @param node The node
     * @return The first incoming port.
     */
    private LPort getFirstIncomingPortOfNode(LNode node) {
        return node.getPorts().stream().filter(p -> !p.getIncomingEdges().isEmpty()).findFirst().orElse(null);
    }

    
    /**
     * Helper method to get the source port of the first incoming port of a node.
     * 
     * @param node The node
     * @return The source port of the first incoming port.
     */
    private LPort getFirstIncomingSourcePortOfNode(LNode node) {
        return getFirstIncomingPortOfNode(node).getIncomingEdges().get(0).getSource();
    }

    /**
     * Helper method to get the first outgoing port of a node.
     * 
     * @param node The node
     * @return The first outgoing port.
     */
    private LPort getFirstOutgoingPortOfNode(LNode node) {
        return node.getPorts().stream().filter(p -> !p.getOutgoingEdges().isEmpty()).findFirst().orElse(null);
    }

    /**
     * Helper method to get the target port of the first outgoing port of a node.
     * 
     * @param node The node
     * @return The target port of the first outgoing port.
     */
    private LPort getFirstOutgoingTargetPortOfNode(LNode node) {
        return getFirstOutgoingPortOfNode(node).getOutgoingEdges().get(0).getTarget();
    }
    
    /**
     * Clears the transitive ordering.
     */
    public void clearTransitiveOrdering() {
        this.biggerThan = new HashMap<>();
        this.smallerThan = new HashMap<>();
    }
}
