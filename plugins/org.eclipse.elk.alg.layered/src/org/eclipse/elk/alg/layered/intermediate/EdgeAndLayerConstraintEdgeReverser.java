/*******************************************************************************
 * Copyright (c) 2010, 2017 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.layered.intermediate;

import java.util.ArrayList;
import java.util.List;

import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LGraphUtil;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LNode.NodeType;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.options.EdgeConstraint;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.LayerConstraint;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.layered.options.PortType;
import org.eclipse.elk.core.alg.ILayoutProcessor;
import org.eclipse.elk.core.options.PortSide;
import org.eclipse.elk.core.util.IElkProgressMonitor;

/**
 * Makes sure nodes with edge or layer constraints have only incoming or only outgoing edges, as appropriate. This is
 * done even before cycle breaking because the result may already break some cycles (or, indeed, cause new ones). This
 * processor is required for {@link LayerConstraintProcessor} to work correctly. If edge constraints are in conflict
 * with layer constraints, the latter take precedence. Furthermore, this processor handles nodes with fixed port sides
 * for which all ports are reversed, i.e. input ports are on the right and output ports are on the left. All incident
 * edges are reversed in such cases.
 * 
 * <p>
 * Special handling applies to nodes that are to be placed in the {@code FIRST} or {@code LAST} layer if they have
 * incoming or outgoing edges with labels, respectively. The labels are represented by label dummy nodes, which will
 * later be placed in a separate layer between {@code FIRST_SEPARATE} and {@code FIRST} (or {@code LAST} and
 * {@code LAST_SEPARATE}).
 * </p>
 * 
 * <p>
 * For this processor, <em>outer nodes</em> are nodes with their layer constraint set to
 * {@link LayerConstraint#FIRST_SEPARATE}, {@link LayerConstraint#FIRST}, {@link LayerConstraint#LAST} or
 * {@link LayerConstraint#LAST_SEPARATE}. Nodes that are not outer nodes are called <em>inner nodes</em>.
 * 
 * <dl>
 *   <dt>Precondition:</dt>
 *     <dd>an unlayered graph.</dd>
 *   <dt>Postcondition:</dt>
 *     <dd>nodes with layer constraints have only incoming or only outgoing edges, as appropriate.</dd>
 *   <dt>Slots:</dt>
 *     <dd>Before phase 1.</dd>
 *   <dt>Same-slot dependencies:</dt>
 *     <dd>None.</dd>
 * </dl>
 * 
 * @see LayerConstraintProcessor
 */
public final class EdgeAndLayerConstraintEdgeReverser implements ILayoutProcessor<LGraph> {

    @Override
    public void process(final LGraph layeredGraph, final IElkProgressMonitor monitor) {
        monitor.begin("Edge and layer constraint edge reversal", 1);

        // Apply edge reversal to all FIRST and FIRST_SEPARATE as well as LAST and LAST_SEPARATE nodes
        List<LNode> remainingNodes = handleOuterNodes(layeredGraph);

        // Apply edge reversal to remaining nodes
        handleInnerNodes(layeredGraph, remainingNodes);

        monitor.done();
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Outer Nodes

    /**
     * Ensures that edges incident to outer nodes are correct.
     * 
     * @return list of those nodes that are not outer nodes.
     */
    private List<LNode> handleOuterNodes(final LGraph layeredGraph) {
        List<LNode> remainingNodes = new ArrayList<>(layeredGraph.getLayerlessNodes().size());
        for (LNode node : layeredGraph.getLayerlessNodes()) {
            // Check if there is a layer constraint
            LayerConstraint layerConstraint = node.getProperty(LayeredOptions.LAYERING_LAYER_CONSTRAINT);
            EdgeConstraint edgeConstraint = null;

            switch (layerConstraint) {
            case FIRST:
            case FIRST_SEPARATE:
                edgeConstraint = EdgeConstraint.OUTGOING_ONLY;
                break;

            case LAST:
            case LAST_SEPARATE:
                edgeConstraint = EdgeConstraint.INCOMING_ONLY;
                break;
            }

            if (edgeConstraint != null) {
                // Set the edge constraint on the node
                node.setProperty(InternalProperties.EDGE_CONSTRAINT, EdgeConstraint.OUTGOING_ONLY);

                if (edgeConstraint == EdgeConstraint.INCOMING_ONLY) {
                    reverseEdges(layeredGraph, node, layerConstraint, PortType.INPUT);
                } else if (edgeConstraint == EdgeConstraint.OUTGOING_ONLY) {
                    reverseEdges(layeredGraph, node, layerConstraint, PortType.OUTPUT);
                }
            } else {
                // Remember this node for the next loop
                remainingNodes.add(node);
            }
        }

        return remainingNodes;
    }
    

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Inner Nodes

    /**
     * @param remainingNodes
     */
    private void handleInnerNodes(final LGraph layeredGraph, final List<LNode> remainingNodes) {
        // Iterate over the remaining nodes
        for (LNode node : remainingNodes) {
            // Check if there is a layer constraint
            LayerConstraint layerConstraint = node.getProperty(LayeredOptions.LAYERING_LAYER_CONSTRAINT);
            EdgeConstraint edgeConstraint = null;

            switch (layerConstraint) {
            case FIRST:
            case FIRST_SEPARATE:
                edgeConstraint = EdgeConstraint.OUTGOING_ONLY;
                break;

            case LAST:
            case LAST_SEPARATE:
                edgeConstraint = EdgeConstraint.INCOMING_ONLY;
                break;
            }

            if (edgeConstraint != null) {
                // Set the edge constraint on the node
                node.setProperty(InternalProperties.EDGE_CONSTRAINT, EdgeConstraint.OUTGOING_ONLY);

                if (edgeConstraint == EdgeConstraint.INCOMING_ONLY) {
                    reverseEdges(layeredGraph, node, layerConstraint, PortType.INPUT);
                } else if (edgeConstraint == EdgeConstraint.OUTGOING_ONLY) {
                    reverseEdges(layeredGraph, node, layerConstraint, PortType.OUTPUT);
                }
            } else {
                // If the port sides are fixed, but all ports are reversed, that probably means that we
                // have a feedback node. Normally, the connected edges would be routed around the node,
                // but that hides the feedback node character. We thus simply reverse all connected
                // edges and thus make ELK Layered think we have a regular node
                //
                // Note that this behavior is only desired if none of the connected nodes have
                // layer constraints set. Otherwise this processing causes issues with an external
                // port dummy with FIRST_SEPARATE and an inverted ports on the target node's EAST side.
                if (node.getProperty(LayeredOptions.PORT_CONSTRAINTS).isSideFixed() && !node.getPorts().isEmpty()) {

                    boolean allPortsReversed = true;

                    for (LPort port : node.getPorts()) {
                        // A non-reversed eastern port has more outgoing than incoming edges (net flow <= 0).
                        // A non-reversed western port has more incoming than outgoing edges (net flow >= 0).
                        // If this port doesn't violate one of these conditions, not all of the node's ports are
                        // reversed.
                        if (!(port.getSide() == PortSide.EAST && port.getNetFlow() > 0
                                || port.getSide() == PortSide.WEST && port.getNetFlow() < 0)) {

                            allPortsReversed = false;
                            break;
                        }

                        // no LAST or LAST_SEPARATE allowed for the target of outgoing edges
                        for (LEdge e : port.getOutgoingEdges()) {
                            LayerConstraint lc = e.getTarget().getNode().getProperty(
                                    LayeredOptions.LAYERING_LAYER_CONSTRAINT);
                            if (lc == LayerConstraint.LAST || lc == LayerConstraint.LAST_SEPARATE) {
                                allPortsReversed = false;
                                break;
                            }
                        }
                        // no FIRST or FIRST_SEPARATE allowed for the source of incoming edgs
                        for (LEdge e : port.getIncomingEdges()) {
                            LayerConstraint lc = e.getSource().getNode().getProperty(
                                    LayeredOptions.LAYERING_LAYER_CONSTRAINT);
                            if (lc == LayerConstraint.FIRST || lc == LayerConstraint.FIRST_SEPARATE) {
                                allPortsReversed = false;
                                break;
                            }

                        }
                    }

                    if (allPortsReversed) {
                        // Reverse all kinds of edges
                        reverseEdges(layeredGraph, node, layerConstraint, PortType.UNDEFINED);
                    }
                }
            }
        }
    }
    

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Edge Reversal

    /**
     * Reverses edges as appropriate.
     * 
     * @param layeredGraph
     *            the layered graph.
     * @param node
     *            the node to place in the layer.
     * @param nodeLayerConstraint
     *            the layer constraint put on the node.
     * @param targetPortType
     *            what ports should be turned into. If this is {@link PortType#OUTPUT}, for example, all incoming edges
     *            are reversed such that the ports become pure output ports.
     */
    private void reverseEdges(final LGraph layeredGraph, final LNode node, final LayerConstraint nodeLayerConstraint,
            final PortType targetPortType) {

        // Iterate through the node's edges and reverse them, if necessary
        for (LPort port : LGraphUtil.toPortArray(node.getPorts())) {
            // Only incoming edges
            if (targetPortType == PortType.INPUT || targetPortType == PortType.UNDEFINED) {
                LEdge[] outgoing = LGraphUtil.toEdgeArray(port.getOutgoingEdges());

                for (LEdge edge : outgoing) {
                    // Reverse the edge if we're allowed to do so
                    if (canReverseOutgoingEdge(nodeLayerConstraint, edge)) {
                        edge.reverse(layeredGraph, true);
                    }
                }
            }

            // Only outgoing edges
            if (targetPortType == PortType.OUTPUT || targetPortType == PortType.UNDEFINED) {
                LEdge[] incoming = LGraphUtil.toEdgeArray(port.getIncomingEdges());

                for (LEdge edge : incoming) {
                    // Reverse the edge if we're allowed to do so
                    if (canReverseIncomingEdge(nodeLayerConstraint, edge)) {
                        edge.reverse(layeredGraph, true);
                    }
                }
            }
        }
    }
    

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Reversal Testing

    /**
     * Checks whether or not a given edge outgoing edge can actually be reversed. It cannot be reversed if it already
     * has been, or if it connects a node in the {@code LAST} layer to either a node in the {@code LAST_SEPARATE} layer
     * or to a label dummy node.
     * 
     * @param sourceNodeLayerConstraint
     *            the source node's layer constraint.
     * @param edge
     *            the edge to possibly be reversed.
     * @return {@code true} if it's okay to reverse the edge.
     */
    private boolean canReverseOutgoingEdge(final LayerConstraint sourceNodeLayerConstraint, final LEdge edge) {
        // The layer constraint that gets passed to us
        assert sourceNodeLayerConstraint == edge.getSource().getNode()
                .getProperty(LayeredOptions.LAYERING_LAYER_CONSTRAINT);

        // If the edge is already reversed, we don't want to reverse it again
        if (edge.getProperty(InternalProperties.REVERSED)) {
            return false;
        }

        LNode targetNode = edge.getTarget().getNode();
        
        // If the node is supposed to be in the lAST layer...
        if (sourceNodeLayerConstraint == LayerConstraint.LAST) {
            // ...and is connected to a label dummy, we won't reverse it since that dummy has to be placed between the
            // LAST and the LAST_SEPARATE layer
            if (targetNode.getType() == NodeType.LABEL) {
                return false;
            }
        }
        
        // If the target is a node in the LAST_SEPARATE layer, we won't reverse it
        LayerConstraint targetLayerConstraint = targetNode.getProperty(LayeredOptions.LAYERING_LAYER_CONSTRAINT);
        if (targetLayerConstraint == LayerConstraint.LAST_SEPARATE) {
            return false;
        }

        return true;
    }

    /**
     * Checks whether or not a given edge incoming edge can actually be reversed. It cannot be reversed if it already
     * has been, or if it connects a node in the {@code FIRST} layer to either a node in the {@code FIRST_SEPARATE}
     * layer or to a label dummy node.
     * 
     * @param targetNodeLayerConstraint
     *            the target node's layer constraint.
     * @param edge
     *            the edge to possibly be reversed.
     * @return {@code true} if it's okay to reverse the edge.
     */
    private boolean canReverseIncomingEdge(final LayerConstraint targetNodeLayerConstraint, final LEdge edge) {
        // The layer constraint that gets passed to us
        assert targetNodeLayerConstraint == edge.getTarget().getNode()
                .getProperty(LayeredOptions.LAYERING_LAYER_CONSTRAINT);

        // If the edge is already reversed, we don't want to reverse it again
        if (edge.getProperty(InternalProperties.REVERSED)) {
            return false;
        }

        LNode sourceNode = edge.getSource().getNode();
        
        // If the node is supposed to be in the FIRST layer...
        if (targetNodeLayerConstraint == LayerConstraint.FIRST) {
            // ...and is connected to a label dummy, we won't reverse it since that dummy has to be placed between the
            // FIRST_SEPARATE and the FIRST layer
            if (sourceNode.getType() == NodeType.LABEL) {
                return false;
            }
        }
        
        // If the source is a node in the FIRST_SEPARATE layer, we won't reverse it
        LayerConstraint sourceLayerConstraint = sourceNode.getProperty(LayeredOptions.LAYERING_LAYER_CONSTRAINT);
        if (sourceLayerConstraint == LayerConstraint.FIRST_SEPARATE) {
            return false;
        }

        return true;
    }

}
