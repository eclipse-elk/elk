/*******************************************************************************
 * Copyright (c) 2020 Kiel University and others.
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
import java.util.ListIterator;

import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LNode.NodeType;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.LayerConstraint;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.core.UnsupportedConfigurationException;
import org.eclipse.elk.core.alg.ILayoutProcessor;
import org.eclipse.elk.core.util.IElkProgressMonitor;
import org.eclipse.elk.graph.properties.IProperty;
import org.eclipse.elk.graph.properties.Property;

/**
 * Hides nodes with {@link LayerConstraint#FIRST_SEPARATE} and {@link LayerConstraint#LAST_SEPARATE} from the graph.
 * These nodes would influence the layering and would cause strange results once they are moved into separate layers, so
 * we simply hide them for the time being and restore them once layering has finished.
 * 
 * <dl>
 *   <dt>Precondition:</dt>
 *     <dd>an unlayered graph.</dd>
 *     <dd>nodes with {@code FIRST_SEPARATE} layer constraint have only outgoing edges.</dd>
 *     <dd>nodes with {@code LAST_SEPARATE} layer constraint have only incoming edges.</dd>
 *   <dt>Postcondition:</dt>
 *     <dd>nodes with layer constraints {@link LayerConstraint#FIRST_SEPARATE} or {@link LayerConstraint#LAST_SEPARATE}
 *       are hidden from the graph.</dd>
 *     <dd>hidden nodes are saved in the graph through {@link InternalProperties#HIDDEN_NODES}.
 *   <dt>Slots:</dt>
 *     <dd>Before phase 2.</dd>
 *   <dt>Same-slot dependencies:</dt>
 *     <dd>{@link LabelDummyInserter}</dd>
 * </dl>
 * 
 * @see EdgeAndLayerConstraintEdgeReverser
 * @see LayerConstraintPostprocessor
 */
public final class LayerConstraintPreprocessor implements ILayoutProcessor<LGraph> {
    
    @Override
    public void process(final LGraph layeredGraph, final IElkProgressMonitor monitor) {
        monitor.begin("Layer constraint preprocessing", 1);
        
        List<LNode> hiddenNodes = new ArrayList<>();
        
        ListIterator<LNode> nodeIterator = layeredGraph.getLayerlessNodes().listIterator();
        while (nodeIterator.hasNext()) {
            LNode lNode = nodeIterator.next();
            
            if (isRelevantNode(lNode)) {
                hide(lNode);
                hiddenNodes.add(lNode);
                nodeIterator.remove();
            }
        }
        
        if (!hiddenNodes.isEmpty()) {
            layeredGraph.setProperty(InternalProperties.HIDDEN_NODES, hiddenNodes);
        }
        
        monitor.done();
    }

    /**
     * Checks if the given node has one of the layer constraints set that we're interested in here.
     */
    private boolean isRelevantNode(final LNode lNode) {
        switch (lNode.getProperty(LayeredOptions.LAYERING_LAYER_CONSTRAINT)) {
        case FIRST_SEPARATE:
        case LAST_SEPARATE:
            return true;
            
        default:
            return false;
        }
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Node Hiding
    
    /**
     * Disconnects the node from its adjacent nodes.
     */
    private void hide(final LNode lNode) {
        ensureNoInacceptableEdges(lNode);
        
        for (LEdge lEdge : lNode.getConnectedEdges()) {
            hide(lNode, lEdge);
        }
    }

    /**
     * Hides the edge of the given node. The node is the one that's currently being hidden.
     */
    private void hide(final LNode lNode, final LEdge lEdge) {
        boolean isOutgoing = lEdge.getSource().getNode() == lNode;
        LPort oppositePort = isOutgoing
                ? lEdge.getTarget()
                : lEdge.getSource();
        
        // Hide the edge and remember its original target or source
        if (isOutgoing) {
            lEdge.setTarget(null);
        } else {
            lEdge.setSource(null);
        }
        
        lEdge.setProperty(InternalProperties.ORIGINAL_OPPOSITE_PORT, oppositePort);
        
        // We may want to control where the other node will end up in the layering...
        updateOppositeNodeLayerConstraints(lNode, oppositePort.getNode());
    }
    
    /**
     * Used to keep track of which kinds of hidden nodes a given node had connections to.
     */
    private static enum HiddenNodeConnections {
        /** No connections to hidden nodes found yet. */
        NONE,
        /** Only connections to {@link LayerConstraint#FIRST_SEPARATE} nodes found so far. */
        FIRST_SEPARATE,
        /** Only connections to {@link LayerConstraint#LAST_SEPARATE} nodes found so far. */
        LAST_SEPARATE,
        /** Connections to both kinds of hidden nodes found. */
        BOTH;
        
        private HiddenNodeConnections combine(final LayerConstraint layerConstraint) {
            switch (this) {
            case NONE:
                return layerConstraint == LayerConstraint.FIRST_SEPARATE
                    ? FIRST_SEPARATE
                    : LAST_SEPARATE;
            
            case FIRST_SEPARATE:
                return layerConstraint == LayerConstraint.FIRST_SEPARATE
                    ? FIRST_SEPARATE
                    : BOTH;
            
            case LAST_SEPARATE:
                return layerConstraint == LayerConstraint.FIRST_SEPARATE
                    ? BOTH
                    : LAST_SEPARATE;
            
            default:
                return BOTH;
            }
        }
    }
    
    /** We'll use this property to keep track of each node's connections to hidden nodes. */
    private static final IProperty<HiddenNodeConnections> HIDDEN_NODE_CONNECTIONS = new Property<>(
            "separateLayerConnections", HiddenNodeConnections.NONE);

    /**
     * Possibly sets additional layer constraints on the opposite node if it doesn't have any yet. If the opposite node
     * only had the hidden node as a connection, the layering may end up not being nice. For instance, if the hidden
     * node is {@link LayerConstraint#FIRST_SEPARATE} and is the only node the opposite node is connected to, then the
     * opposite node will now be all alone. This may mean that it will end up in the last layer, far away from the
     * hidden node. We would thus set {@link LayerConstraint#FIRST} on the opposite node to keep that from happening.
     * 
     * <p>
     * We only do this, however, if the opposite node had connections only to {@link LayerConstraint#FIRST_SEPARATE} or
     * only to {@link LayerConstraint#LAST_SEPARATE} nodes. If it had connections to both, we don't do anything.
     */
    private void updateOppositeNodeLayerConstraints(final LNode hiddenNode, final LNode oppositeNode) {
        // If the node already has a layer constraint, don't bother
        if (oppositeNode.hasProperty(LayeredOptions.LAYERING_LAYER_CONSTRAINT)) {
            return;
        }
        
        // Update which kinds of hidden nodes the opposite node is connected to
        HiddenNodeConnections connections = oppositeNode.getProperty(HIDDEN_NODE_CONNECTIONS)
                .combine(hiddenNode.getProperty(LayeredOptions.LAYERING_LAYER_CONSTRAINT));
        oppositeNode.setProperty(HIDDEN_NODE_CONNECTIONS, connections);
        
        // If the node is still connected to things, don't bother
        if (oppositeNode.getConnectedEdges().iterator().hasNext()) {
            return;
        }
        
        // Our hidden node is the last node the opposite node was connected to. Check if we should add a constraint
        switch (connections) {
        case FIRST_SEPARATE:
            oppositeNode.setProperty(LayeredOptions.LAYERING_LAYER_CONSTRAINT, LayerConstraint.FIRST);
            break;
            
        case LAST_SEPARATE:
            oppositeNode.setProperty(LayeredOptions.LAYERING_LAYER_CONSTRAINT, LayerConstraint.LAST);
            break;
        }
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Checks

    /**
     * Ensures that the given node does not have incident edges that it must not have.
     */
    private void ensureNoInacceptableEdges(final LNode lNode) {
        LayerConstraint layerConstraint = lNode.getProperty(LayeredOptions.LAYERING_LAYER_CONSTRAINT);
        
        // Ensure the node only has acceptable incident edges
        if (layerConstraint == LayerConstraint.FIRST_SEPARATE) {
            for (LEdge inEdge : lNode.getIncomingEdges()) {
                if (!isAcceptableIncidentEdge(inEdge)) {
                    throw new UnsupportedConfigurationException("Node '" + lNode.getDesignation()
                            + "' has its layer constraint set to FIRST_SEPARATE, but has at least one incoming edge. "
                            + "FIRST_SEPARATE nodes must not have incoming edges.");
                }
            }
        } else if (layerConstraint == LayerConstraint.LAST_SEPARATE) {
            for (LEdge outEdge : lNode.getOutgoingEdges()) {
                if (!isAcceptableIncidentEdge(outEdge)) {
                    throw new UnsupportedConfigurationException("Node '" + lNode.getDesignation()
                            + "' has its layer constraint set to LAST_SEPARATE, but has at least one outgoing edge. "
                            + "LAST_SEPARATE nodes must not have outgoing edges.");
                }
            }
        }
    }
    
    /**
     * Checks whether the given edge is allowed.
     */
    private boolean isAcceptableIncidentEdge(final LEdge edge) {
        LNode sourceNode = edge.getSource().getNode();
        LNode targetNode = edge.getTarget().getNode();
        
        // If both nodes are external port dummies, that's fine
        return sourceNode.getType() == NodeType.EXTERNAL_PORT && targetNode.getType() == NodeType.EXTERNAL_PORT;
    }

}
