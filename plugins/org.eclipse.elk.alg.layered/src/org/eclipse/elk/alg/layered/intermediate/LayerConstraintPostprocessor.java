/*******************************************************************************
 * Copyright (c) 2010, 2020 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.layered.intermediate;

import java.util.List;
import java.util.ListIterator;

import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LGraphUtil;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.graph.LNode.NodeType;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.LayerConstraint;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.core.UnsupportedConfigurationException;
import org.eclipse.elk.core.alg.ILayoutProcessor;
import org.eclipse.elk.core.util.IElkProgressMonitor;

/**
 * Moves nodes with layer constraints to the appropriate layers. To meet the preconditions of
 * this processor, the {@link EdgeAndLayerConstraintEdgeReverser} can be used.
 * 
 * <dl>
 *   <dt>Precondition:</dt>
 *     <dd>a layered graph.</dd>
 *     <dd>nodes with {@link LayerConstraint#FIRST_SEPARATE} and {@link LayerConstraint#LAST_SEPARATE} were hidden from
 *       the graph and saved in {@link InternalProperties#HIDDEN_NODES}.</dd>
 *     <dd>nodes with {@link LayerConstraint#FIRST} have only outgoing edges, except for edges incoming from label
 *       dummy nodes.</dd>
 *     <dd>nodes with {@link LayerConstraint#LAST} have only incoming edges, except for edges going to label dummy
 *       nodes.</dd>
 *   <dt>Postcondition:</dt>
 *     <dd>nodes hidden by {@link LayerConstraintPreprocessor} are restored.</dd>
 *     <dd>nodes with layer constraints are placed in the appropriate layers.</dd>
 *   <dt>Slots:</dt>
 *     <dd>Before phase 3.</dd>
 *   <dt>Same-slot dependencies:</dt>
 *     <dd>{@link HierarchicalPortConstraintProcessor}</dd>
 *     <dd>{@link NodePromotion}</dd>
 * </dl>
 * 
 * @see EdgeAndLayerConstraintEdgeReverser
 * @see LayerConstraintPreprocessor
 */
public final class LayerConstraintPostprocessor implements ILayoutProcessor<LGraph> {

    @Override
    public void process(final LGraph layeredGraph, final IElkProgressMonitor monitor) {
        monitor.begin("Layer constraint postprocessing", 1);
        
        List<Layer> layers = layeredGraph.getLayers();
        
        // If the list of layers is non-empty, we may have FIRST and LAST nodes to move around
        if (!layers.isEmpty()) {
            // Retrieve the current first and last layers
            Layer firstLayer = layers.get(0);
            Layer lastLayer = layers.get(layers.size() - 1);
            
            // We may also need label dummy layers between the very first / last layers and the first / last layers
            Layer firstLabelLayer = new Layer(layeredGraph);
            Layer lastLabelLayer = new Layer(layeredGraph);

            // Move our nodes!
            moveFirstAndLastNodes(layeredGraph, firstLayer, lastLayer, firstLabelLayer, lastLabelLayer);
            
            // Add non-empty label layers
            if (!firstLabelLayer.getNodes().isEmpty()) {
                layers.add(0, firstLabelLayer);
            }
            
            if (!lastLabelLayer.getNodes().isEmpty()) {
                layers.add(lastLabelLayer);
            }
        }
        
        // Restore hidden FIRST_SEPARATE and LAST_SEPARATE nodes, if any
        if (layeredGraph.hasProperty(InternalProperties.HIDDEN_NODES)) {
            // Create the new first and last layers, in case they will be needed
            Layer firstSeparateLayer = new Layer(layeredGraph);
            Layer lastSeparateLayer = new Layer(layeredGraph);
            
            // Restore our nodes
            restoreHiddenNodes(layeredGraph, firstSeparateLayer, lastSeparateLayer);
            
            // Add non-empty separate layers
            if (!firstSeparateLayer.getNodes().isEmpty()) {
                layers.add(0, firstSeparateLayer);
            }
            
            if (!lastSeparateLayer.getNodes().isEmpty()) {
                layers.add(lastSeparateLayer);
            }
        }
        
        monitor.done();
    }

    /**
     * Moves nodes with {@link LayerConstraint#FIRST} and {@link LayerConstraint#LAST} to their respective layers and
     * removes layers that have become empty.
     */
    private void moveFirstAndLastNodes(final LGraph layeredGraph, final Layer firstLayer, final Layer lastLayer,
            final Layer firstLabelLayer, final Layer lastLabelLayer) {
        
        // We'll start by moving FIRST and LAST nodes (as well as connected label dummies) to their proper layers
        for (Layer layer : layeredGraph) {
            // Iterate through a node array to avoid ConcurrentModificationExceptions
            LNode[] nodes = LGraphUtil.toNodeArray(layer.getNodes());
            
            for (LNode node : nodes) {
                switch (node.getProperty(LayeredOptions.LAYERING_LAYER_CONSTRAINT)) {
                case FIRST:
                    throwUpUnlessNoIncomingEdges(node);
                    node.setLayer(firstLayer);
                    moveLabelsToLabelLayer(node, true, firstLabelLayer);
                    break;
                
                case LAST:
                    throwUpUnlessNoOutgoingEdges(node);
                    node.setLayer(lastLayer);
                    moveLabelsToLabelLayer(node, false, lastLabelLayer);
                    break;
                
                }
            }
        }
        
        // All the movement can cause layers to be empty. Instead of going to great lengths trying to prevent that in
        // the code above, we simply iterate over all layers and remove the empty ones
        ListIterator<Layer> layerIter = layeredGraph.getLayers().listIterator();
        while (layerIter.hasNext()) {
            if (layerIter.next().getNodes().isEmpty()) {
                layerIter.remove();
            }
        }
    }
    
    
    /**
     * Moves the label dummies coming in to the given node or going out from the given node to the given label dummy
     * layer.
     * 
     * @param node
     *            the node whose adjacent label dummies to move.
     * @param incoming
     *            {@code true} if label dummies on incoming edges should be moved, {@code false} if label dummies on
     *            outgoing edges should be moved.
     * @param labelLayer
     *            the layer to move the label dummies to.
     */
    private void moveLabelsToLabelLayer(final LNode node, final boolean incoming, final Layer labelLayer) {
        for (LEdge edge : incoming ? node.getIncomingEdges() : node.getOutgoingEdges()) {
            LNode possibleLableDummy = incoming
                    ? edge.getSource().getNode()
                    : edge.getTarget().getNode();
            
            if (possibleLableDummy.getType() == NodeType.LABEL) {
                possibleLableDummy.setLayer(labelLayer);
            }
        }
    }

    /**
     * Restores nodes hidden by {@link LayerConstraintPreprocessor}.
     */
    private void restoreHiddenNodes(final LGraph layeredGraph, final Layer firstSeparateLayer,
            final Layer lastSeparateLayer) {
        
        for (LNode hiddenNode : layeredGraph.getProperty(InternalProperties.HIDDEN_NODES)) {
            // Add the node to the appropriate layer
            switch (hiddenNode.getProperty(LayeredOptions.LAYERING_LAYER_CONSTRAINT)) {
            case FIRST_SEPARATE:
                hiddenNode.setLayer(firstSeparateLayer);
                break;
                
            case LAST_SEPARATE:
                hiddenNode.setLayer(lastSeparateLayer);
                break;
                
            default:
                // Only *_SEPARATE nodes should be in the list
                assert false;
            }
            
            // Restore the node's edges
            for (LEdge hiddenEdge : hiddenNode.getConnectedEdges()) {
                // A hidden node will usually be incident only to edges that were hidden together with the node.
                // There is one exception, though: if an edge connects two hidden nodes, only one will have been
                // responsible for hiding it. Thus, we may find an edge here that is already restored!
                if (hiddenEdge.getSource() != null && hiddenEdge.getTarget() != null) {
                    continue;
                }
                
                // One end point is the hidden node, the other is null
                boolean isOutgoing = hiddenEdge.getTarget() == null;
                
                LPort originalOppositePort = hiddenEdge.getProperty(InternalProperties.ORIGINAL_OPPOSITE_PORT);
                if (isOutgoing) {
                    hiddenEdge.setTarget(originalOppositePort);
                } else {
                    hiddenEdge.setSource(originalOppositePort);
                }
            }
        }
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Checks
    
    /**
     * Check that the node has no incoming edges except from label dummies, and fail if it has any.
     */
    private void throwUpUnlessNoIncomingEdges(final LNode node) {
        for (LEdge incoming : node.getIncomingEdges()) {
            if (incoming.getSource().getNode().getType() != NodeType.LABEL) {
                throw new UnsupportedConfigurationException("Node '" + node.getDesignation()
                        + "' has its layer constraint set to FIRST, but has at least one incoming edge that "
                        + " does not come from a FIRST_SEPARATE node. That must not happen.");
            }
        }
    }
    
    /**
     * Check that the node has no outgoing edges except to label dummies, and fail if it has any.
     */
    private void throwUpUnlessNoOutgoingEdges(final LNode node) {
        for (LEdge outgoing : node.getOutgoingEdges()) {
            if (outgoing.getTarget().getNode().getType() != NodeType.LABEL) {
                throw new UnsupportedConfigurationException("Node '" + node.getDesignation()
                        + "' has its layer constraint set to LAST, but has at least one outgoing edge that "
                        + " does not go to a LAST_SEPARATE node. That must not happen.");
            }
        }
    }

}
