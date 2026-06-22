/*******************************************************************************
 * Copyright (c) 2012, 2017 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.alg.layered.intermediate;

import java.util.List;
import java.util.ListIterator;

import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LLabel;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LNode.NodeType;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.core.alg.ILayoutProcessor;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.options.Direction;
import org.eclipse.elk.core.options.EdgeLabelPlacement;
import org.eclipse.elk.core.options.PortConstraints;
import org.eclipse.elk.core.util.IElkProgressMonitor;

import java.util.function.Predicate;
import java.util.ArrayList;

/**
 * Processor that inserts dummy nodes into edges that have center labels to reserve space for them.
 * At most one dummy node is created for each edge.
 * 
 * <dl>
 *   <dt>Precondition:</dt>
 *     <dd>an unlayered acyclic graph.</dd>
 *   <dt>Postcondition:</dt>
 *     <dd>dummy nodes are inserted that represent center labels.</dd>
 *     <dd>center labels are removed from their respective edges and attached to the dummy nodes.</dd>
 *   <dt>Slots:</dt>
 *     <dd>Before phase 2.</dd>
 *   <dt>Same-slot dependencies:</dt>
 *     <dd>None.</dd>
 * </dl>
 */
public final class LabelDummyInserter implements ILayoutProcessor<LGraph> {
    
    /** Predicate that checks for center labels. */
    private static final Predicate<LLabel> CENTER_LABEL =
            label -> label.getProperty(LayeredOptions.EDGE_LABELS_PLACEMENT) == EdgeLabelPlacement.CENTER;
    
    @Override
    public void process(final LGraph layeredGraph, final IElkProgressMonitor monitor) {
        monitor.begin("Label dummy insertions", 1);
        
        // We cannot add the nodes to the graph while we're iterating over it, so remember the dummy nodes we create
        List<LNode> newDummyNodes = new ArrayList<>();
        
        double edgeLabelSpacing = layeredGraph.getProperty(LayeredOptions.SPACING_EDGE_LABEL);
        double labelLabelSpacing = layeredGraph.getProperty(LayeredOptions.SPACING_LABEL_LABEL);
        Direction layoutDirection = layeredGraph.getProperty(LayeredOptions.DIRECTION);

        for (LNode node : layeredGraph.getLayerlessNodes()) {
            for (LEdge edge : node.getOutgoingEdges()) {
                if (edgeNeedsToBeProcessed(edge)) {
                    double thickness = retrieveThickness(edge);
                    
                    // Create dummy node and remember represented labels (to be filled below)
                    List<LLabel> representedLabels = new ArrayList<>(edge.getLabels().size());
                    LNode dummyNode = createLabelDummy(layeredGraph, edge, thickness, representedLabels);
                    newDummyNodes.add(dummyNode);
                    
                    // Determine the size of the dummy node and move labels over to it
                    KVector dummySize = dummyNode.getSize();
                    
                    ListIterator<LLabel> iterator = edge.getLabels().listIterator();
                    while (iterator.hasNext()) {
                        LLabel label = iterator.next();
                        
                        if (label.getProperty(LayeredOptions.EDGE_LABELS_PLACEMENT) == EdgeLabelPlacement.CENTER) {
                            // The way we stack labels depends on the layout direction
                            if (layoutDirection.isVertical()) {
                                dummySize.x += label.getSize().x + labelLabelSpacing;
                                dummySize.y = Math.max(dummySize.y, label.getSize().y);
                            } else {
                                dummySize.x = Math.max(dummySize.x, label.getSize().x);
                                dummySize.y += label.getSize().y + labelLabelSpacing;
                            }
                            
                            // Move the label over to the dummy node's REPRESENTED_LABELS property
                            representedLabels.add(label);
                            iterator.remove();
                        }
                    }
                    
                    // The dummy node now contains a superfluous label-label spacing and does not include the
                    // edge-label spacing yet
                    if (layoutDirection.isVertical()) {
                        dummySize.x -= labelLabelSpacing;
                        dummySize.y += edgeLabelSpacing + thickness;
                    } else {
                        dummySize.y += edgeLabelSpacing - labelLabelSpacing + thickness;
                    }
                }
            }
        }

        // Add created dummies to graph
        layeredGraph.getLayerlessNodes().addAll(newDummyNodes);
        
        monitor.done();
    }
    
    /**
     * Checks whether the given edge needs to be processed. That's the case if the edge is not a self-loop
     * and if it has center edge labels in the first place.
     */
    private boolean edgeNeedsToBeProcessed(final LEdge edge) {
        return edge.getSource().getNode() != edge.getTarget().getNode()
                && java.util.stream.StreamSupport.stream(edge.getLabels().spliterator(), false).anyMatch(CENTER_LABEL);
    }
    
    /**
     * Retrieves the given edge's thickness. If this is a negative value, zero is returned and set on the edge.
     */
    private double retrieveThickness(final LEdge edge) {
        double thickness = edge.getProperty(LayeredOptions.EDGE_THICKNESS);
        if (thickness < 0) {
            thickness = 0;
            edge.setProperty(LayeredOptions.EDGE_THICKNESS, thickness);
        }
        
        return thickness;
    }

    /**
     * Creates a label dummy for the given edge.
     * 
     * @param layeredGraph
     *            graph the dummy will later be placed in.
     * @param edge
     *            the edge the label dummy is created for.
     * @param thickness
     *            the edge's thickness.
     * @param representedLabels
     *            currently empty list of labels represented by the new label dummy. This is set on the edge as a
     *            property and will later be filled with the represented labels by the calling method.
     */
    private LNode createLabelDummy(final LGraph layeredGraph, final LEdge edge, final double thickness,
            final List<LLabel> representedLabels) {
        
        LNode dummyNode = new LNode(layeredGraph);
        dummyNode.setType(NodeType.LABEL);
        
        dummyNode.setProperty(InternalProperties.ORIGIN, edge);
        dummyNode.setProperty(InternalProperties.REPRESENTED_LABELS, representedLabels);
        dummyNode.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_POS);
        dummyNode.setProperty(InternalProperties.LONG_EDGE_SOURCE, edge.getSource());
        dummyNode.setProperty(InternalProperties.LONG_EDGE_TARGET, edge.getTarget());
        
        
        // Actually split the edge
        LongEdgeSplitter.splitEdge(edge, dummyNode);
        
        // Place ports at the edge's center
        double portPos = Math.floor(thickness / 2);
        for (LPort dummyPort : dummyNode.getPorts()) {
            dummyPort.getPosition().y = portPos;
        }
        
        return dummyNode;
    }

}
