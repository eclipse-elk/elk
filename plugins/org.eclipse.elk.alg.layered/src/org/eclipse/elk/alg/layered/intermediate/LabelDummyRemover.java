/*******************************************************************************
 * Copyright (c) 2012, 2020 Kiel University and others.
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
package org.eclipse.elk.alg.layered.intermediate;

import java.util.List;
import java.util.ListIterator;

import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LLabel;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LNode.NodeType;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.core.alg.ILayoutProcessor;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.options.Direction;
import org.eclipse.elk.core.options.EdgeRouting;
import org.eclipse.elk.core.options.LabelSide;
import org.eclipse.elk.core.util.IElkProgressMonitor;
import org.eclipse.elk.alg.layered.graph.LGraphUtil;


/**
 * <p>Processor that removes the inserted center label dummies and places the labels on their
 * position.</p>
 * 
 * <dl>
 *   <dt>Precondition:</dt>
 *     <dd>a layered graph<dd>
 *     <dd>nodes are placed</dd>
 *     <dd>edges are routed</dd>
 *     <dd>center labels are represented by and attached to center label dummy nodes.</dd>
 *   <dt>Postcondition:</dt>
 *     <dd>labels are placed</dd>
 *     <dd>there are no dummy nodes of type
 *       {@link org.eclipse.elk.alg.layered.options.NodeType#LABEL}.</dd>
 *     <dd>center labels are attached to their original edges again.</dd>
 *   <dt>Slots:</dt>
 *     <dd>After phase 5.</dd>
 *   <dt>Same-slot dependencies:</dt>
 *     <dd>{@link HierarchicalPortOrthogonalEdgeRouter}</dd>
 * </dl>
 * 
 * @see LabelDummyInserter
 */
public final class LabelDummyRemover implements ILayoutProcessor<LGraph> {

    @Override
    public void process(final LGraph layeredGraph, final IElkProgressMonitor monitor) {
        monitor.begin("Label dummy removal", 1);
        
        double edgeLabelSpacing = layeredGraph.getProperty(LayeredOptions.SPACING_EDGE_LABEL);
        double labelLabelSpacing = layeredGraph.getProperty(LayeredOptions.SPACING_LABEL_LABEL);
        Direction layoutDirection = layeredGraph.getProperty(LayeredOptions.DIRECTION);
        
        for (Layer layer : layeredGraph.getLayers()) {
            // An iterator is necessary for traversing nodes, since dummy nodes might be removed
            ListIterator<LNode> nodeIterator = layer.getNodes().listIterator();
            
            while (nodeIterator.hasNext()) {
                LNode node = nodeIterator.next();
                
                if (node.getType() == NodeType.LABEL) {
                    // First, place labels on position of dummy node 
                    LEdge originEdge = (LEdge) node.getProperty(InternalProperties.ORIGIN);
                    double thickness = originEdge.getProperty(LayeredOptions.EDGE_THICKNESS).doubleValue();
                    boolean labelsBelowEdge = node.getProperty(InternalProperties.LABEL_SIDE) == LabelSide.BELOW;
                    
                    KVector currLabelPos = new KVector(node.getPosition());
                    
                    // If the labels are to be placed below their edge, we need to move the first label's
                    // position down a bit to respect the label spacing
                    if (labelsBelowEdge) {
                        currLabelPos.y += thickness + edgeLabelSpacing;
                    }
                    
                    // Calculate the space available for the placement of labels
                    KVector labelSpace = new KVector(
                            node.getSize().x,
                            node.getSize().y
                                + (node.isInlineEdgeLabel() ? 0 : -thickness - edgeLabelSpacing));
                    
                    // Place labels
                    List<LLabel> representedLabels = node.getProperty(InternalProperties.REPRESENTED_LABELS);
                    
                    if (layoutDirection.isVertical()) {
                        placeLabelsForVerticalLayout(
                                representedLabels,
                                currLabelPos,
                                labelLabelSpacing,
                                labelSpace,
                                labelsBelowEdge,
                                layoutDirection);
                    } else {
                        placeLabelsForHorizontalLayout(
                                representedLabels,
                                currLabelPos,
                                labelLabelSpacing,
                                labelSpace);
                    }
                    
                    // Add represented labels back to the original edge
                    originEdge.getLabels().addAll(representedLabels);
                    
                    // Whether we need to add unnecessary bend points around the label dummy depends on the edge
                    // router. For orthogonal edge routing, they are not necessary. For splines, they may even be
                    // harmful. For polylines, they are necessary to keep the routes edges take close to their labels
                    LongEdgeJoiner.joinAt(node,
                            layeredGraph.getProperty(LayeredOptions.EDGE_ROUTING) == EdgeRouting.POLYLINE);
                    
                    // Remove the node
                    nodeIterator.remove();
                }
            }
        }
        
        monitor.done();
    }

    private void placeLabelsForHorizontalLayout(final List<LLabel> labels, final KVector labelPos,
            final double labelSpacing, final KVector labelSpace) {
        
        for (LLabel label : labels) {
            label.getPosition().x = labelPos.x + (labelSpace.x - label.getSize().x) / 2.0;
            label.getPosition().y = labelPos.y;
            
            labelPos.y += label.getSize().y + labelSpacing;
        }
    }

    private void placeLabelsForVerticalLayout(final List<LLabel> labels, final KVector labelPos,
            final double labelSpacing, final KVector labelSpace, final boolean leftAligned,
            final Direction layoutDirection) {
        
        // We may have to override the alignment if all labels here are inline labels
        boolean inline = labels.stream().allMatch(label -> label.getProperty(LayeredOptions.EDGE_LABELS_INLINE));
        
        // Due to the way layout directions work, we need to pay attention to the order in which we place labels. While
        // we can simply place them as they come for the DOWN direction, doing the same for the UP direction will
        // reverse the label order in the final result. Thus, in that case we iterate over the reversed label list
        List<LLabel> effectiveLabels = labels;
        if (layoutDirection == Direction.UP) {
            effectiveLabels = LGraphUtil.reversed(effectiveLabels);
        }
        
        for (LLabel label : effectiveLabels) {
            label.getPosition().x = labelPos.x;
            
            if (inline) {
                label.getPosition().y = labelPos.y + (labelSpace.y - label.getSize().y) / 2;
            } else if (leftAligned) {
                label.getPosition().y = labelPos.y;
            } else {
                label.getPosition().y = labelPos.y + labelSpace.y - label.getSize().y;
            }
            
            labelPos.x += label.getSize().x + labelSpacing;
        }
    }

}
