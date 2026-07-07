/*******************************************************************************
 * Copyright (c) 2011, 2020 Kiel University and others.
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

import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LMargin;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.layered.options.Spacings;
import org.eclipse.elk.core.alg.ILayoutProcessor;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.util.IElkProgressMonitor;

import com.google.common.collect.Lists;

/**
 * Post-processor for comment boxes. If any comments are found that were removed by the
 * {@link CommentPreprocessor}, they are reinserted and placed above or below their
 * corresponding connected node. This requires the margin around the node to be large
 * enough to hold all comments, which is ensured by the {@link InnermostNodeMarginCalculator}.
 * 
 * <dl>
 *   <dt>Precondition:</dt>
 *      <dd>Comments have been processed by {@link CommentPreprocessor}.</dd>
 *      <dd>Nodes are organized in layers.</dd>
 *      <dd>Nodes have been placed with enough spacing to hold their connected comment boxes.</dd>
 *   <dt>Postcondition:</dt>
 *      <dd>Comments that have been removed by pre-processing are reinserted properly in the graph.</dd>
 *   <dt>Slots:</dt>
 *      <dd>After phase 5.</dd>
 * </dl>
 */
public final class CommentPostprocessor implements ILayoutProcessor<LGraph> {

    @Override
    public void process(final LGraph layeredGraph, final IElkProgressMonitor monitor) {
        monitor.begin("Comment post-processing", 1);
        
        for (Layer layer : layeredGraph) {
            List<LNode> boxes = Lists.newArrayList();
            for (LNode node : layer) {
                List<LNode> topBoxes = node.getProperty(InternalProperties.TOP_COMMENTS);
                List<LNode> bottomBoxes = node.getProperty(InternalProperties.BOTTOM_COMMENTS);
                
                if (topBoxes != null || bottomBoxes != null) {
                    process(node, topBoxes, bottomBoxes);
                    
                    if (topBoxes != null) {
                        boxes.addAll(topBoxes);
                    }
                    
                    if (bottomBoxes != null) {
                        boxes.addAll(bottomBoxes);
                    }
                }
            }
            
            layer.getNodes().addAll(boxes);
        }
        
        monitor.done();
    }
    
    /**
     * Process a node with its connected comment boxes.
     * 
     * @param node a normal node
     * @param topBoxes a list of boxes to be placed on top, or {@code null}
     * @param bottomBoxes a list of boxes to be placed in the bottom, or {@code null}
     */
    private void process(final LNode node, final List<LNode> topBoxes, final List<LNode> bottomBoxes) {
        KVector nodePos = node.getPosition();
        KVector nodeSize = node.getSize();
        LMargin margin = node.getMargin();

        double commentCommentSpacing = Spacings.getIndividualOrDefault(node, LayeredOptions.SPACING_COMMENT_COMMENT);
        
        if (topBoxes != null) {
            // determine the total width and maximal height of the top boxes
            double boxesWidth = commentCommentSpacing * (topBoxes.size() - 1);
            double maxHeight = 0;
            for (LNode box : topBoxes) {
                boxesWidth += box.getSize().x;
                maxHeight = Math.max(maxHeight, box.getSize().y);
            }
            
            // place the boxes on top of the node, horizontally centered around the node itself
            double x = nodePos.x - (boxesWidth - nodeSize.x) / 2;
            double baseLine = nodePos.y - margin.top + maxHeight;
            double anchorInc = nodeSize.x / (topBoxes.size() + 1);
            double anchorX = anchorInc;
            for (LNode box : topBoxes) {
                box.getPosition().x = x;
                box.getPosition().y = baseLine - box.getSize().y;
                x += box.getSize().x + commentCommentSpacing;
                // set source and target point for the connecting edge
                LPort boxPort = getBoxPort(box);
                boxPort.getPosition().x = box.getSize().x / 2 - boxPort.getAnchor().x;
                boxPort.getPosition().y = box.getSize().y;
                LPort nodePort = box.getProperty(InternalProperties.COMMENT_CONN_PORT);
                if (nodePort.getDegree() == 1) {
                    nodePort.getPosition().x = anchorX - nodePort.getAnchor().x;
                    nodePort.getPosition().y = 0;
                    nodePort.setNode(node);
                }
                anchorX += anchorInc;
            }
        }

        if (bottomBoxes != null) {
            // determine the total width and maximal height of the bottom boxes
            double boxesWidth = commentCommentSpacing * (bottomBoxes.size() - 1);
            double maxHeight = 0;
            for (LNode box : bottomBoxes) {
                boxesWidth += box.getSize().x;
                maxHeight = Math.max(maxHeight, box.getSize().y);
            }
            
            // place the boxes in the bottom of the node, horizontally centered around the node itself
            double x = nodePos.x - (boxesWidth - nodeSize.x) / 2;
            double baseLine = nodePos.y + nodeSize.y + margin.bottom - maxHeight;
            double anchorInc = nodeSize.x / (bottomBoxes.size() + 1);
            double anchorX = anchorInc;
            for (LNode box : bottomBoxes) {
                box.getPosition().x = x;
                box.getPosition().y = baseLine;
                x += box.getSize().x + commentCommentSpacing;
                // set source and target point for the connecting edge
                LPort boxPort = getBoxPort(box);
                boxPort.getPosition().x = box.getSize().x / 2 - boxPort.getAnchor().x;
                boxPort.getPosition().y = 0;
                LPort nodePort = box.getProperty(InternalProperties.COMMENT_CONN_PORT);
                if (nodePort.getDegree() == 1) {
                    nodePort.getPosition().x = anchorX - nodePort.getAnchor().x;
                    nodePort.getPosition().y = nodeSize.y;
                    nodePort.setNode(node);
                }
                anchorX += anchorInc;
            }
        }
    }
    
    /**
     * Retrieve the port of the given comment box that connects it with the
     * corresponding node. The box is expected to have exactly one connection, which
     * is checked by the {@link CommentPreprocessor}.
     * 
     * @param commentBox a comment box
     * @return the port that connects the box
     */
    private LPort getBoxPort(final LNode commentBox) {
        LPort nodePort = commentBox.getProperty(InternalProperties.COMMENT_CONN_PORT);
        for (LPort port : commentBox.getPorts()) {
            for (LEdge edge : port.getOutgoingEdges()) {
                // reconnect the edge (has been disconnected by pre-processor)
                edge.setTarget(nodePort);
                return port;
            }
            for (LEdge edge : port.getIncomingEdges()) {
                // reconnect the edge (has been disconnected by pre-processor)
                edge.setSource(nodePort);
                return port;
            }
        }
        return null;
    }

}
