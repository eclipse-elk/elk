/*******************************************************************************
 * Copyright (c) 2010, 2020 Kiel University and others.
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

import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LMargin;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.layered.options.Spacings;
import org.eclipse.elk.core.alg.ILayoutProcessor;
import org.eclipse.elk.core.util.IElkProgressMonitor;

/**
 * Computes and sets the node margins required to place comment boxes.
 * 
 * <dl>
 *   <dt>Preconditions:</dt>
 *     <dd>A layered graph.</dd>
 *     <dd>Node margins include space for ports, port labels, and self loops.</dd>
 *   <dt>Postcondition:</dt>
 *     <dd>Node margins are extended to include space for comment boxes.</dd>
 *   <dt>Slots:</dt>
 *     <dd>Before phase 4.</dd>
 *   <dt>Same-slot dependencies:</dt>
 *     <dd>{@link SelfLoopRouter}</dd>
 * </dl>
 * 
 * @see InnermostNodeMarginCalculator
 */
public final class CommentNodeMarginCalculator implements ILayoutProcessor<LGraph> {

    @Override
    public void process(final LGraph layeredGraph, final IElkProgressMonitor monitor) {
        monitor.begin("Node margin calculation", 1);

        // Iterate through the layers to additionally handle comments
        layeredGraph.getLayers().stream()
            .flatMap(layer -> layer.getNodes().stream())
            .forEach(lnode -> processComments(lnode));

        monitor.done();
    }

    /**
     * Make some extra space for comment boxes that are placed near the given node.
     */
    private void processComments(final LNode node) {
        LMargin margin = node.getMargin();
        
        List<LNode> topBoxes = node.getProperty(InternalProperties.TOP_COMMENTS);
        List<LNode> bottomBoxes = node.getProperty(InternalProperties.BOTTOM_COMMENTS);
        
        if (topBoxes == null && bottomBoxes == null) {
            // Shortcut if there are no attached comments
            return;
        }
        
        // Retrieve the spacings that apply to this node
        double commentCommentSpacing = Spacings.getIndividualOrDefault(node, LayeredOptions.SPACING_COMMENT_COMMENT);
        double commentNodeSpacing = Spacings.getIndividualOrDefault(node, LayeredOptions.SPACING_COMMENT_NODE);
        
        // Consider comment boxes that are put on top of the node
        double topWidth = 0;
        if (topBoxes != null) {
            double maxHeight = 0;
            for (LNode commentBox : topBoxes) {
                maxHeight = Math.max(maxHeight, commentBox.getSize().y);
                topWidth += commentBox.getSize().x;
            }
            topWidth += commentCommentSpacing * (topBoxes.size() - 1);
            margin.top += maxHeight + commentNodeSpacing;
        }

        // Consider comment boxes that are put in the bottom of the node
        double bottomWidth = 0;
        if (bottomBoxes != null) {
            double maxHeight = 0;
            for (LNode commentBox : bottomBoxes) {
                maxHeight = Math.max(maxHeight, commentBox.getSize().y);
                bottomWidth += commentBox.getSize().x;
            }
            bottomWidth += commentCommentSpacing * (bottomBoxes.size() - 1);
            margin.bottom += maxHeight + commentNodeSpacing;
        }

        // Check if the maximum width of the comments is wider than the node itself, which the comments
        // are centered on
        double maxCommentWidth = Math.max(topWidth, bottomWidth);
        if (maxCommentWidth > node.getSize().x) {
            double protrusion = (maxCommentWidth - node.getSize().x) / 2;
            margin.left = Math.max(margin.left, protrusion);
            margin.right = Math.max(margin.right, protrusion);
        }
    }

}
