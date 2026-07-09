/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
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
package org.eclipse.elk.alg.common.overlaps;

import java.util.Set;

import org.eclipse.elk.alg.common.overlaps.RectangleStripOverlapRemover.RectangleNode;
import org.eclipse.elk.core.math.ElkRectangle;

import com.google.common.collect.Sets;

/**
 * Removes rectangle overlaps by greedily choosing the smallest y position that won't cause overlaps. Fast algorithm,
 * but not necessarily optimal in terms of resulting rectangle strip size.
 */
public final class GreedyRectangleStripOverlapRemover implements IRectangleStripOverlapRemovalStrategy {

    @Override
    public double removeOverlaps(final RectangleStripOverlapRemover overlapRemover) {
        final double verticalGap = overlapRemover.getVerticalGap();
        Set<RectangleNode> alreadyPlacedNodes = Sets.newHashSet();
        double stripSize = 0;
        
        for (RectangleNode currNode : overlapRemover.getRectangleNodes()) {
            // We start with an initial y coordinate of zero
            double yPos = 0;
            
            // Sort the node's list of overlapping nodes by y coordinate
            currNode.getOverlappingNodes().sort(GreedyRectangleStripOverlapRemover::compareByYCoordinate);
            
            // We now iterate over every conflicting rectangle node that we have already placed to see if the current
            // coordinate would cause an overlap
            for (RectangleNode overlapNode : currNode.getOverlappingNodes()) {
                if (alreadyPlacedNodes.contains(overlapNode)) {
                    ElkRectangle currRect = currNode.getRectangle();
                    ElkRectangle overlapRect = overlapNode.getRectangle();
                    
                    // Check if the current y coordinate would cause an overlap with the overlap node
                    if (yPos < overlapRect.y + overlapRect.height + verticalGap
                            && yPos + currRect.height + verticalGap > overlapRect.y) {
                        
                        yPos = overlapRect.y + overlapRect.height + verticalGap;
                    }
                }
            }
            
            // Apply the y coordinate and remember that this node is now placed
            currNode.getRectangle().y = yPos;
            alreadyPlacedNodes.add(currNode);
            
            // Update the strip size
            stripSize = Math.max(stripSize, currNode.getRectangle().y + currNode.getRectangle().height);
        }
        
        return stripSize;
    }
    
    /**
     * Compares two rectangle nodes by the y coordinates of their rectangles.
     * 
     * @param node1
     *            first node.
     * @param node2
     *            second node.
     */
    public static int compareByYCoordinate(final RectangleNode node1, final RectangleNode node2) {
        return Double.compare(node1.getRectangle().y, node2.getRectangle().y);
    }

}
