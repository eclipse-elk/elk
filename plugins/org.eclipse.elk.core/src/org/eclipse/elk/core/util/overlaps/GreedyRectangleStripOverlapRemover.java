/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the Eclipse Public License v1.0
 * which accompanies this distribution, and is available at
 * http://www.eclipse.org/legal/epl-v10.html
 *
 * Contributors:
 *    Kiel University - initial API and implementation
 *******************************************************************************/
package org.eclipse.elk.core.util.overlaps;

import java.util.Set;

import org.eclipse.elk.core.math.ElkRectangle;
import org.eclipse.elk.core.util.overlaps.RectangleStripOverlapRemover.RectangleNode;

import com.google.common.collect.Sets;

/**
 * Removes rectangle overlaps by greedily choosing the smallest y position that won't cause overlaps. Fast algorithm,
 * but not necessarily optimal in terms of resulting rectangle strip size.
 */
public final class GreedyRectangleStripOverlapRemover implements IRectangleStripOverlapRemovalStrategy {

    @Override
    public double removeOverlaps(final RectangleStripOverlapRemover overlapRemover) {
        final double gap = overlapRemover.getGap();
        Set<RectangleNode> alreadyPlacedNodes = Sets.newHashSet();
        double stripSize = 0;
        
        for (RectangleNode currNode : overlapRemover.getRectangleNodes()) {
            // We start with an initial y coordinate of zero
            double yPos = 0;
            
            // We now iterate over every conflicting rectangle node that we have already placed to see if the current
            // coordinate would cause an overlap
            for (RectangleNode overlapNode : currNode.getOverlappingNodes()) {
                if (alreadyPlacedNodes.contains(overlapNode)) {
                    ElkRectangle currRect = currNode.getRectangle();
                    ElkRectangle overlapRect = overlapNode.getRectangle();
                    
                    // Check if the current y coordinate would cause an overlap with the overlap node
                    if (yPos < overlapRect.y + overlapRect.height + gap
                            && yPos + currRect.height + gap > overlapRect.y) {
                        
                        yPos = overlapRect.y + overlapRect.height + gap;
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

}
