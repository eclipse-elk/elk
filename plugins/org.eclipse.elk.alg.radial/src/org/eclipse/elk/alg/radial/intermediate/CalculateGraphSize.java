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
package org.eclipse.elk.alg.radial.intermediate;

import org.eclipse.elk.alg.radial.InternalProperties;
import org.eclipse.elk.alg.radial.options.RadialOptions;
import org.eclipse.elk.core.alg.ILayoutProcessor;
import org.eclipse.elk.core.math.ElkMargin;
import org.eclipse.elk.core.math.ElkPadding;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.util.IElkProgressMonitor;
import org.eclipse.elk.graph.ElkNode;

/**
 * Calculate the size of the graph and shift nodes into positive coordinates if necessary.
 *
 */
public class CalculateGraphSize implements ILayoutProcessor<ElkNode> {

    /** Shift the nodes such that each nodes has x and y coordinates bigger 0. */
    public void process(final ElkNode graph, final IElkProgressMonitor progressMonitor) {
        progressMonitor.begin("Calculate Graph Size", 1);
        progressMonitor.logGraph(graph, "Before");
        // calculate the offset from border spacing and node distribution
        double minXPos = Double.MAX_VALUE;
        double minYPos = Double.MAX_VALUE;
        double maxXPos = Double.MIN_VALUE;
        double maxYPos = Double.MIN_VALUE;

        for (ElkNode node : graph.getChildren()) {
            double posX = node.getX();
            double posY = node.getY();
            double width = node.getWidth();
            double height = node.getHeight();
            ElkMargin margins = node.getProperty(CoreOptions.MARGINS);

            minXPos = Math.min(minXPos, posX - margins.left);
            minYPos = Math.min(minYPos, posY - margins.top);
            maxXPos = Math.max(maxXPos, posX + width + margins.right);
            maxYPos = Math.max(maxYPos, posY + height + margins.bottom);
        }

        ElkPadding padding = graph.getProperty(CoreOptions.PADDING);
        KVector offset = new KVector(minXPos - padding.getLeft(), minYPos - padding.getTop());
        
        
        double width = maxXPos - minXPos + padding.getHorizontal();
        double height = maxYPos - minYPos + padding.getVertical();
        
        if (graph.getProperty(RadialOptions.CENTER_ON_ROOT)) {
            ElkNode root = graph.getProperty(InternalProperties.ROOT_NODE);
            ElkMargin rootMargins = root.getProperty(CoreOptions.MARGINS);
            // calculate the current midpoint of the root, taking into account the defined margins and the already
            // calculated offset necessary to shift the graph into the positive quadrant of the coordinate system
            double rootX = root.getX() + root.getWidth()/2 + (rootMargins.left + rootMargins.right)/2 - offset.x;
            double rootY = root.getY() + root.getHeight()/2 + (rootMargins.top + rootMargins.bottom)/2 - offset.y;
            
            double dx = width - rootX;
            double dy = height - rootY;
            
            if (dx < width / 2) {
                // need to add additional space on the left
                double additionalX = dx - rootX;
                width += additionalX;
                offset.x -= additionalX;
            } else {
                // add addtional space on the right
                double additionalX = rootX - dx;
                width += additionalX;
            }
            
            if (dy < height / 2) {
                //need to add additional space on the top
                double additionalY = dy - rootY;
                height += additionalY;
                offset.y -= additionalY;
            } else {
                // add addtional space on the bottom
                double additionalY = rootY - dy;
                height += additionalY;
            }
            
        }

        // process the nodes
        for (ElkNode node : graph.getChildren()) {
            // set the node position
            node.setX(node.getX() - offset.x);
            node.setY(node.getY() - offset.y);
        }

        // set up the graph
        if (!graph.getProperty(CoreOptions.NODE_SIZE_FIXED_GRAPH_SIZE)) {
            graph.setWidth(width);
            graph.setHeight(height);
        }

        // store child area info
        graph.setProperty(CoreOptions.CHILD_AREA_WIDTH, width - padding.getHorizontal());
        graph.setProperty(CoreOptions.CHILD_AREA_HEIGHT, height - padding.getVertical());
        progressMonitor.logGraph(graph, "After");
    }
}
