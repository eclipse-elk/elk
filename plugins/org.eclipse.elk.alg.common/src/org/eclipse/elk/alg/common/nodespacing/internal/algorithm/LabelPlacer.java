/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.common.nodespacing.internal.algorithm;

import org.eclipse.elk.alg.common.nodespacing.cellsystem.StripContainerCell;
import org.eclipse.elk.alg.common.nodespacing.internal.NodeContext;
import org.eclipse.elk.core.math.ElkRectangle;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.options.PortSide;
import org.eclipse.elk.core.options.SizeOptions;

/**
 * Knows how to properly size and position outer node label containers and to place node and port labels.
 */
public final class LabelPlacer {
    
    /**
     * No instance required.
     */
    private LabelPlacer() {
        
    }
    
    
    /**
     * Places outer node label containers as well as all labels. ALL OF THEM!!!
     */
    public static void placeLabels(final NodeContext nodeContext) {
        // Properly place all label cells for outer node labels
        placeOuterNodeLabelContainers(nodeContext);
        
        // Tell all node label cells to place their labels
        nodeContext.nodeLabelCells.values().stream()
                .forEach(labelCell -> labelCell.applyLabelLayout());
        
        // Tell all port label cells to place their labels
        nodeContext.portContexts.values().stream()
                .filter(portContext -> portContext.portLabelCell != null)
                .forEach(portContext -> portContext.portLabelCell.applyLabelLayout());
    }

    /**
     * Places the outer node label containers according to size constraints and the current node size.
     */
    private static void placeOuterNodeLabelContainers(final NodeContext nodeContext) {
        boolean outerNodeLabelsOverhang = nodeContext.sizeOptions.contains(SizeOptions.OUTSIDE_NODE_LABELS_OVERHANG);
        
        placeHorizontalOuterNodeLabelContainer(nodeContext, outerNodeLabelsOverhang, PortSide.NORTH);
        placeHorizontalOuterNodeLabelContainer(nodeContext, outerNodeLabelsOverhang, PortSide.SOUTH);
        placeVerticalOuterNodeLabelContainer(nodeContext, outerNodeLabelsOverhang, PortSide.EAST);
        placeVerticalOuterNodeLabelContainer(nodeContext, outerNodeLabelsOverhang, PortSide.WEST);
    }

    /**
     * Places a horizontal outer node label container on the given side.
     */
    private static void placeHorizontalOuterNodeLabelContainer(final NodeContext nodeContext,
            final boolean outerNodeLabelsOverhang, final PortSide portSide) {
        
        KVector nodeSize = nodeContext.nodeSize;
        StripContainerCell nodeLabelContainer = nodeContext.outsideNodeLabelContainers.get(portSide);
        ElkRectangle nodeLabelContainerRect = nodeLabelContainer.getCellRectangle();
        
        // Set the container's width and height to its minimum width and height
        nodeLabelContainerRect.width = nodeLabelContainer.getMinimumWidth();
        nodeLabelContainerRect.height = nodeLabelContainer.getMinimumHeight();
        
        // The container must be at least as wide as the node is
        nodeLabelContainerRect.width = Math.max(nodeLabelContainerRect.width, nodeSize.x);
        
        // If node labels are not allowed to overhang and if they would do so right now, make the container smaller
        if (nodeLabelContainerRect.width > nodeSize.x && !outerNodeLabelsOverhang) {
            nodeLabelContainerRect.width = nodeSize.x;
        }
        
        // Container's x coordinate
        nodeLabelContainerRect.x = -(nodeLabelContainerRect.width - nodeSize.x) / 2;
        
        // Container's y coordinate depends on whether we place the thing on the northern or southern side
        switch (portSide) {
        case NORTH:
            nodeLabelContainerRect.y = -nodeLabelContainerRect.height;
            break;
            
        case SOUTH:
            nodeLabelContainerRect.y = nodeSize.y;
            break;
        }
        
        // Layout the container's children
        nodeLabelContainer.layoutChildrenHorizontally();
        nodeLabelContainer.layoutChildrenVertically();
    }


    /**
     * Places a vertical outer node label container on the given side.
     */
    private static void placeVerticalOuterNodeLabelContainer(final NodeContext nodeContext,
            final boolean outerNodeLabelsOverhang, final PortSide portSide) {
        
        KVector nodeSize = nodeContext.nodeSize;
        StripContainerCell nodeLabelContainer = nodeContext.outsideNodeLabelContainers.get(portSide);
        ElkRectangle nodeLabelContainerRect = nodeLabelContainer.getCellRectangle();
        
        // Set the container's width and height to its minimum width and height
        nodeLabelContainerRect.width = nodeLabelContainer.getMinimumWidth();
        nodeLabelContainerRect.height = nodeLabelContainer.getMinimumHeight();
        
        // The container must be at least as high as the node is
        nodeLabelContainerRect.height = Math.max(nodeLabelContainerRect.height, nodeSize.y);
        
        // If node labels are not allowed to overhang and if they would do so right now, make the container smaller
        if (nodeLabelContainerRect.height > nodeSize.y && !outerNodeLabelsOverhang) {
            nodeLabelContainerRect.height = nodeSize.y;
        }
        
        // Container's y coordinate
        nodeLabelContainerRect.y = -(nodeLabelContainerRect.height - nodeSize.y) / 2;
        
        // Container's x coordinate depends on whether we place the thing on the eastern or western side
        switch (portSide) {
        case WEST:
            nodeLabelContainerRect.x = -nodeLabelContainerRect.width;
            break;
            
        case EAST:
            nodeLabelContainerRect.x = nodeSize.x;
            break;
        }
        
        // Layout the container's children
        nodeLabelContainer.layoutChildrenHorizontally();
        nodeLabelContainer.layoutChildrenVertically();
    }
    
}
