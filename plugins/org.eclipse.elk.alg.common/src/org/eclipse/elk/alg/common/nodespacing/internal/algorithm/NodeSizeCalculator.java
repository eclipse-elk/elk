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

import org.eclipse.elk.alg.common.nodespacing.internal.NodeContext;
import org.eclipse.elk.core.math.ElkRectangle;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.options.PortConstraints;
import org.eclipse.elk.core.options.PortSide;
import org.eclipse.elk.core.options.SizeConstraint;
import org.eclipse.elk.core.options.SizeOptions;

/**
 * Configures the cell system according to the node size constraints and determines the ultimate node size.
 */
public final class NodeSizeCalculator {
    
    /**
     * No instance required.
     */
    private NodeSizeCalculator() {
        
    }
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Node Width
    
    /**
     * Sets the node's width according to the active node size constraints. Also sets that width on the cell system and
     * tells it to compute a horizontal layout.
     */
    public static void setNodeWidth(final NodeContext nodeContext) {
        KVector nodeSize = nodeContext.nodeSize;
        double width;
        
        if (NodeLabelAndSizeUtilities.areSizeConstraintsFixed(nodeContext)) {
            // Simply use the node's current width
            width = nodeSize.x;
        } else {
            // Ask the cell system how wide it would like to be or take the node's width if it has already been set to a
            // greater value
            if (nodeContext.topdownLayout) {
                width = Math.max(nodeSize.x, nodeContext.nodeContainer.getMinimumWidth());
            } else {
                width = nodeContext.nodeContainer.getMinimumWidth();
            }
            
            // If we include node labels and outside node labels are not to overhang, we need to include those as well
            if (nodeContext.sizeConstraints.contains(SizeConstraint.NODE_LABELS)
                    && !nodeContext.sizeOptions.contains(SizeOptions.OUTSIDE_NODE_LABELS_OVERHANG)) {
                
                width = Math.max(width, nodeContext.outsideNodeLabelContainers.get(PortSide.NORTH).getMinimumWidth());
                width = Math.max(width, nodeContext.outsideNodeLabelContainers.get(PortSide.SOUTH).getMinimumWidth());
            }
            
            // The node might have a minimum size set...
            KVector minNodeSize = NodeLabelAndSizeUtilities.getMinimumNodeSize(nodeContext);
            if (minNodeSize != null) {
                width = Math.max(width, minNodeSize.x);
            }
        }
        
        // Set the node's width
        if (nodeContext.node.getGraph().getProperty(CoreOptions.NODE_SIZE_FIXED_GRAPH_SIZE)) {
            nodeSize.x = Math.max(nodeSize.x, width);
        } else {
            nodeSize.x = width;
        }
        
        // Set the cell system's width and tell it to compute horizontal coordinates and widths
        ElkRectangle nodeCellRectangle = nodeContext.nodeContainer.getCellRectangle();
        nodeCellRectangle.x = 0;
        nodeCellRectangle.width = width;
        
        nodeContext.nodeContainer.layoutChildrenHorizontally();
    }
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Node Height
    
    /**
     * Sets the node's height according to the active node size constraints. Also sets that height on the cell system
     * and tells it to compute a horizontal layout.
     */
    public static void setNodeHeight(final NodeContext nodeContext) {
        KVector nodeSize = nodeContext.nodeSize;
        double height;
        
        if (NodeLabelAndSizeUtilities.areSizeConstraintsFixed(nodeContext)) {
            // Simply use the node's current height
            height = nodeSize.y;
        } else {
            // Ask the cell system how high it would like to be or take the node's height if it has already been set to 
            // a greater value
            if (nodeContext.topdownLayout) {
                height = Math.max(nodeSize.y, nodeContext.nodeContainer.getMinimumHeight());
            } else {
                height = nodeContext.nodeContainer.getMinimumHeight();
            }
            
            // If we include node labels and outside node labels are not to overhang, we need to include those as well
            if (nodeContext.sizeConstraints.contains(SizeConstraint.NODE_LABELS)
                    && !nodeContext.sizeOptions.contains(SizeOptions.OUTSIDE_NODE_LABELS_OVERHANG)) {
                
                height = Math.max(height, nodeContext.outsideNodeLabelContainers.get(PortSide.EAST).getMinimumHeight());
                height = Math.max(height, nodeContext.outsideNodeLabelContainers.get(PortSide.WEST).getMinimumHeight());
            }
            
            // The node might have a minimum size set...
            KVector minNodeSize = NodeLabelAndSizeUtilities.getMinimumNodeSize(nodeContext);
            if (minNodeSize != null) {
                height = Math.max(height, minNodeSize.y);
            }
            
            // If size constraints are set to include ports, but port constraints are FIXED_POS or FIXED_RATIO, we need
            // to manually apply the height required to place eastern and western ports because those heights don't
            // come out of the cell system
            if (nodeContext.sizeConstraints.contains(SizeConstraint.PORTS)) {
                if (nodeContext.portConstraints == PortConstraints.FIXED_RATIO
                        || nodeContext.portConstraints == PortConstraints.FIXED_POS) {
                    
                    height = Math.max(height, nodeContext.insidePortLabelCells.get(PortSide.EAST).getMinimumHeight());
                    height = Math.max(height, nodeContext.insidePortLabelCells.get(PortSide.WEST).getMinimumHeight());
                }
            }
        }
        
        // Set the node's height
        if (nodeContext.node.getGraph().getProperty(CoreOptions.NODE_SIZE_FIXED_GRAPH_SIZE)) {
            nodeSize.y = Math.max(nodeSize.y, height);
        } else {
            nodeSize.y = height;
        }
        
        // Set the cell system's height and tell it to compute vertical coordinates and heights
        ElkRectangle nodeCellRectangle = nodeContext.nodeContainer.getCellRectangle();
        nodeCellRectangle.y = 0;
        nodeCellRectangle.height = height;
        
        nodeContext.nodeContainer.layoutChildrenVertically();
    }
}
