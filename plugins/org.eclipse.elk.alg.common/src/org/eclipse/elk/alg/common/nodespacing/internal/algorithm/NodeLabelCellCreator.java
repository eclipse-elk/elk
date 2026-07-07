/*******************************************************************************
 * Copyright (c) 2017, 2019 Kiel University and others.
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
package org.eclipse.elk.alg.common.nodespacing.internal.algorithm;

import java.util.Set;

import org.eclipse.elk.alg.common.nodespacing.cellsystem.ContainerArea;
import org.eclipse.elk.alg.common.nodespacing.cellsystem.GridContainerCell;
import org.eclipse.elk.alg.common.nodespacing.cellsystem.LabelCell;
import org.eclipse.elk.alg.common.nodespacing.cellsystem.StripContainerCell;
import org.eclipse.elk.alg.common.nodespacing.cellsystem.StripContainerCell.Strip;
import org.eclipse.elk.alg.common.nodespacing.internal.NodeContext;
import org.eclipse.elk.alg.common.nodespacing.internal.NodeLabelLocation;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.options.NodeLabelPlacement;
import org.eclipse.elk.core.options.PortSide;
import org.eclipse.elk.core.options.SizeOptions;
import org.eclipse.elk.core.util.adapters.GraphAdapters.LabelAdapter;

/**
 * Knows how to take all of a node's labels and create the appropriate grid cells.
 */
public final class NodeLabelCellCreator {
    
    /**
     * No instances required.
     */
    private NodeLabelCellCreator() {
        
    }
    
    
    /**
     * Iterates over all of the node's labels and creates all required cell containers and label cells.
     * 
     * @param nodeContext
     *            the node context to work with.
     * @param onlyInside
     *            {@code true} if only inside node labels should be assigned to their respective label cells. Used for
     *            computing the space required for inside node label cells, but nothing else.
     * @param horizontalLayoutMode 
     *            Whether we operate in horizontal or vertical layout mode.
     */
    public static void createNodeLabelCells(final NodeContext nodeContext, final boolean onlyInside,
            final boolean horizontalLayoutMode) {
        // Make sure all the relevant containers exist
        createNodeLabelCellContainers(nodeContext, onlyInside);
        
        // Handle each of the node's labels
        nodeContext.node.getLabels().forEach(label -> handleNodeLabel(nodeContext, label, onlyInside, 
                horizontalLayoutMode));
    }
    
    /**
     * Handles the given node label by adding it to the corresponding node label cell. If such a cell doesn't exist
     * yet, it is created and added to the correct container cell.
     */
    private static void handleNodeLabel(final NodeContext nodeContext, final LabelAdapter<?> label,
            final boolean onlyInside, final boolean horizontalLayoutMode) {
        
        // Find the effective label location
        Set<NodeLabelPlacement> labelPlacement = label.hasProperty(CoreOptions.NODE_LABELS_PLACEMENT)
                ? label.getProperty(CoreOptions.NODE_LABELS_PLACEMENT)
                : nodeContext.nodeLabelPlacement;
        NodeLabelLocation labelLocation = NodeLabelLocation.fromNodeLabelPlacement(labelPlacement);
        
        // If the label has its location fixed, we will ignore it
        if (labelLocation == NodeLabelLocation.UNDEFINED) {
            return;
        }
        
        // If the label's location is on the node's outside but we only want inside node labels, we will ignore it
        if (onlyInside && !labelLocation.isInsideLocation()) {
            return;
        }
        
        retrieveNodeLabelCell(nodeContext, labelLocation, horizontalLayoutMode).addLabel(label);
    }
    

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Cell Creation and Retrieval

    /**
     * Creates all node label containers.
     */
    private static void createNodeLabelCellContainers(final NodeContext nodeContext, final boolean onlyInside) {
        boolean symmetry = !nodeContext.sizeOptions.contains(SizeOptions.ASYMMETRICAL);
        boolean tabularNodeLabels = nodeContext.sizeOptions.contains(SizeOptions.FORCE_TABULAR_NODE_LABELS);
        
        // Inside container
        nodeContext.insideNodeLabelContainer = new GridContainerCell(
                tabularNodeLabels, symmetry, nodeContext.labelCellSpacing);
        
        if (nodeContext.nodeLabelsPadding != null) {
            nodeContext.insideNodeLabelContainer.getPadding().copy(nodeContext.nodeLabelsPadding);
        }
        nodeContext.nodeContainerMiddleRow.setCell(ContainerArea.CENTER, nodeContext.insideNodeLabelContainer);
        
        // Outside containers, if requested
        if (!onlyInside) {
            StripContainerCell northContainer = new StripContainerCell(
                    Strip.HORIZONTAL, symmetry, nodeContext.labelCellSpacing);
            northContainer.getPadding().bottom = nodeContext.nodeLabelSpacing;
            nodeContext.outsideNodeLabelContainers.put(PortSide.NORTH, northContainer);
            
            StripContainerCell southContainer = new StripContainerCell(
                    Strip.HORIZONTAL, symmetry, nodeContext.labelCellSpacing);
            southContainer.getPadding().top = nodeContext.nodeLabelSpacing;
            nodeContext.outsideNodeLabelContainers.put(PortSide.SOUTH, southContainer);
            
            StripContainerCell westContainer = new StripContainerCell(
                    Strip.VERTICAL, symmetry, nodeContext.labelCellSpacing);
            westContainer.getPadding().right = nodeContext.nodeLabelSpacing;
            nodeContext.outsideNodeLabelContainers.put(PortSide.WEST, westContainer);
            
            StripContainerCell eastContainer = new StripContainerCell(
                    Strip.VERTICAL, symmetry, nodeContext.labelCellSpacing);
            eastContainer.getPadding().left = nodeContext.nodeLabelSpacing;
            nodeContext.outsideNodeLabelContainers.put(PortSide.EAST, eastContainer);
        }
    }
    
    /**
     * Retrieves the node label cell for the given location. If it doesn't exist yet, it is created.
     */
    private static LabelCell retrieveNodeLabelCell(final NodeContext nodeContext,
            final NodeLabelLocation nodeLabelLocation, final boolean horizontalLayoutMode) {
        
        LabelCell nodeLabelCell = nodeContext.nodeLabelCells.get(nodeLabelLocation);
        
        if (nodeLabelCell == null) {
            // The node label cell doesn't exist yet, so create one and add it to the relevant container
            nodeLabelCell = new LabelCell(nodeContext.labelLabelSpacing, nodeLabelLocation, horizontalLayoutMode);
            nodeContext.nodeLabelCells.put(nodeLabelLocation, nodeLabelCell);
            
            // Find the correct container and add the cell to it
            if (nodeLabelLocation.isInsideLocation()) {
                nodeContext.insideNodeLabelContainer.setCell(
                        nodeLabelLocation.getContainerRow(),
                        nodeLabelLocation.getContainerColumn(),
                        nodeLabelCell);
            } else {
                PortSide outsideSide = nodeLabelLocation.getOutsideSide();
                StripContainerCell containerCell = nodeContext.outsideNodeLabelContainers.get(outsideSide);
                
                switch (outsideSide) {
                case NORTH:
                case SOUTH:
                    nodeLabelCell.setContributesToMinimumHeight(true);
                    containerCell.setCell(nodeLabelLocation.getContainerColumn(), nodeLabelCell);
                    break;
                    
                case WEST:
                case EAST:
                    nodeLabelCell.setContributesToMinimumWidth(true);
                    containerCell.setCell(nodeLabelLocation.getContainerRow(), nodeLabelCell);
                    break;
                }
            }
        }
        
        return nodeLabelCell;
    }
    
}
