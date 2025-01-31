/*******************************************************************************
 * Copyright (c) 2017, 2020 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.alg.common.nodespacing.internal.algorithm;

import java.util.Collection;

import org.eclipse.elk.alg.common.nodespacing.cellsystem.AtomicCell;
import org.eclipse.elk.alg.common.nodespacing.cellsystem.HorizontalLabelAlignment;
import org.eclipse.elk.alg.common.nodespacing.cellsystem.LabelCell;
import org.eclipse.elk.alg.common.nodespacing.cellsystem.VerticalLabelAlignment;
import org.eclipse.elk.alg.common.nodespacing.internal.NodeContext;
import org.eclipse.elk.alg.common.nodespacing.internal.PortContext;
import org.eclipse.elk.alg.common.overlaps.RectangleStripOverlapRemover;
import org.eclipse.elk.alg.common.overlaps.RectangleStripOverlapRemover.OverlapRemovalDirection;
import org.eclipse.elk.core.math.ElkMath;
import org.eclipse.elk.core.math.ElkRectangle;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.options.PortConstraints;
import org.eclipse.elk.core.options.PortLabelPlacement;
import org.eclipse.elk.core.options.PortSide;
import org.eclipse.elk.core.options.SizeConstraint;

/**
 * Knows how to place port labels.
 */
public final class PortLabelPlacementCalculator {

    /**
     * No instance required.
     */
    private PortLabelPlacementCalculator() {
        
    }
    
    
    /**
     * Places port labels for northern and southern ports. If port labels are placed on the inside, the height required
     * for the placement is set as the height of the content area of northern and southern inside port label cells.
     */
    public static void placeHorizontalPortLabels(final NodeContext nodeContext) {
        placePortLabels(nodeContext, PortSide.NORTH);
        placePortLabels(nodeContext, PortSide.SOUTH);
    }
    
    /**
     * Places port labels for eastern and western ports.
     */
    public static void placeVerticalPortLabels(final NodeContext nodeContext) {
        placePortLabels(nodeContext, PortSide.EAST);
        placePortLabels(nodeContext, PortSide.WEST);
    }
    
    /**
     * Places port labels for ports on the given side.
     */
    private static void placePortLabels(final NodeContext nodeContext, final PortSide portSide) {
        // If port labels were not taken into account when calculating the node size or if port placement was set to
        // fixed positions, we don't have an arbitrary amount of freedom to place our labels
        boolean constrainedPlacement = !nodeContext.sizeConstraints.contains(SizeConstraint.PORT_LABELS)
                || nodeContext.portConstraints == PortConstraints.FIXED_POS;

        if (nodeContext.portLabelsPlacement.contains(PortLabelPlacement.INSIDE)) {
            if (constrainedPlacement) {
                constrainedInsidePortLabelPlacement(nodeContext, portSide);
            } else {
                simpleInsidePortLabelPlacement(nodeContext, portSide);
            }
        } else if (nodeContext.portLabelsPlacement.contains(PortLabelPlacement.OUTSIDE)) {
            if (constrainedPlacement) {
                constrainedOutsidePortLabelPlacement(nodeContext, portSide);
            } else {
                simpleOutsidePortLabelPlacement(nodeContext, portSide);
            }
        }
    }
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Simple Inside Port Labels
    
    /**
     * Place the port label cells on the node's insides.
     */
    private static void simpleInsidePortLabelPlacement(final NodeContext nodeContext, final PortSide portSide) {
        // For northern and southern port labels, we need to set the inside port label cell's client height later
        double insideNorthOrSouthPortLabelAreaHeight = 0;
        
        // Some spacings we may need later
        double labelBorderOffset = portLabelBorderOffsetForPortSide(nodeContext, portSide);
        double portLabelSpacingHorizontal = nodeContext.portLabelSpacingHorizontal;
        double portLabelSpacingVertical = nodeContext.portLabelSpacingVertical;
        
        for (PortContext portContext : nodeContext.portContexts.get(portSide)) {
            // If the port doesn't have labels, skip
            if (portContext.portLabelCell == null || !portContext.portLabelCell.hasLabels()) {
                continue;
            }
            
            // Retrieve information about the port itself
            KVector portSize = portContext.port.getSize();
            double portBorderOffset = portContext.port.hasProperty(CoreOptions.PORT_BORDER_OFFSET)
                    ? portContext.port.getProperty(CoreOptions.PORT_BORDER_OFFSET)
                    : 0;
            
            // Retrieve the label cell and its rectangle and set the rectangle's size (we will use the rectangle to
            // place the cell relative to the port below)
            LabelCell portLabelCell = portContext.portLabelCell;
            ElkRectangle portLabelCellRect = portLabelCell.getCellRectangle();
            portLabelCellRect.width = portLabelCell.getMinimumWidth();
            portLabelCellRect.height = portLabelCell.getMinimumHeight();
            
            // Calculate the position of the port's label cell. If the node is a compound node, we make an effort to
            // place port labels such that edges won't cross them
            switch (portSide) {
            case NORTH:
                portLabelCellRect.x = portContext.labelsNextToPort
                        ? (portSize.x - portLabelCellRect.width) / 2
                        : portSize.x + portLabelSpacingHorizontal;
                portLabelCellRect.y = portSize.y + portBorderOffset + labelBorderOffset;
                portLabelCell.setHorizontalAlignment(HorizontalLabelAlignment.CENTER);
                portLabelCell.setVerticalAlignment(VerticalLabelAlignment.TOP);
                break;
                
            case SOUTH:
                portLabelCellRect.x = portContext.labelsNextToPort
                        ? (portSize.x - portLabelCellRect.width) / 2
                        : portSize.x + portLabelSpacingHorizontal;
                portLabelCellRect.y = -portBorderOffset - labelBorderOffset - portLabelCellRect.height;
                portLabelCell.setHorizontalAlignment(HorizontalLabelAlignment.CENTER);
                portLabelCell.setVerticalAlignment(VerticalLabelAlignment.BOTTOM);
                break;
                
            case EAST:
                portLabelCellRect.x = -portBorderOffset - labelBorderOffset - portLabelCellRect.width;
                if (portContext.labelsNextToPort) {
                    double labelHeight = nodeContext.portLabelsTreatAsGroup
                            ? portLabelCellRect.height
                            : portLabelCell.getLabels().get(0).getSize().y;
                    portLabelCellRect.y = (portSize.y - labelHeight) / 2;
                } else {
                    portLabelCellRect.y = portSize.y + portLabelSpacingVertical;
                }
                portLabelCell.setHorizontalAlignment(HorizontalLabelAlignment.RIGHT);
                portLabelCell.setVerticalAlignment(VerticalLabelAlignment.CENTER);
                break;
                
            case WEST:
                portLabelCellRect.x = portSize.x + portBorderOffset + labelBorderOffset;
                if (portContext.labelsNextToPort) {
                    double labelHeight = nodeContext.portLabelsTreatAsGroup
                            ? portLabelCellRect.height
                            : portLabelCell.getLabels().get(0).getSize().y;
                    portLabelCellRect.y = (portSize.y - labelHeight) / 2;
                } else {
                    portLabelCellRect.y = portSize.y + portLabelSpacingVertical;
                }
                portLabelCell.setHorizontalAlignment(HorizontalLabelAlignment.LEFT);
                portLabelCell.setVerticalAlignment(VerticalLabelAlignment.CENTER);
                break;
            }
            
            // If we have a north or south port, update our port label area height
            if (portSide == PortSide.NORTH || portSide == PortSide.SOUTH) {
                insideNorthOrSouthPortLabelAreaHeight = Math.max(
                        insideNorthOrSouthPortLabelAreaHeight,
                        portLabelCellRect.height);
            }
        }
        
        // If we have a northern or southern label area height, apply it
        if (insideNorthOrSouthPortLabelAreaHeight > 0) {
            nodeContext.insidePortLabelCells.get(portSide).getMinimumContentAreaSize().y =
                    insideNorthOrSouthPortLabelAreaHeight;
        }
    }
    
    /**
     * Returns the amount of space port labels are offset from the node border on the given port side. This takes the
     * port border padding into account as well as the port-label spacing.
     */
    private static double portLabelBorderOffsetForPortSide(final NodeContext nodeContext, final PortSide portSide) {
        switch (portSide) {
        case NORTH:
            return nodeContext.nodeContainer.getPadding().top + nodeContext.portLabelSpacingVertical;

        case SOUTH:
            return nodeContext.nodeContainer.getPadding().bottom + nodeContext.portLabelSpacingVertical;

        case EAST:
            return nodeContext.nodeContainer.getPadding().right + nodeContext.portLabelSpacingHorizontal;

        case WEST:
            return nodeContext.nodeContainer.getPadding().left + nodeContext.portLabelSpacingHorizontal;
        
        default:
            assert false;
            return 0;
        }
    }
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Constrained Inside Port Labels

    /**
     * Place the port label cells outside of the node in the knowledge that there might not be enough space to place
     * them without overlaps.
     */
    private static void constrainedInsidePortLabelPlacement(final NodeContext nodeContext, final PortSide portSide) {
        Collection<PortContext> portContexts = nodeContext.portContexts.get(portSide);
    
        // If it's neither the northern nor the southern side, simply revert to simple port label placement
        if (portSide == PortSide.EAST || portSide == PortSide.WEST) {
            simpleInsidePortLabelPlacement(nodeContext, portSide);
            return;
        }
        
        // Prepare things
        OverlapRemovalDirection overlapRemovalDirection = portSide == PortSide.NORTH
                ? OverlapRemovalDirection.DOWN
                : OverlapRemovalDirection.UP;
        VerticalLabelAlignment verticalLabelAlignment = portSide == PortSide.NORTH
                ? VerticalLabelAlignment.TOP
                : VerticalLabelAlignment.BOTTOM;
        
        // To keep labels from extending over the content area of the inside port label container, we need to know
        // where its content area's left and right boundaries are. We also make sure to always keep a bit of space to
        // the node border
        AtomicCell insidePortLabelContainer = nodeContext.insidePortLabelCells.get(portSide);
        ElkRectangle labelContainerRect = insidePortLabelContainer.getCellRectangle();
        double leftBorder = labelContainerRect.x + ElkMath.maxd(
                insidePortLabelContainer.getPadding().left,
                nodeContext.surroundingPortMargins.left,
                nodeContext.nodeLabelSpacing);
        double rightBorder = labelContainerRect.x + labelContainerRect.width - ElkMath.maxd(
                insidePortLabelContainer.getPadding().right,
                nodeContext.surroundingPortMargins.right,
                nodeContext.nodeLabelSpacing);
        
        // Obtain a rectangle strip overlap remover, which will actually do most of the work
        RectangleStripOverlapRemover overlapRemover = RectangleStripOverlapRemover
                .createForDirection(overlapRemovalDirection)
                .withGap(nodeContext.portLabelSpacingHorizontal, nodeContext.portLabelSpacingVertical);
        
        // Iterate over our ports and add rectangles to the overlap remover. Also, calculate the start coordinate
        double startCoordinate = portSide == PortSide.NORTH
                ? Double.NEGATIVE_INFINITY
                : Double.POSITIVE_INFINITY;
        
        for (PortContext portContext : portContexts) {
            if (portContext.portLabelCell == null || !portContext.portLabelCell.hasLabels()) {
                continue;
            }
            
            KVector portSize = portContext.port.getSize();
            KVector portPosition = portContext.portPosition;
            LabelCell portLabelCell = portContext.portLabelCell;
            ElkRectangle portLabelCellRect = portLabelCell.getCellRectangle();
            
            // Setup the less interesting cell properties
            portLabelCellRect.width = portLabelCell.getMinimumWidth();
            portLabelCellRect.height = portLabelCell.getMinimumHeight();
            
            portLabelCell.setVerticalAlignment(verticalLabelAlignment);
            portLabelCell.setHorizontalAlignment(HorizontalLabelAlignment.RIGHT);
            
            // Center the label, but make sure it doesn't hang over the node boundaries
            centerPortLabel(portLabelCellRect, portPosition, portSize, leftBorder, rightBorder);
            
            // Add the rectangle to the overlap remover
            overlapRemover.addRectangle(portLabelCellRect);
            
            // Update start coordinate
            startCoordinate = portSide == PortSide.NORTH
                    ? Math.max(startCoordinate, portPosition.y + portContext.port.getSize().y)
                    : Math.min(startCoordinate, portPosition.y);
        }
        
        // The start coordinate needs to be offset by the port-label space
        startCoordinate += portSide == PortSide.NORTH
                ? nodeContext.portLabelSpacingVertical
                : -nodeContext.portLabelSpacingVertical;
        
        // Invoke the overlap remover
        double stripHeight = overlapRemover
            .withStartCoordinate(startCoordinate)
            .removeOverlaps();
        
        if (stripHeight > 0) {
            nodeContext.insidePortLabelCells.get(portSide).getMinimumContentAreaSize().y = stripHeight;
        }
        
        // We need to update the label cell's coordinates to be relative to the ports
        for (PortContext portContext : portContexts) {
            if (portContext.portLabelCell == null || !portContext.portLabelCell.hasLabels()) {
                continue;
            }
            
            ElkRectangle portLabelCellRect = portContext.portLabelCell.getCellRectangle();
            
            // Setup the label cell's cell rectangle
            portLabelCellRect.x -= portContext.portPosition.x;
            portLabelCellRect.y -= portContext.portPosition.y;
        }
    }
    
    
    /**
     * Centers the given label under its port, but makes an effort to keep it from hanging over the given minimum and
     * maximum coordinates. The label position is absolute, not relative to the port.
     */
    private static void centerPortLabel(final ElkRectangle portLabelCellRect, final KVector portPosition,
            final KVector portSize, final double minX, final double maxX) {
        
        // Center the label
        portLabelCellRect.x = portPosition.x - (portLabelCellRect.width - portSize.x) / 2;
        
        // Make sure that the label won't slide past the port
        double actualMinX = Math.min(minX, portPosition.x);
        double actualMaxX = Math.max(maxX, portPosition.x + portSize.x);
        
        // Make sure that the label stays inside the boundaries, but only correct in one of the two possible directions
        if (portLabelCellRect.x < actualMinX) {
            portLabelCellRect.x = actualMinX;
        } else if (portLabelCellRect.x + portLabelCellRect.width > actualMaxX) {
            portLabelCellRect.x = actualMaxX - portLabelCellRect.width;
        }
    }


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Simple Outside Port Labels

    /**
     * Place the port label cells outside of the node.
     */
    private static void simpleOutsidePortLabelPlacement(final NodeContext nodeContext, final PortSide portSide) {
        Collection<PortContext> portContexts = nodeContext.portContexts.get(portSide);
        
        // If there are only two ports on a side, we place the first port's label on its other side to make it
        // especially clear which port it belongs to. The same applies if the user requested space-efficient mode
        boolean placeFirstPortDifferently = NodeLabelAndSizeUtilities.isFirstOutsidePortLabelPlacedDifferently(
                nodeContext, portSide);
        
        boolean alwaysAbove = nodeContext.portLabelsPlacement.contains(PortLabelPlacement.ALWAYS_OTHER_SAME_SIDE);
        
        for (PortContext portContext : portContexts) {
            // If the port doesn't have labels, skip
            if (portContext.portLabelCell == null || !portContext.portLabelCell.hasLabels()) {
                continue;
            }

            // Retrieve information about the port itself
            KVector portSize = portContext.port.getSize();
            
            // Retrieve the label cell and its rectangle and set the rectangle's size (we will use the rectangle to
            // place the cell relative to the port below)
            LabelCell portLabelCell = portContext.portLabelCell;
            ElkRectangle portLabelCellRect = portLabelCell.getCellRectangle();
            portLabelCellRect.width = portLabelCell.getMinimumWidth();
            portLabelCellRect.height = portLabelCell.getMinimumHeight();
            
            // Calculate the position of the port's label space
            switch (portSide) {
            case NORTH:
                if (portContext.labelsNextToPort) {
                    portLabelCellRect.x = (portSize.x - portLabelCellRect.width) / 2;
                    portLabelCell.setHorizontalAlignment(HorizontalLabelAlignment.CENTER);
                } else if (placeFirstPortDifferently || alwaysAbove) {
                    portLabelCellRect.x = -portLabelCellRect.width - nodeContext.portLabelSpacingHorizontal;
                    portLabelCell.setHorizontalAlignment(HorizontalLabelAlignment.RIGHT);
                } else {
                    portLabelCellRect.x = portSize.x + nodeContext.portLabelSpacingHorizontal;
                    portLabelCell.setHorizontalAlignment(HorizontalLabelAlignment.LEFT);
                }
                portLabelCellRect.y = -portLabelCellRect.height - nodeContext.portLabelSpacingVertical;
                portLabelCell.setVerticalAlignment(VerticalLabelAlignment.BOTTOM);
                break;
                
            case SOUTH:
                if (portContext.labelsNextToPort) {
                    portLabelCellRect.x = (portSize.x - portLabelCellRect.width) / 2;
                    portLabelCell.setHorizontalAlignment(HorizontalLabelAlignment.CENTER);
                } else if (placeFirstPortDifferently || alwaysAbove) {
                    portLabelCellRect.x = -portLabelCellRect.width - nodeContext.portLabelSpacingHorizontal;
                    portLabelCell.setHorizontalAlignment(HorizontalLabelAlignment.RIGHT);
                } else {
                    portLabelCellRect.x = portSize.x + nodeContext.portLabelSpacingHorizontal;
                    portLabelCell.setHorizontalAlignment(HorizontalLabelAlignment.LEFT);
                }
                portLabelCellRect.y = portSize.y + nodeContext.portLabelSpacingVertical;
                portLabelCell.setVerticalAlignment(VerticalLabelAlignment.TOP);
                break;
                
            case EAST:
                if (portContext.labelsNextToPort) {
                    double labelHeight = nodeContext.portLabelsTreatAsGroup
                            ? portLabelCellRect.height
                            : portLabelCell.getLabels().get(0).getSize().y;
                    portLabelCellRect.y = (portSize.y - labelHeight) / 2;
                    portLabelCell.setVerticalAlignment(VerticalLabelAlignment.CENTER);
                } else if (placeFirstPortDifferently || alwaysAbove) {
                    portLabelCellRect.y = -portLabelCellRect.height - nodeContext.portLabelSpacingVertical;
                    portLabelCell.setVerticalAlignment(VerticalLabelAlignment.BOTTOM);
                } else {
                    portLabelCellRect.y = portSize.y + nodeContext.portLabelSpacingVertical;
                    portLabelCell.setVerticalAlignment(VerticalLabelAlignment.TOP);
                }
                portLabelCellRect.x = portSize.x + nodeContext.portLabelSpacingHorizontal;
                portLabelCell.setHorizontalAlignment(HorizontalLabelAlignment.LEFT);
                break;
                
            case WEST:
                if (portContext.labelsNextToPort) {
                    double labelHeight = nodeContext.portLabelsTreatAsGroup
                            ? portLabelCellRect.height
                            : portLabelCell.getLabels().get(0).getSize().y;
                    portLabelCellRect.y = (portSize.y - labelHeight) / 2;
                    portLabelCell.setVerticalAlignment(VerticalLabelAlignment.CENTER);
                } else if (placeFirstPortDifferently || alwaysAbove) {
                    portLabelCellRect.y = -portLabelCellRect.height - nodeContext.portLabelSpacingVertical;
                    portLabelCell.setVerticalAlignment(VerticalLabelAlignment.BOTTOM);
                } else {
                    portLabelCellRect.y = portSize.y + nodeContext.portLabelSpacingVertical;
                    portLabelCell.setVerticalAlignment(VerticalLabelAlignment.TOP);
                }
                portLabelCellRect.x = -portLabelCellRect.width - nodeContext.portLabelSpacingHorizontal;
                portLabelCell.setHorizontalAlignment(HorizontalLabelAlignment.RIGHT);
                break;
            }
            
            // The next port definitely doesn't have special needs anymore
            placeFirstPortDifferently = false;
        }
    }
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Constrained Outside Port Labels

    /**
     * Place the port label cells outside of the node in the knowledge that there might not be enough space to place
     * them without overlaps.
     */
    private static void constrainedOutsidePortLabelPlacement(final NodeContext nodeContext, final PortSide portSide) {
        Collection<PortContext> portContexts = nodeContext.portContexts.get(portSide);
    
        // If there are at most two ports on this port side, or if it's neither the northern nor the southern side,
        // simply revert to simple port label placement
        if (portContexts.size() <= 2 || portSide == PortSide.EAST || portSide == PortSide.WEST) {
            simpleOutsidePortLabelPlacement(nodeContext, portSide);
            return;
        }
        
        // If space-efficient port labels are active, the leftmost / topmost port's label must be placed to its left /
        // above it
        boolean portWithSpecialNeeds = nodeContext.portLabelsPlacement.contains(PortLabelPlacement.SPACE_EFFICIENT);

        // Prepare things
        OverlapRemovalDirection overlapRemovalDirection = portSide == PortSide.NORTH
                ? OverlapRemovalDirection.UP
                : OverlapRemovalDirection.DOWN;
        VerticalLabelAlignment verticalLabelAlignment = portSide == PortSide.NORTH
                ? VerticalLabelAlignment.BOTTOM
                : VerticalLabelAlignment.TOP;
        
        // Obtain a rectangle strip overlap remover, which will actually do most of the work
        RectangleStripOverlapRemover overlapRemover = RectangleStripOverlapRemover
                .createForDirection(overlapRemovalDirection)
                .withGap(nodeContext.portLabelSpacingVertical, nodeContext.portLabelSpacingHorizontal);
        
        // Iterate over our ports and add rectangles to the overlap remover. Also, calculate the start coordinate
        double startCoordinate = portSide == PortSide.NORTH
                ? Double.POSITIVE_INFINITY
                : Double.NEGATIVE_INFINITY;
        
        for (PortContext portContext : portContexts) {
            if (portContext.portLabelCell == null || !portContext.portLabelCell.hasLabels()) {
                continue;
            }
            
            KVector portSize = portContext.port.getSize();
            KVector portPosition = portContext.portPosition;
            LabelCell portLabelCell = portContext.portLabelCell;
            ElkRectangle portLabelCellRect = portLabelCell.getCellRectangle();
            
            // Setup the label cell's cell rectangle
            portLabelCellRect.width = portLabelCell.getMinimumWidth();
            portLabelCellRect.height = portLabelCell.getMinimumHeight();
            if (portWithSpecialNeeds) {
                portLabelCellRect.x =
                        portPosition.x - portLabelCell.getMinimumWidth() - nodeContext.portLabelSpacingHorizontal;
                portWithSpecialNeeds = false;
            } else {
                portLabelCellRect.x = portPosition.x + portSize.x + nodeContext.portLabelSpacingHorizontal;
            }
            
            portLabelCell.setVerticalAlignment(verticalLabelAlignment);
            portLabelCell.setHorizontalAlignment(HorizontalLabelAlignment.RIGHT);
            
            // Add the rectangle to the overlap remover
            overlapRemover.addRectangle(portLabelCellRect);
            
            // Update start coordinate
            startCoordinate = portSide == PortSide.NORTH
                    ? Math.min(startCoordinate, portPosition.y)
                    : Math.max(startCoordinate, portPosition.y + portContext.port.getSize().y);
        }
        
        // The start coordinate needs to be offset by the port-label space
        startCoordinate += portSide == PortSide.NORTH
                ? -nodeContext.portLabelSpacingVertical
                : nodeContext.portLabelSpacingVertical;
        
        // Invoke the overlap remover
        overlapRemover
            .withStartCoordinate(startCoordinate)
            .removeOverlaps();
        
        // We need to update the label cell's coordinates to be relative to the ports
        for (PortContext portContext : portContexts) {
            if (portContext.portLabelCell == null || !portContext.portLabelCell.hasLabels()) {
                continue;
            }
            
            ElkRectangle portLabelCellRect = portContext.portLabelCell.getCellRectangle();
            
            // Setup the label cell's cell rectangle
            portLabelCellRect.x -= portContext.portPosition.x;
            portLabelCellRect.y -= portContext.portPosition.y;
        }
    }
    
}
