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
package org.eclipse.elk.alg.common.nodespacing.internal.algorithm;

import org.eclipse.elk.alg.common.nodespacing.cellsystem.AtomicCell;
import org.eclipse.elk.alg.common.nodespacing.internal.NodeContext;
import org.eclipse.elk.alg.common.nodespacing.internal.PortContext;
import org.eclipse.elk.core.math.ElkPadding;
import org.eclipse.elk.core.math.ElkRectangle;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.options.PortAlignment;
import org.eclipse.elk.core.options.PortConstraints;
import org.eclipse.elk.core.options.PortSide;
import org.eclipse.elk.core.options.SizeOptions;
import org.eclipse.elk.core.util.adapters.GraphAdapters.PortAdapter;
import org.eclipse.elk.graph.properties.IProperty;
import org.eclipse.elk.graph.properties.Property;

/**
 * Actually places ports.
 */
public final class PortPlacementCalculator {
    
    /**
     * Copy of the {@link de.cau.cs.kieler.klay.layered.properties.InternalProperties#PORT_RATIO_OR_POSITION}
     * option. For further information see the documentation found there. We added this copy here to
     * allow a generic treatment of spacing calculations for graph elements. See the
     * {@link org.eclipse.elk.alg.common.nodespacing} package. [programmatically set]
     * 
     * TODO This should be replaced by a proper property declared in a MELK file.
     */
    public static final IProperty<Double> PORT_RATIO_OR_POSITION = new Property<Double>(
            "portRatioOrPosition", 0.0);
    

    /**
     * No instance required.
     */
    private PortPlacementCalculator() {
        
    }
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Horizontal Port Placement
    
    /**
     * Places horizontal ports.
     */
    public static void placeHorizontalPorts(final NodeContext nodeContext) {
     // How we are going to place the ports depends on their constraints
        switch (nodeContext.portConstraints) {
        case FIXED_POS:
            placeHorizontalFixedPosPorts(nodeContext, PortSide.NORTH);
            placeHorizontalFixedPosPorts(nodeContext, PortSide.SOUTH);
            break;
            
        case FIXED_RATIO:
            placeHorizontalFixedRatioPorts(nodeContext, PortSide.NORTH);
            placeHorizontalFixedRatioPorts(nodeContext, PortSide.SOUTH);
            break;
            
        default:
            placeHorizontalFreePorts(nodeContext, PortSide.NORTH);
            placeHorizontalFreePorts(nodeContext, PortSide.SOUTH);
            break;
        }
    }

    /**
     * Places ports on the northern and southern side for the {@link PortConstraints#FIXED_POS} case.
     */
    private static void placeHorizontalFixedPosPorts(final NodeContext nodeContext, final PortSide portSide) {
        for (PortContext portContext : nodeContext.portContexts.get(portSide)) {
            // The port's x coordinate is already fixed anyway, so simply adjust its y coordinate according to its
            // offset, if any
            portContext.portPosition.y = calculateHorizontalPortYCoordinate(portContext);
        }
    }

    /**
     * Places ports on the northern and southern side for the {@link PortConstraints#FIXED_RATIO} case.
     */
    private static void placeHorizontalFixedRatioPorts(final NodeContext nodeContext, final PortSide portSide) {
        final double nodeWidth = nodeContext.nodeSize.x;
        
        for (PortContext portContext : nodeContext.portContexts.get(portSide)) {
            // The x coordinate is a function of the node's width and the port's position ratio
            portContext.portPosition.x = nodeWidth * portContext.port.getProperty(PORT_RATIO_OR_POSITION);
            portContext.portPosition.y = calculateHorizontalPortYCoordinate(portContext);
        }
    }

    /**
     * Places ports on the northern and southern side for the unconstrained cases.
     */
    private static void placeHorizontalFreePorts(final NodeContext nodeContext, final PortSide portSide) {
        // If there are no ports on the given side, abort
        if (nodeContext.portContexts.get(portSide).isEmpty()) {
            return;
        }
        
        // Retrieve the proper inside port label cell, which will give us hints as to where to place our ports
        AtomicCell insidePortLabelCell = nodeContext.insidePortLabelCells.get(portSide);
        ElkRectangle insidePortLabelCellRectangle = insidePortLabelCell.getCellRectangle();
        ElkPadding insidePortLabelCellPadding = insidePortLabelCell.getPadding();
        
        // Note that we don't have to distinguish any cases here because the port margins already include space
        // required for labels, if such space is to be reserved. Yay!
        PortAlignment portAlignment = nodeContext.getPortAlignment(portSide);
        double availableSpace = insidePortLabelCellRectangle.width - insidePortLabelCellPadding.left
                - insidePortLabelCellPadding.right;
        double calculatedPortPlacementWidth = insidePortLabelCell.getMinimumContentAreaSize().x;
        double currentXPos = insidePortLabelCellRectangle.x + insidePortLabelCellPadding.left;
        double spaceBetweenPorts = nodeContext.portPortSpacing;
        
        // If the port alignment is distributed or justified, but there's only a single port, we change the alignment
        // to center to keep things from looking stupid
        if ((portAlignment == PortAlignment.DISTRIBUTED || portAlignment == PortAlignment.JUSTIFIED)
                && nodeContext.portContexts.get(portSide).size() == 1) {
            
            calculatedPortPlacementWidth = modifiedPortPlacementSize(
                    nodeContext, portAlignment, calculatedPortPlacementWidth);
            portAlignment = PortAlignment.CENTER;
        }
        
        if (availableSpace < calculatedPortPlacementWidth
                && !nodeContext.sizeOptions.contains(SizeOptions.PORTS_OVERHANG)) {

            // There is not enough space available for the ports, but they are not allowed to overhang either. Reduce
            // the space between them to cram them into the available space.
            if (portAlignment == PortAlignment.DISTRIBUTED) {
                spaceBetweenPorts += (availableSpace - calculatedPortPlacementWidth)
                        / (nodeContext.portContexts.get(portSide).size() + 1);
                currentXPos += spaceBetweenPorts;
                
            } else {
                spaceBetweenPorts += (availableSpace - calculatedPortPlacementWidth)
                        / (nodeContext.portContexts.get(portSide).size() - 1);
            }
        } else {
            // We are allowed to overhang. Yay. However, if we use distributed or justified alignment, this is another
            // case where we should fall back to centered alignment
            if (availableSpace < calculatedPortPlacementWidth) {
                calculatedPortPlacementWidth = modifiedPortPlacementSize(
                        nodeContext, portAlignment, calculatedPortPlacementWidth);
                portAlignment = PortAlignment.CENTER;
            }
            
            // Calculate where we need to start placing ports (note that the node size required by the port placement
            // includes left and right surrounding port margins, which changes the formulas a bit from what you'd
            // otherwise expect)
            switch (portAlignment) {
            case BEGIN:
                // There's nothing to do here
                break;
                
            case CENTER:
                currentXPos += (availableSpace - calculatedPortPlacementWidth) / 2;
                break;
                
            case END:
                currentXPos += availableSpace - calculatedPortPlacementWidth;
                break;
            
            case DISTRIBUTED:
                // In this case, if there is not enough space available to place the ports, we are allowed to overhang.
                // We thus need to ensure that we're only ever increasing the port spacing here
                double additionalSpaceBetweenPorts = (availableSpace - calculatedPortPlacementWidth)
                        / (nodeContext.portContexts.get(portSide).size() + 1);
                spaceBetweenPorts += Math.max(0, additionalSpaceBetweenPorts);
                currentXPos += spaceBetweenPorts;
                break;
                
            case JUSTIFIED:
                // In this case, if there is not enough space available to place the ports, we are allowed to overhang.
                // We thus need to ensure that we're only ever increasing the port spacing here
                additionalSpaceBetweenPorts = (availableSpace - calculatedPortPlacementWidth)
                        / (nodeContext.portContexts.get(portSide).size() - 1);
                spaceBetweenPorts += Math.max(0, additionalSpaceBetweenPorts);
                break;
            }
        }
        
        // Iterate over all ports and place them
        for (PortContext portContext : nodeContext.portContexts.get(portSide)) {
            portContext.portPosition.x = currentXPos + portContext.portMargin.left;
            portContext.portPosition.y = calculateHorizontalPortYCoordinate(portContext);
            
            // Update the x coordinate for the next port
            currentXPos += portContext.portMargin.left
                    + portContext.port.getSize().x
                    + portContext.portMargin.right
                    + spaceBetweenPorts;
        }
    }

    /**
     * Returns the y coordinate for the given port.
     */
    private static double calculateHorizontalPortYCoordinate(final PortContext portContext) {
        PortAdapter<?> port = portContext.port;
        
        // The y coordinate is set according to the port's offset, if any
        if (port.hasProperty(CoreOptions.PORT_BORDER_OFFSET)) {
            return port.getSide() == PortSide.NORTH
                    ? -port.getSize().y - port.getProperty(CoreOptions.PORT_BORDER_OFFSET)
                    : port.getProperty(CoreOptions.PORT_BORDER_OFFSET);
        } else {
            return port.getSide() == PortSide.NORTH
                    ? -port.getSize().y
                    : 0;
        }
    }
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Vertical Port Placement
    
    /**
     * Places vertical ports.
     */
    public static void placeVerticalPorts(final NodeContext nodeContext) {
        // How we are going to place the ports depends on their constraints
        switch (nodeContext.portConstraints) {
        case FIXED_POS:
            placeVerticalFixedPosPorts(nodeContext, PortSide.EAST);
            placeVerticalFixedPosPorts(nodeContext, PortSide.WEST);
            break;
            
        case FIXED_RATIO:
            placeVerticalFixedRatioPorts(nodeContext, PortSide.EAST);
            placeVerticalFixedRatioPorts(nodeContext, PortSide.WEST);
            break;
            
        default:
            placeVerticalFreePorts(nodeContext, PortSide.EAST);
            placeVerticalFreePorts(nodeContext, PortSide.WEST);
            break;
        }
    }

    /**
     * Places ports on the eastern and western side for the {@link PortConstraints#FIXED_POS} case.
     */
    private static void placeVerticalFixedPosPorts(final NodeContext nodeContext, final PortSide portSide) {
        double nodeWidth = nodeContext.nodeSize.x;
        
        for (PortContext portContext : nodeContext.portContexts.get(portSide)) {
            // The port's y coordinate is already fixed anyway, so simply adjust its x coordinate according to its
            // offset, if any
            portContext.portPosition.x = calculateVerticalPortXCoordinate(portContext, nodeWidth);
        }
    }

    /**
     * Places ports on the eastern and western side for the {@link PortConstraints#FIXED_RATIO} case.
     */
    private static void placeVerticalFixedRatioPorts(final NodeContext nodeContext, final PortSide portSide) {
        final KVector nodeSize = nodeContext.nodeSize;
        
        for (PortContext portContext : nodeContext.portContexts.get(portSide)) {
            // The y coordinate is a function of the node's width and the port's position ratio
            portContext.portPosition.x = calculateVerticalPortXCoordinate(portContext, nodeSize.x);
            portContext.portPosition.y = nodeSize.y * portContext.port.getProperty(PORT_RATIO_OR_POSITION);
        }
    }

    /**
     * Places ports on the eastern and western side for the unconstrained cases.
     */
    private static void placeVerticalFreePorts(final NodeContext nodeContext, final PortSide portSide) {
        // If there are no ports on the given side, abort
        if (nodeContext.portContexts.get(portSide).isEmpty()) {
            return;
        }
        
        // Retrieve the proper inside port label cell, which will give us hints as to where to place our ports
        AtomicCell insidePortLabelCell = nodeContext.insidePortLabelCells.get(portSide);
        ElkRectangle insidePortLabelCellRectangle = insidePortLabelCell.getCellRectangle();
        ElkPadding insidePortLabelCellPadding = insidePortLabelCell.getPadding();
        
        // Note that we don't have to distinguish any cases here because the port margins already include space
        // required for labels, if such space is to be reserved. Yay!
        PortAlignment portAlignment = nodeContext.getPortAlignment(portSide);
        double availableSpace = insidePortLabelCellRectangle.height - insidePortLabelCellPadding.top
                - insidePortLabelCellPadding.bottom;
        double calculatedPortPlacementHeight = insidePortLabelCell.getMinimumContentAreaSize().y;
        double currentYPos = insidePortLabelCellRectangle.y + insidePortLabelCellPadding.top;
        double spaceBetweenPorts = nodeContext.portPortSpacing;
        double nodeWidth = nodeContext.nodeSize.x;
        
        // If the port alignment is distributed or justified, but there's only a single port, we change the alignment
        // to center to keep things from looking stupid
        if ((portAlignment == PortAlignment.DISTRIBUTED || portAlignment == PortAlignment.JUSTIFIED)
                && nodeContext.portContexts.get(portSide).size() == 1) {
            
            calculatedPortPlacementHeight = modifiedPortPlacementSize(
                    nodeContext, portAlignment, calculatedPortPlacementHeight);
            portAlignment = PortAlignment.CENTER;
        }
        
        if (availableSpace < calculatedPortPlacementHeight
                && !nodeContext.sizeOptions.contains(SizeOptions.PORTS_OVERHANG)) {

            // There is not enough space available for the ports, but they are not allowed to overhang either. Reduce
            // the space between them to cram them into the available space.
            if (portAlignment == PortAlignment.DISTRIBUTED) {
                spaceBetweenPorts += (availableSpace - calculatedPortPlacementHeight)
                        / (nodeContext.portContexts.get(portSide).size() + 1);
                currentYPos += spaceBetweenPorts;
                
            } else {
                spaceBetweenPorts += (availableSpace - calculatedPortPlacementHeight)
                        / (nodeContext.portContexts.get(portSide).size() - 1);
            }
        } else {
            // We are allowed to overhang. Yay. However, if we use distributed or justified alignment, this is another
            // case where we should fall back to centered alignment
            if (availableSpace < calculatedPortPlacementHeight) {
                calculatedPortPlacementHeight = modifiedPortPlacementSize(
                        nodeContext, portAlignment, calculatedPortPlacementHeight);
                portAlignment = PortAlignment.CENTER;
            }
            
            // Calculate where we need to start placing ports (note that the node size required by the port placement
            // includes left and right surrounding port margins, which changes the formulas a bit from what you'd
            // otherwise expect)
            switch (portAlignment) {
            case BEGIN:
                // There's nothing to do here
                break;
                
            case CENTER:
                currentYPos += (availableSpace - calculatedPortPlacementHeight) / 2;
                break;
                
            case END:
                currentYPos += availableSpace - calculatedPortPlacementHeight;
                break;
                
            case DISTRIBUTED:
                // In this case, if there is not enough space available to place the ports, we are allowed to overhang.
                // We thus need to ensure that we're only ever increasing the port spacing here
                double additionalSpaceBetweenPorts = (availableSpace - calculatedPortPlacementHeight)
                        / (nodeContext.portContexts.get(portSide).size() + 1);
                spaceBetweenPorts += Math.max(0, additionalSpaceBetweenPorts);
                currentYPos += spaceBetweenPorts;
                break;
                
            case JUSTIFIED:
                // In this case, if there is not enough space available to place the ports, we are allowed to overhang.
                // We thus need to ensure that we're only ever increasing the port spacing here
                additionalSpaceBetweenPorts = (availableSpace - calculatedPortPlacementHeight)
                        / (nodeContext.portContexts.get(portSide).size() - 1);
                spaceBetweenPorts += Math.max(0, additionalSpaceBetweenPorts);
                break;
            }
        }
        
        // Iterate over all ports and place them
        for (PortContext portContext : nodeContext.portContexts.get(portSide)) {
            portContext.portPosition.x = calculateVerticalPortXCoordinate(portContext, nodeWidth);
            portContext.portPosition.y = currentYPos + portContext.portMargin.top;
            
            // Update the y coordinate for the next port
            currentYPos += portContext.portMargin.top
                    + portContext.port.getSize().y
                    + portContext.portMargin.bottom
                    + spaceBetweenPorts;
        }
    }

    /**
     * Returns the x coordinate for the given port.
     */
    private static double calculateVerticalPortXCoordinate(final PortContext portContext, final double nodeWidth) {
        PortAdapter<?> port = portContext.port;
        // The x coordinate is set according to the port's offset, if any
        if (port.hasProperty(CoreOptions.PORT_BORDER_OFFSET)) {
            return port.getSide() == PortSide.WEST
                    ? -port.getSize().x - port.getProperty(CoreOptions.PORT_BORDER_OFFSET)
                    : nodeWidth + port.getProperty(CoreOptions.PORT_BORDER_OFFSET);
        } else {
            return port.getSide() == PortSide.WEST
                    ? -port.getSize().x
                    : nodeWidth;
        }
    }
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Utilities
    
    /**
     * If we switch from distributed or justified alignment back to centered alignment, this may require us to modify
     * the required port placement size calculated in a previous phase. This is the case if the old alignment was
     * distributed, in which case the placement size includes two port-port spacings we don't need anymore.
     */
    private static double modifiedPortPlacementSize(final NodeContext nodeContext,
            final PortAlignment oldPortAlignment, final double currentPortPlacementSize) {
        
        if (oldPortAlignment == PortAlignment.DISTRIBUTED) {
            return currentPortPlacementSize - 2 * nodeContext.portPortSpacing;
        } else {
            return currentPortPlacementSize;
        }
    }

}
