/*******************************************************************************
 * Copyright (c) 2017, 2020 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.common.nodespacing.internal.algorithm;

import java.util.Collection;
import java.util.Iterator;

import org.eclipse.elk.alg.common.nodespacing.cellsystem.AtomicCell;
import org.eclipse.elk.alg.common.nodespacing.internal.NodeContext;
import org.eclipse.elk.alg.common.nodespacing.internal.PortContext;
import org.eclipse.elk.core.math.ElkRectangle;
import org.eclipse.elk.core.options.PortAlignment;
import org.eclipse.elk.core.options.PortLabelPlacement;
import org.eclipse.elk.core.options.PortSide;
import org.eclipse.elk.core.options.SizeConstraint;
import org.eclipse.elk.core.options.SizeOptions;
import org.eclipse.elk.core.util.ElkUtil;

/**
 * Calculates the space required to setup port labels.
 */
public final class VerticalPortPlacementSizeCalculator {

    /**
     * No instance required.
     */
    private VerticalPortPlacementSizeCalculator() {
        
    }
    
    
    /**
     * Calculates the space required for vertical port placements. If the port placement is not fixed, this will
     * also setup the top and bottom padding of the inside port label placement cells to reflect the amount of space
     * to be left around the port placement.
     */
    public static void calculateVerticalPortPlacementSize(final NodeContext nodeContext) {
        // We need to calculate the space required by the ports even if ports are not part of the size constraints,
        // since we will use that later to place them
        switch (nodeContext.portConstraints) {
        case FIXED_POS:
            // We don't have any freedom at all, so simply calculate where the bottommost port is on each side
            calculateVerticalNodeSizeRequiredByFixedPosPorts(nodeContext, PortSide.EAST);
            calculateVerticalNodeSizeRequiredByFixedPosPorts(nodeContext, PortSide.WEST);
            break;
            
        case FIXED_RATIO:
            // We can require the node to be large enough to avoid spacing violations with fixed ratio ports
            calculateVerticalNodeSizeRequiredByFixedRatioPorts(nodeContext, PortSide.EAST);
            calculateVerticalNodeSizeRequiredByFixedRatioPorts(nodeContext, PortSide.WEST);
            break;
            
        default:
            // If we are free to place things, make the node large enough to place everything properly
            calculateVerticalNodeSizeRequiredByFreePorts(nodeContext, PortSide.EAST);
            calculateVerticalNodeSizeRequiredByFreePorts(nodeContext, PortSide.WEST);
            break;
        }
    }
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Fixed Position
    
    /**
     * Sets up inside port label spaces for fixed position ports. If the node size constraints include ports, we also
     * mess with the cell padding.
     */
    private static void calculateVerticalNodeSizeRequiredByFixedPosPorts(final NodeContext nodeContext,
            final PortSide portSide) {
        
        double bottommostPortBorder = 0.0;
        
        // Check all ports on the correct side
        for (PortContext portContext : nodeContext.portContexts.get(portSide)) {
            bottommostPortBorder = Math.max(
                    bottommostPortBorder,
                    portContext.portPosition.y + portContext.port.getSize().y);
        }
        
        // Set the cell size and remove top padding since the cell size itself already includes all the space we need
        AtomicCell cell = nodeContext.insidePortLabelCells.get(portSide);
        cell.getPadding().top = 0;
        cell.getMinimumContentAreaSize().y = bottommostPortBorder;
    }
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Fixed Ratio
    
    /**
     * Sets up inside port label spaces for fixed ratio ports. If the node size constraints include ports, we also
     * mess with the cell padding.
     */
    private static void calculateVerticalNodeSizeRequiredByFixedRatioPorts(final NodeContext nodeContext,
            final PortSide portSide) {
        
        AtomicCell cell = nodeContext.insidePortLabelCells.get(portSide);
        
        // Fetch the port contexts on the given side and abort if there are none
        Collection<PortContext> portContexts = nodeContext.portContexts.get(portSide);
        if (portContexts.isEmpty()) {
            cell.getPadding().top = 0;
            cell.getPadding().bottom = 0;
            return;
        }

        boolean portLabelsInside = nodeContext.portLabelsPlacement.contains(PortLabelPlacement.INSIDE);
        double minHeight = 0;
        
        // If port labels are to be respected, we need to calculate the port's margins to do so (we ignore the
        // first port in this process if space-efficient mode is active, port labels are to be placed outside
        // and the port labels are not placed next to the port)
        if (nodeContext.sizeConstraints.contains(SizeConstraint.PORT_LABELS)) {
            setupPortMargins(nodeContext, portSide);
        }
        
        // Go over all pairs of consecutive ports
        Iterator<PortContext> portContextIterator = portContexts.iterator();
        PortContext previousPortContext = null;
        double previousPortRatio = 0;
        double previousPortHeight = 0;
        
        while (portContextIterator.hasNext()) {
            // Get the next port and find out things about it
            PortContext currentPortContext = portContextIterator.next();
            double currentPortRatio = currentPortContext.port.getProperty(
                    PortPlacementCalculator.PORT_RATIO_OR_POSITION);
            double currentPortHeight = currentPortContext.port.getSize().y;
            
            if (previousPortContext == null) {
                // This is the first port, so find out how high the node needs to be to respect the top surrounding
                // port margins, if any
                if (nodeContext.surroundingPortMargins != null && nodeContext.surroundingPortMargins.top > 0) {
                    minHeight = Math.max(
                            minHeight,
                            HorizontalPortPlacementSizeCalculator.minSizeRequiredToRespectSpacing(
                                    nodeContext.surroundingPortMargins.top + currentPortContext.portMargin.top,
                                    0,
                                    currentPortRatio));
                }
            } else {
                double requiredSpace = previousPortHeight + previousPortContext.portMargin.bottom
                        + nodeContext.portPortSpacing + currentPortContext.portMargin.top;
                minHeight = Math.max(
                        minHeight,
                        HorizontalPortPlacementSizeCalculator.minSizeRequiredToRespectSpacing(
                                requiredSpace,
                                previousPortRatio,
                                currentPortRatio));
            }
            
            // Our current port is going to be the previous port during the next iteration
            previousPortContext = currentPortContext;
            previousPortRatio = currentPortRatio;
            previousPortHeight = currentPortHeight;
        }
        
        // if there are bottom surrounding port margins, apply those as well
        if (nodeContext.surroundingPortMargins != null && nodeContext.surroundingPortMargins.bottom > 0) {
            // We're using the port's bare width here because we don't care about label sizes on the bottommost port
            double requiredSpace = previousPortHeight + nodeContext.surroundingPortMargins.bottom;
            
            // We're only interested in the port's bottom margin if it's label is placed inside
            if (portLabelsInside) {
                requiredSpace += previousPortContext.portMargin.bottom;
            }
            
            minHeight = Math.max(
                    minHeight,
                    HorizontalPortPlacementSizeCalculator.minSizeRequiredToRespectSpacing(
                            requiredSpace,
                            previousPortRatio,
                            1));
        }
        
        // Set the cell size and remove top padding since the cell size itself already includes all the space we need
        cell.getPadding().top = 0;
        cell.getMinimumContentAreaSize().y = minHeight;
    }
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Free
    
    /**
     * Sets up inside port label spaces for free ports.
     */
    private static void calculateVerticalNodeSizeRequiredByFreePorts(final NodeContext nodeContext,
            final PortSide portSide) {

        AtomicCell cell = nodeContext.insidePortLabelCells.get(portSide);
        
        // Handle the common case first: if there are no ports, set everything to zero
        if (nodeContext.portContexts.get(portSide).isEmpty()) {
            cell.getPadding().top = 0;
            cell.getPadding().bottom = 0;
            return;
        }
        
        // Set the padding to match the surrounding port space
        cell.getPadding().top = nodeContext.surroundingPortMargins.top;
        cell.getPadding().bottom = nodeContext.surroundingPortMargins.bottom;
        
        // If we are to take labels into account, we need to setup the port margins such that they include the space
        // required for their labels
        if (nodeContext.sizeConstraints.contains(SizeConstraint.PORT_LABELS)) {
            setupPortMargins(nodeContext, portSide);
        }
        
        double height = portHeightPlusPortPortSpacing(nodeContext, portSide);

        // For distributed port alignment, we need to surround the ports by a port-port spacing on each side
        if (nodeContext.getPortAlignment(portSide) == PortAlignment.DISTRIBUTED) {
            height += 2 * nodeContext.portPortSpacing;
        }
        
        // Set the cell size 
        cell.getMinimumContentAreaSize().y = height;
    }

    /**
     * Sets up the port margins of all ports such that they include any space required for labels, subject to any other
     * options that might affect label placement.
     */
    private static void setupPortMargins(final NodeContext nodeContext, final PortSide portSide) {
        Collection<PortContext> portContexts = nodeContext.portContexts.get(portSide);

        boolean portLabelsOutside = nodeContext.portLabelsPlacement.contains(PortLabelPlacement.OUTSIDE);
        boolean alwaysSameSide = nodeContext.portLabelsPlacement.contains(PortLabelPlacement.ALWAYS_SAME_SIDE);
        boolean alwaysSameSideAbove =
                nodeContext.portLabelsPlacement.contains(PortLabelPlacement.ALWAYS_OTHER_SAME_SIDE);
        boolean spaceEfficient = nodeContext.portLabelsPlacement.contains(PortLabelPlacement.SPACE_EFFICIENT);
        boolean uniformPortSpacing = nodeContext.sizeOptions.contains(SizeOptions.UNIFORM_PORT_SPACING);

        boolean spaceEfficientPortLabels = !alwaysSameSide && !alwaysSameSideAbove
                && (spaceEfficient || portContexts.size() == 2);

        // Set the vertical port margins of all ports according to how their labels will be placed. We'll be
        // modifying the margins soon enough.
        computeVerticalPortMargins(nodeContext, portSide, portLabelsOutside);
        
        // The topmost and bottommost ports are possibly required
        PortContext topmostPortContext = null;
        PortContext bottommostPortContext = null;
        
        // If port labels are placed outside, there's stuff we can do
        if (portLabelsOutside) {
            // Find topmost and bottommost ports
            Iterator<PortContext> portContextIterator = portContexts.iterator();
            
            topmostPortContext = portContextIterator.next();
            bottommostPortContext = topmostPortContext;
            
            while (portContextIterator.hasNext()) {
                bottommostPortContext = portContextIterator.next();
            }
            
            // The topmost and bottommost ports don't need their top and bottom margin, respectively
            topmostPortContext.portMargin.top = 0;
            bottommostPortContext.portMargin.bottom = 0;
            
            // If we place port labels space-efficiently and the topmost port doesn't have its label placed right
            // next to it, it doesn't need its bottom margin either since its label will be placed above
            if (spaceEfficientPortLabels && !topmostPortContext.labelsNextToPort) {
                topmostPortContext.portMargin.bottom = 0;
            }
        }
        
        // If ports are placed uniformly, we reflect that here by equalizing all port margins
        if (uniformPortSpacing) {
            unifyPortMargins(portContexts);
            
            // Uniforming may have reset the topmost port's top and bottommost port's bottom margins 
            if (portLabelsOutside) {
                topmostPortContext.portMargin.top = 0;
                bottommostPortContext.portMargin.bottom = 0;
            }
        }
    }
    
    /**
     * Sets the vertical margins of all ports such that they include the space necessary to place their labels.
     */
    private static void computeVerticalPortMargins(final NodeContext nodeContext, final PortSide portSide,
            final boolean portLabelsOutside) {
        
        for (PortContext portContext : nodeContext.portContexts.get(portSide)) {
            double labelHeight = portContext.portLabelCell != null
                    ? portContext.portLabelCell.getMinimumHeight()
                    : 0;
            
            if (labelHeight > 0) {
                if (portContext.labelsNextToPort) {
                    // The label is placed next to the port
                    double portHeight = portContext.port.getSize().y;
                    if (labelHeight > portHeight) {
                        if (nodeContext.portLabelsTreatAsGroup || portContext.portLabelCell.getLabels().size() == 1) {
                            // We are to center all of the labels
                            double overhang = (labelHeight - portHeight) / 2;
                            portContext.portMargin.top = overhang;
                            portContext.portMargin.bottom = overhang;
                            
                        } else {
                            // Simulate centering the first port label
                            double firstLabelHeight = portContext.portLabelCell.getLabels().get(0).getSize().y;
                            double firstLabelOverhang = (firstLabelHeight - portHeight) / 2;
                            
                            portContext.portMargin.top = Math.max(0, firstLabelOverhang);
                            portContext.portMargin.bottom = labelHeight - firstLabelOverhang - portHeight;
                        }
                    }
                    
                } else {
                    // The label is either placed outside (below the port) or possibly inside, but for a compound node,
                    // which means that it is placed below the port as well to keep it from overlapping with inside
                    // edges
                    portContext.portMargin.bottom = nodeContext.portLabelSpacingVertical + labelHeight;
                }
            } else if (PortLabelPlacement.isFixed(nodeContext.portLabelsPlacement)) {
                // The fixed port label is not considered with portContext.portLabelCell. Nevertheless, a port margin
                // must be added if necessary.
                ElkRectangle labelsBounds = ElkUtil.getLabelsBounds(portContext.port);
                if (labelsBounds.y < 0) {
                 // Add the part of the label that is above the port to the top margin
                    portContext.portMargin.top = -labelsBounds.y;
                }
                if (labelsBounds.y + labelsBounds.height > portContext.port.getSize().y) {
                 // Add the part of the label that is below the port to the bottom margin
                    portContext.portMargin.bottom = labelsBounds.y + labelsBounds.height - portContext.port.getSize().y;
                }
            }
        }
    }


    /**
     * Sets all port margins to be equal to the maximum margins.
     */
    private static void unifyPortMargins(final Collection<PortContext> portContexts) {
        double maxTop = 0;
        double maxBottom = 0;
        
        // Find maximum
        for (PortContext portContext : portContexts) {
            maxTop = Math.max(maxTop, portContext.portMargin.top);
            maxBottom = Math.max(maxBottom, portContext.portMargin.bottom);
        }
        
        // Apply maximum
        for (PortContext portContext : portContexts) {
            portContext.portMargin.top = maxTop;
            portContext.portMargin.bottom = maxBottom;
        }
    }

    /**
     * Takes all ports on the given side, sums up their height, and inserts port-port spacings between them. The margins
     * of ports can be included as well to include space used for labels.
     */
    private static double portHeightPlusPortPortSpacing(final NodeContext nodeContext, final PortSide portSide) {
        double result = 0;
        
        Iterator<PortContext> portContextIterator = nodeContext.portContexts.get(portSide).iterator();
        while (portContextIterator.hasNext()) {
            PortContext portContext = portContextIterator.next();
            
            result += portContext.portMargin.top + portContext.port.getSize().y + portContext.portMargin.bottom;
            
            if (portContextIterator.hasNext()) {
                result += nodeContext.portPortSpacing;
            }
        }
        
        return result;
    }
    
}
