/*******************************************************************************
 * Copyright (c) 2017, 2020 Kiel University and others.
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

import com.google.common.math.DoubleMath;

/**
 * Calculates the space required to setup port labels.
 */
public final class HorizontalPortPlacementSizeCalculator {

    /**
     * No instance required.
     */
    private HorizontalPortPlacementSizeCalculator() {
        
    }
    
    
    /**
     * Calculates the space required for horizontal port placements. If the port placement is not fixed, this will
     * also setup the left and right padding of the inside port label placement cells to reflect the amount of space
     * to be left around the port placement.
     */
    public static void calculateHorizontalPortPlacementSize(final NodeContext nodeContext) {
        // We need to calculate the space required by the ports even if ports are not part of the size constraints,
        // since we will use that later to place them
        switch (nodeContext.portConstraints) {
        case FIXED_POS:
            // We don't have any freedom at all, so simply calculate where the rightmost port is on each side
            calculateHorizontalNodeSizeRequiredByFixedPosPorts(nodeContext, PortSide.NORTH);
            calculateHorizontalNodeSizeRequiredByFixedPosPorts(nodeContext, PortSide.SOUTH);
            break;
            
        case FIXED_RATIO:
            // We can require the node to be large enough to avoid spacing violations with fixed ratio ports
            calculateHorizontalNodeSizeRequiredByFixedRatioPorts(nodeContext, PortSide.NORTH);
            calculateHorizontalNodeSizeRequiredByFixedRatioPorts(nodeContext, PortSide.SOUTH);
            break;
            
        default:
            // If we are free to place things, make the node large enough to place everything properly
            calculateHorizontalNodeSizeRequiredByFreePorts(nodeContext, PortSide.NORTH);
            calculateHorizontalNodeSizeRequiredByFreePorts(nodeContext, PortSide.SOUTH);
            break;
        }
    }
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Fixed Position
    
    /**
     * Sets up inside port label spaces for fixed position ports. If the node size constraints include ports, we also
     * mess with the cell padding.
     */
    private static void calculateHorizontalNodeSizeRequiredByFixedPosPorts(final NodeContext nodeContext,
            final PortSide portSide) {
        
        double rightmostPortBorder = 0.0;
        
        // Check all ports on the correct side
        for (PortContext portContext : nodeContext.portContexts.get(portSide)) {
            rightmostPortBorder = Math.max(
                    rightmostPortBorder,
                    portContext.portPosition.x + portContext.port.getSize().x);
        }
        
        // Set the cell size and remove left padding since the cell size itself already includes all the space we need
        AtomicCell cell = nodeContext.insidePortLabelCells.get(portSide);
        cell.getPadding().left = 0;
        cell.getMinimumContentAreaSize().x = rightmostPortBorder;
    }
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Fixed Ratio
    
    /**
     * Sets up inside port label spaces for fixed ratio ports. If the node size constraints include ports, we also
     * mess with the cell padding.
     */
    private static void calculateHorizontalNodeSizeRequiredByFixedRatioPorts(final NodeContext nodeContext,
            final PortSide portSide) {
        
        AtomicCell cell = nodeContext.insidePortLabelCells.get(portSide);
        
        // Fetch the port contexts on the given side and abort if there are none
        Collection<PortContext> portContexts = nodeContext.portContexts.get(portSide);
        if (portContexts.isEmpty()) {
            cell.getPadding().left = 0;
            cell.getPadding().right = 0;
            return;
        }

        boolean portLabelsInside = nodeContext.portLabelsPlacement.contains(PortLabelPlacement.INSIDE);
        double minWidth = 0;
        
        // Go over all pairs of consecutive ports
        Iterator<PortContext> portContextIterator = portContexts.iterator();
        PortContext previousPortContext = null;
        double previousPortRatio = 0;
        double previousPortWidth = 0;
        
        while (portContextIterator.hasNext()) {
            // Get the next port and find out things about it
            PortContext currentPortContext = portContextIterator.next();
            double currentPortRatio = currentPortContext.port.getProperty(
                    PortPlacementCalculator.PORT_RATIO_OR_POSITION);
            double currentPortWidth = currentPortContext.port.getSize().x;
            
            // If port labels are to be respected, we need to calculate the port's margins to do so (we ignore the
            // first port in this process if space-efficient mode is active and port labels are to be placed outside)
            if (nodeContext.sizeConstraints.contains(SizeConstraint.PORT_LABELS)) {
                setupPortMargins(nodeContext, portSide);
            }
            
            if (previousPortContext == null) {
                // This is the first port, so find out how wide the node needs to be to respect the left surrounding
                // port margins, if any
                if (nodeContext.surroundingPortMargins != null && nodeContext.surroundingPortMargins.left > 0) {
                    minWidth = Math.max(
                            minWidth,
                            minSizeRequiredToRespectSpacing(
                                    nodeContext.surroundingPortMargins.left + currentPortContext.portMargin.left,
                                    0,
                                    currentPortRatio));
                }
            } else {
                double requiredSpace = previousPortWidth + previousPortContext.portMargin.right
                        + nodeContext.portPortSpacing + currentPortContext.portMargin.left;
                minWidth = Math.max(
                        minWidth,
                        minSizeRequiredToRespectSpacing(
                                requiredSpace,
                                previousPortRatio,
                                currentPortRatio));
            }
            
            // Our current port is going to be the previous port during the next iteration
            previousPortContext = currentPortContext;
            previousPortRatio = currentPortRatio;
            previousPortWidth = currentPortWidth;
        }
        
        // if there are right surrounding port margins, apply those as well
        if (nodeContext.surroundingPortMargins != null && nodeContext.surroundingPortMargins.right > 0) {
            double requiredSpace = previousPortWidth + nodeContext.surroundingPortMargins.right;
            
            // We're only interested in the port's right margin if it's label is placed inside
            if (portLabelsInside) {
                requiredSpace += previousPortContext.portMargin.right;
            }
            
            minWidth = Math.max(
                    minWidth,
                    minSizeRequiredToRespectSpacing(
                            requiredSpace,
                            previousPortRatio,
                            1));
        }
        
        // Set the cell size and remove left padding since the cell size itself already includes all the space we need
        cell.getPadding().left = 0;
        cell.getMinimumContentAreaSize().x = minWidth;
    }
    
    /** Fuzzyness allowed to still consider two double values to be equal. */
    private static final double EQUALITY_TOLERANCE = 0.01;
    
    /**
     * Returns the minimum width that will satisfy the given spacing between the two ratios multiplied by the result
     * width.
     * 
     * @param spacing
     *            the spacing to respect.
     * @param firstRatio
     *            position ratio of the first element.
     * @param secondRatio
     *            position ratio of the second element.
     * @return minimum width required to satisfy the spacing.
     */
    static double minSizeRequiredToRespectSpacing(final double spacing, final double firstRatio,
            final double secondRatio) {
        
        assert secondRatio >= firstRatio;
        
        // Some failsafing
        if (DoubleMath.fuzzyEquals(firstRatio, secondRatio, EQUALITY_TOLERANCE)) {
            return 0;
        } else {
            return spacing / (secondRatio - firstRatio);
        }
    }
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Free
    
    /**
     * Sets up inside port label spaces for free ports.
     */
    private static void calculateHorizontalNodeSizeRequiredByFreePorts(final NodeContext nodeContext,
            final PortSide portSide) {

        AtomicCell cell = nodeContext.insidePortLabelCells.get(portSide);
        
        // Handle the common case first: if there are no ports, set everything to zero
        if (nodeContext.portContexts.get(portSide).isEmpty()) {
            cell.getPadding().left = 0;
            cell.getPadding().right = 0;
            return;
        }
        
        // Set the padding to match the surrounding port space
        cell.getPadding().left = nodeContext.surroundingPortMargins.left;
        cell.getPadding().right = nodeContext.surroundingPortMargins.right;
        
        // If we are to take labels into account, we need to setup the port margins such that they include the space
        // required for their labels
        if (nodeContext.sizeConstraints.contains(SizeConstraint.PORT_LABELS)) {
            setupPortMargins(nodeContext, portSide);
        }
        
        double width = portWidthPlusPortPortSpacing(nodeContext, portSide);

        // For distributed port alignment, we need to surround the ports by a port-port spacing on each side
        if (nodeContext.getPortAlignment(portSide) == PortAlignment.DISTRIBUTED) {
            width += 2 * nodeContext.portPortSpacing;
        }
        
        // Set the cell size 
        cell.getMinimumContentAreaSize().x = width;
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

        // Set the horizontal port margins of all ports according to how their labels will be placed. We'll be
        // modifying the margins soon enough.
        computeHorizontalPortMargins(nodeContext, portSide, portLabelsOutside);
        
        // The leftmost and rightmost ports are possibly required
        PortContext leftmostPortContext = null;
        PortContext rightmostPortContext = null;
        
        // If port labels are placed outside, there's stuff we can do
        if (portLabelsOutside) {
            // Find leftmost and rightmost ports
            Iterator<PortContext> portContextIterator = portContexts.iterator();
            
            leftmostPortContext = portContextIterator.next();
            rightmostPortContext = leftmostPortContext;
            
            while (portContextIterator.hasNext()) {
                rightmostPortContext = portContextIterator.next();
            }
            
            // The leftmost and rightmost ports don't need their left and right margin, respectively
            leftmostPortContext.portMargin.left = 0;
            rightmostPortContext.portMargin.right = 0;
            
            // If we place port labels space-efficiently and the leftmost port doesn't have its label placed right
            // next to it, it doesn't need its right margin either since its label will be placed to its left
            if (spaceEfficientPortLabels && !leftmostPortContext.labelsNextToPort) {
                leftmostPortContext.portMargin.right = 0;
            }
        }
        
        // If ports are placed uniformly, we reflect that here by equalizing all port margins
        if (uniformPortSpacing) {
            unifyPortMargins(portContexts);
            
            // Uniforming may have reset the leftmost port's left and rightmost port's right margins 
            if (portLabelsOutside) {
                leftmostPortContext.portMargin.left = 0;
                rightmostPortContext.portMargin.right = 0;
            }
        }
    }
    
    /**
     * Sets the horizontal margins of all ports such that they include the space necessary to place their labels.
     */
    private static void computeHorizontalPortMargins(final NodeContext nodeContext, final PortSide portSide,
            final boolean portLabelsOutside) {
        
        for (PortContext portContext : nodeContext.portContexts.get(portSide)) {
            double labelWidth = portContext.portLabelCell != null
                    ? portContext.portLabelCell.getMinimumWidth()
                    : 0;
            
            if (labelWidth > 0) {
                if (portContext.labelsNextToPort) {
                    // The label is placed next to the port
                    double portWidth = portContext.port.getSize().x;
                    if (labelWidth > portWidth) {
                        double overhang = (labelWidth - portWidth) / 2;
                        portContext.portMargin.left = overhang;
                        portContext.portMargin.right = overhang;
                    }
                    
                } else {
                    // The label is either placed outside (right to the port) or possibly inside, but for a compound
                    // node, which means that it is placed right of the port as well to keep it from overlapping with
                    // inside edges
                    portContext.portMargin.right = nodeContext.portLabelSpacingHorizontal + labelWidth;
                }
            } else if (PortLabelPlacement.isFixed(nodeContext.portLabelsPlacement)) {
                // The fixed port label is not considered with portContext.portLabelCell. Nevertheless, a port margin
                // must be added if necessary.
                ElkRectangle labelsBounds = ElkUtil.getLabelsBounds(portContext.port);
                if (labelsBounds.x < 0) {
                    // Add the part of the label that is on the left of the port to the left margin
                    portContext.portMargin.left = -labelsBounds.x;
                }
                if (labelsBounds.x + labelsBounds.width > portContext.port.getSize().x) {
                    // Add the part of the label that is on the right of the port to the right margin
                    portContext.portMargin.right = labelsBounds.x + labelsBounds.width - portContext.port.getSize().x;
                }
            }
        }
    }


    /**
     * Sets all port margins to be equal to the maximum margins.
     */
    private static void unifyPortMargins(final Collection<PortContext> portContexts) {
        double maxLeft = 0;
        double maxRight = 0;
        
        // Find maximum
        for (PortContext portContext : portContexts) {
            maxLeft = Math.max(maxLeft, portContext.portMargin.left);
            maxRight = Math.max(maxRight, portContext.portMargin.right);
        }
        
        // Apply maximum
        for (PortContext portContext : portContexts) {
            portContext.portMargin.left = maxLeft;
            portContext.portMargin.right = maxRight;
        }
    }

    /**
     * Takes all ports on the given side, sums up their width, and inserts port-port spacings between them. The margins
     * of ports can be included as well to include space used for labels.
     */
    private static double portWidthPlusPortPortSpacing(final NodeContext nodeContext, final PortSide portSide) {
        double result = 0;
        
        Iterator<PortContext> portContextIterator = nodeContext.portContexts.get(portSide).iterator();
        while (portContextIterator.hasNext()) {
            PortContext portContext = portContextIterator.next();
            
            // Add the current port's width to our result
            result += portContext.portMargin.left + portContext.port.getSize().x + portContext.portMargin.right;
            
            // If there is another port after this one, include the required spacing
            if (portContextIterator.hasNext()) {
                result += nodeContext.portPortSpacing;
            }
        }
        
        return result;
    }
    
}
