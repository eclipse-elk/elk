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
import java.util.EnumSet;

import org.eclipse.elk.alg.common.nodespacing.internal.NodeContext;
import org.eclipse.elk.alg.common.nodespacing.internal.PortContext;
import org.eclipse.elk.core.math.ElkPadding;
import org.eclipse.elk.core.math.ElkRectangle;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.options.PortLabelPlacement;
import org.eclipse.elk.core.options.PortSide;
import org.eclipse.elk.core.options.SizeConstraint;
import org.eclipse.elk.core.options.SizeOptions;
import org.eclipse.elk.core.util.ElkUtil;

/**
 * Various little methods that didn't quite fit into any of the other classes.
 */
public final class NodeLabelAndSizeUtilities {

    /**
     * No instance required.
     */
    private NodeLabelAndSizeUtilities() {
        
    }
    

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Algorithm Phase Implementation Fragments
    
    /**
     * If the client area should have a minimum size, the inside node label container is setup accordingly. This does
     * not yet set any size contribution constraints, but only applies the minimum size itself.
     */
    public static void setupMinimumClientAreaSize(final NodeContext nodeContext) {
        KVector minSize = getMinimumClientAreaSize(nodeContext);
        if (minSize != null) {
            nodeContext.insideNodeLabelContainer.setCenterCellMinimumSize(minSize);
        }
    }
    
    /**
     * Sets up the padding of the main node cell to account for ports that extend inside the node.
     */
    public static void setupNodePaddingForPortsWithOffset(final NodeContext nodeContext) {
        ElkPadding nodeCellPadding = nodeContext.nodeContainer.getPadding();
        
        for (PortContext portContext : nodeContext.portContexts.values()) {
            // If the port extends into the node, ensure the inside port space is enough
            double portBorderOffset = 0;
            if (portContext.port.hasProperty(CoreOptions.PORT_BORDER_OFFSET)) {
                portBorderOffset = portContext.port.getProperty(CoreOptions.PORT_BORDER_OFFSET);
                
                if (portBorderOffset < 0) {
                    // The port does extend into the node, by -portBorderOffset
                    switch (portContext.port.getSide()) {
                    case NORTH:
                        nodeCellPadding.top = Math.max(nodeCellPadding.top, -portBorderOffset);
                        break;
                        
                    case SOUTH:
                        nodeCellPadding.bottom = Math.max(nodeCellPadding.bottom, -portBorderOffset);
                        break;
                        
                    case EAST:
                        nodeCellPadding.right = Math.max(nodeCellPadding.right, -portBorderOffset);
                        break;
                        
                    case WEST:
                        nodeCellPadding.left = Math.max(nodeCellPadding.left, -portBorderOffset);
                        break;
                    }
                }
            }
            if (PortLabelPlacement.isFixed(nodeContext.portLabelsPlacement)) {
                // Ensure that the part of the fixed label, that is inside the node, has enough space.
                double insidePart =
                        ElkUtil.computeInsidePart(portContext.port, portBorderOffset);
                boolean symmetry =
                        !nodeContext.node.getProperty(CoreOptions.NODE_SIZE_OPTIONS).contains(SizeOptions.ASYMMETRICAL);
                boolean insidePartIsBigger = false;
                switch (portContext.port.getSide()) {
                case NORTH:
                    insidePartIsBigger = insidePart > nodeCellPadding.top;
                    nodeCellPadding.top = Math.max(nodeCellPadding.top, insidePart);
                    if (symmetry && insidePartIsBigger) {
                        nodeCellPadding.top = Math.max(nodeCellPadding.top, nodeCellPadding.bottom);
                        // For symmetry, the portBorderOffset is not considered (only label + label padding)
                        nodeCellPadding.bottom = nodeCellPadding.top + portBorderOffset;
                    }
                    break;

                case SOUTH:
                    insidePartIsBigger = insidePart > nodeCellPadding.bottom;
                    nodeCellPadding.bottom = Math.max(nodeCellPadding.bottom, insidePart);
                    if (symmetry && insidePartIsBigger) {
                        nodeCellPadding.bottom = Math.max(nodeCellPadding.bottom, nodeCellPadding.top);
                        // For symmetry, the portBorderOffset is not considered (only label + label padding)
                        nodeCellPadding.top = nodeCellPadding.bottom + portBorderOffset;
                    }
                    break;

                case EAST:
                    insidePartIsBigger = insidePart > nodeCellPadding.right;
                    nodeCellPadding.right = Math.max(nodeCellPadding.right, insidePart);
                    if (symmetry && insidePartIsBigger) {
                        nodeCellPadding.right = Math.max(nodeCellPadding.left, nodeCellPadding.right);
                        // For symmetry, the portBorderOffset is not considered (only label + label padding)
                        nodeCellPadding.left = nodeCellPadding.right + portBorderOffset;

                    }
                    break;

                case WEST:
                    insidePartIsBigger = insidePart > nodeCellPadding.left;
                    nodeCellPadding.left = Math.max(nodeCellPadding.left, insidePart);
                    if (symmetry && insidePartIsBigger) {
                        nodeCellPadding.left = Math.max(nodeCellPadding.left, nodeCellPadding.right);
                        // For symmetry, the portBorderOffset is not considered (only label + label padding)
                        nodeCellPadding.right = nodeCellPadding.left + portBorderOffset;
                    }
                    break;
                }
            }
        }
    }
    
    /**
     * Offsets all southern ports, previously placed relative to vertical coordinate 0, by the node's size to place
     * them along the node's southern border.
     */
    public static void offsetSouthernPortsByNodeSize(final NodeContext nodeContext) {
        double nodeHeight = nodeContext.nodeSize.y;
        
        for (PortContext portContext : nodeContext.portContexts.get(PortSide.SOUTH)) {
            KVector portPosition = portContext.portPosition;
            portPosition.y += nodeHeight;
        }
    }
    
    /**
     * Calculates and stores the node padding, if requested by layout options.
     */
    public static void setNodePadding(final NodeContext nodeContext) {
        if (!nodeContext.sizeOptions.contains(SizeOptions.COMPUTE_PADDING)) {
            return;
        }
        
        ElkRectangle nodeRect = nodeContext.nodeContainer.getCellRectangle();
        ElkRectangle clientArea = nodeContext.insideNodeLabelContainer.getCenterCellRectangle();
        ElkPadding nodePadding = new ElkPadding();
        
        // The following code assumes that the client area rectangle lies fully inside the node rectangle, which should
        // always be the case because of how the client area rectangle is computed
        nodePadding.left = clientArea.x - nodeRect.x;
        nodePadding.top = clientArea.y - nodeRect.y;
        nodePadding.right = (nodeRect.x + nodeRect.width) - (clientArea.x + clientArea.width);
        nodePadding.bottom = (nodeRect.y + nodeRect.height) - (clientArea.y + clientArea.height);
        
        nodeContext.node.setPadding(nodePadding);
    }

    /**
     * Applies the calculated node size to the node and the calculated port positions to the ports.
     */
    public static void applyStuff(final NodeContext nodeContext) {
        nodeContext.applyNodeSize();
        nodeContext.portContexts.values().stream().forEach(pc -> pc.applyPortPosition());
    }
    

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Minimum Size Things
    
    /**
     * Returns the client area's minimum size if size constraints are configured such that the client area has a
     * minimum size. Otherwise, returns {@code null}.
     */
    public static KVector getMinimumClientAreaSize(final NodeContext nodeContext) {
        if (nodeContext.sizeConstraints.contains(SizeConstraint.MINIMUM_SIZE)
                && nodeContext.sizeOptions.contains(SizeOptions.MINIMUM_SIZE_ACCOUNTS_FOR_PADDING)) {
            
            return getMinimumNodeOrClientAreaSize(nodeContext);
            
        } else {
            return null;
        }
    }
    
    /**
     * Returns the node's minimum size if size constraints are configured such that the node as a whole has a minimum
     * size. Otherwise, returns {@code null}.
     */
    public static KVector getMinimumNodeSize(final NodeContext nodeContext) {
        if (nodeContext.sizeConstraints.contains(SizeConstraint.MINIMUM_SIZE)) {
            if (!nodeContext.sizeOptions.contains(SizeOptions.MINIMUM_SIZE_ACCOUNTS_FOR_PADDING)) {
                return getMinimumNodeOrClientAreaSize(nodeContext);
            }
        }
        
        return null;
    }
    
    /**
     * Returns the minimum size configured for the node, without regard for size constraints.
     */
    public static KVector getMinimumNodeOrClientAreaSize(final NodeContext nodeContext) {
        // Retrieve the minimum size
        KVector minSize = new KVector(nodeContext.node.getProperty(CoreOptions.NODE_SIZE_MINIMUM));
        
        // If we are instructed to revert to a default minimum size, we check whether we need to revert to that
        if (nodeContext.sizeOptions.contains(SizeOptions.DEFAULT_MINIMUM_SIZE)) {
            if (minSize.x <= 0) {
                minSize.x = ElkUtil.DEFAULT_MIN_WIDTH;
            }

            if (minSize.y <= 0) {
                minSize.y = ElkUtil.DEFAULT_MIN_HEIGHT;
            }
        }
        
        return minSize;
    }
    

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Utilities
    
    /** Size constraints that only contain port labels should not cause a node to resize. */
    private static final EnumSet<SizeConstraint> EFFECTIVELY_FIXED_SIZE_CONSTRAINTS =
            EnumSet.of(SizeConstraint.PORT_LABELS);
    
    /**
     * Checks if the size constraints, even if not empty, should cause the node not to be resized. This is the case if
     * the constraints are either empty or contain only {@link SizeConstraint#PORT_LABELS}.
     */
    public static boolean areSizeConstraintsFixed(final NodeContext nodeContext) {
        return nodeContext.sizeConstraints.isEmpty()
                || (nodeContext.sizeConstraints.equals(EFFECTIVELY_FIXED_SIZE_CONSTRAINTS));
    }
    
    /**
     * Outside ports usually have their labels placed below or to the right. The first port, however, may have its label
     * placed on the other side. This is true if several conditions apply:
     * <ul>
     *   <li>The port does <b>not</b> want its label to be placed next to it (see
     *       {@link PortLabelPlacement#NEXT_TO_PORT_IF_POSSIBLE}.</li>
     *   <li>{@link PortLabelPlacement#ALWAYS_SAME_SIDE} is <b>not</b> set.</li>
     *   <li>There are only two ports or space-efficient port label placement is turned on (see
     *       {@link PortLabelPlacement#SPACE_EFFICIENT}).</li>
     * </ul>
     */
    public static boolean isFirstOutsidePortLabelPlacedDifferently(final NodeContext nodeContext,
            final PortSide portSide) {
        
        Collection<PortContext> portContexts = nodeContext.portContexts.get(portSide);
        if (portContexts.size() >= 2) {
            PortContext firstPort = portContexts.iterator().next();
            
            final boolean alwaysSameSide =
                    nodeContext.portLabelsPlacement.contains(PortLabelPlacement.ALWAYS_SAME_SIDE);
            final boolean spaceEfficient = nodeContext.portLabelsPlacement.contains(PortLabelPlacement.SPACE_EFFICIENT);

            return !firstPort.labelsNextToPort && !alwaysSameSide
                    && (portContexts.size() == 2 || spaceEfficient);
        } else {
            return false;
        }
    }
    
}
