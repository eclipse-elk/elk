/*******************************************************************************
 * Copyright (c) 2019 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.alg.layered.intermediate.loops.routing;

import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LGraphUtil;
import org.eclipse.elk.alg.layered.graph.LLabel;
import org.eclipse.elk.alg.layered.graph.LMargin;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.intermediate.loops.SelfHyperLoop;
import org.eclipse.elk.alg.layered.intermediate.loops.SelfHyperLoopLabels;
import org.eclipse.elk.alg.layered.intermediate.loops.SelfLoopEdge;
import org.eclipse.elk.alg.layered.intermediate.loops.SelfLoopHolder;
import org.eclipse.elk.alg.layered.intermediate.loops.SelfLoopPort;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.math.KVectorChain;
import org.eclipse.elk.core.options.PortSide;

/**
 * Routes self loops orthogonally. The routing functionality is exposed to subclasses. If a subclass wants to base its
 * routing on the orthogonal routing, it can simply override {@link #modifyBendPoints(SelfLoopEdge, KVectorChain)} to
 * turn an orthogonal routing into some kind of a strange, exotic, exciting new routing style!
 */
public class OrthogonalSelfLoopRouter extends AbstractSelfLoopRouter {
    
    /**
     * The direction in which an edge goes around a node in its quest to reach its target.
     */
    protected enum EdgeRoutingDirection {
        /** Clockwise direction. */
        CLOCKWISE,
        /** Counter-clockwise direction. */
        COUNTER_CLOCKWISE;
    }
    
    @Override
    public void routeSelfLoops(final SelfLoopHolder slHolder) {
        LNode lNode = slHolder.getLNode();
        
        KVector nodeSize = lNode.getSize();
        LMargin nodeMargins = lNode.getMargin();
        
        double edgeEdgeDistance = LGraphUtil.getIndividualOrInherited(lNode, LayeredOptions.SPACING_EDGE_EDGE);
        double edgeLabelDistance = LGraphUtil.getIndividualOrInherited(lNode, LayeredOptions.SPACING_EDGE_LABEL);
        double nodeSLDistance = LGraphUtil.getIndividualOrInherited(lNode, LayeredOptions.SPACING_NODE_SELF_LOOP);
        
        LMargin newNodeMargins = new LMargin();
        newNodeMargins.set(nodeMargins);
        
        // Compute how far away from the node each routing slot on each side is (this takes labels into account)
        double[][] routingSlotPositions = computeRoutingSlotPositions(
                slHolder, edgeEdgeDistance, edgeLabelDistance, nodeSLDistance);
        
        for (SelfHyperLoop slLoop : slHolder.getSLHyperLoops()) {
            for (SelfLoopEdge slEdge : slLoop.getSLEdges()) {
                LEdge lEdge = slEdge.getLEdge();
                
                EdgeRoutingDirection routingDirection = computeEdgeRoutingDirection(slEdge);
                
                // Compute orthogonal bend points and give subclasses a chance to modify them to suit their particular
                // routing style. The default implementation in this class won't change the bend points.
                KVectorChain bendPoints = computeOrthogonalBendPoints(slEdge, routingDirection, routingSlotPositions);
                bendPoints = modifyBendPoints(slEdge, routingDirection, bendPoints);
                
                lEdge.getBendPoints().clear();
                lEdge.getBendPoints().addAll(bendPoints);
                
                bendPoints.stream().forEach(bp -> updateNewNodeMargins(nodeSize, newNodeMargins, bp));
            }
            
            // Place the self loop's labels (the edges were routed such that there is enough space available)
            SelfHyperLoopLabels slLabels = slLoop.getSLLabels();
            if (slLabels != null) {
                placeLabels(slLoop, slLabels, routingSlotPositions, edgeLabelDistance);
                updateNewNodeMargins(nodeSize, newNodeMargins, slLabels);
            }
            
        }
        
        // Update the node's margins to include the space required for self loops
        nodeMargins.set(newNodeMargins);
    }

    /**
     * Computes how the edge reaches its target.
     */
    private EdgeRoutingDirection computeEdgeRoutingDirection(final SelfLoopEdge slEdge) {
        LPort sourceLPort = slEdge.getSLSource().getLPort();
        PortSide sourcePortSide = sourceLPort.getSide();
        
        LPort targetLPort = slEdge.getSLTarget().getLPort();
        PortSide targetPortSide = targetLPort.getSide();
        
        if (sourcePortSide == targetPortSide) {
            // If the port sides are equal, it doesn't really matter what we return, but we'll make a best effort based
            // on port IDs
            return sourceLPort.id < targetLPort.id
                    ? EdgeRoutingDirection.CLOCKWISE
                    : EdgeRoutingDirection.COUNTER_CLOCKWISE;
        
        } else if (sourcePortSide.right() == targetPortSide) {
            return EdgeRoutingDirection.CLOCKWISE;
            
        } else if (sourcePortSide.left() == targetPortSide) {
            return EdgeRoutingDirection.COUNTER_CLOCKWISE;
            
        } else {
            assert sourcePortSide.opposed() == targetPortSide;
            
            // What we do here totally depends on the port sides occupied by the self hyper loop. We prefer clockwise
            // routing.
            SelfHyperLoop slLoop = slEdge.getSLHyperLoop();
            if (slLoop.getOccupiedPortSides().contains(sourcePortSide.right())) {
                return EdgeRoutingDirection.CLOCKWISE;
            } else {
                assert slLoop.getOccupiedPortSides().contains(sourcePortSide.left());
                return EdgeRoutingDirection.COUNTER_CLOCKWISE;
            }
        }
    }

    /**
     * Places any labels of the given self loop.
     */
    private void placeLabels(final SelfHyperLoop slLoop, final SelfHyperLoopLabels slLabels,
            final double[][] routingSlotPositions, double edgeLabelDistance) {
        
        // Find the baseline of the routing slot (we need to offset this by the spacing to be left between label and
        // edge)
        PortSide labelSide = slLabels.getSide();
        double labelPosition = routingSlotPositions[labelSide.ordinal()][slLoop.getRoutingSlot(labelSide)];
        boolean inline = false;
        for (LLabel label : slLabels.getLLabels()) {
            if (label.getProperty(LayeredOptions.EDGE_LABELS_INLINE)) {
                inline = true;
                break;
            }
        }
        // Do not use edgeLabelDistance if the labels are inline
        if (inline) {
            edgeLabelDistance = 0;
        }
        
        switch (labelSide) {
        case NORTH:
            labelPosition -= edgeLabelDistance + slLabels.getSize().y;
            slLabels.getPosition().y = labelPosition;
            break;
            
        case SOUTH:
            labelPosition += edgeLabelDistance;
            slLabels.getPosition().y = labelPosition;
            break;
            
        case WEST:
            labelPosition -= edgeLabelDistance + slLabels.getSize().x;
            slLabels.getPosition().x = labelPosition;
            break;
            
        case EAST:
            labelPosition += edgeLabelDistance;
            slLabels.getPosition().x = labelPosition;
            break;
        
        default:
            assert false;
        }
    }

    /**
     * Extends the node margins to include the given bend point.
     */
    private void updateNewNodeMargins(final KVector nodeSize, final LMargin newNodeMargins, final KVector bendPoint) {
        newNodeMargins.left = Math.max(newNodeMargins.left, -bendPoint.x);
        newNodeMargins.right = Math.max(newNodeMargins.right, bendPoint.x - nodeSize.x);
        
        newNodeMargins.top = Math.max(newNodeMargins.top, -bendPoint.y);
        newNodeMargins.bottom = Math.max(newNodeMargins.bottom, bendPoint.y - nodeSize.y);
    }
    
    /**
     * Ensures that self loop labels are fully covered by the node margins.
     */
    private void updateNewNodeMargins(final KVector nodeSize, final LMargin newNodeMargins,
            final SelfHyperLoopLabels slLabels) {
        
        KVector pos = new KVector(slLabels.getPosition());
        updateNewNodeMargins(nodeSize, newNodeMargins, pos);
        
        pos.add(slLabels.getSize());
        updateNewNodeMargins(nodeSize, newNodeMargins, pos);
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Routing Slot Positions
    
    /**
     * Computes the position of each routing slot on each side, taking care to leave enough space for labels between
     * adjacent routing slots on the north and south sides.
     */
    private double[][] computeRoutingSlotPositions(final SelfLoopHolder slHolder, final double edgeEdgeDistance,
            final double edgeLabelDistance, final double nodeSLDistance) {
        
        // Initialize array
        double[][] positions = new double[PortSide.values().length][];
        for (PortSide portSide : PortSide.values()) {
            positions[portSide.ordinal()] = new double[slHolder.getRoutingSlotCount()[portSide.ordinal()]];
        }
        
        // To know how much space we need to leave between adjacent routing slots, we have to find the size of labels
        // first (for north and south sides, that is)
        initializeWithMaxLabelHeight(positions, slHolder, PortSide.NORTH);
        initializeWithMaxLabelHeight(positions, slHolder, PortSide.SOUTH);
        
        // Compute the positions for each side
        computePositions(positions, slHolder, PortSide.NORTH, edgeEdgeDistance, edgeLabelDistance, nodeSLDistance);
        computePositions(positions, slHolder, PortSide.EAST, edgeEdgeDistance, edgeLabelDistance, nodeSLDistance);
        computePositions(positions, slHolder, PortSide.SOUTH, edgeEdgeDistance, edgeLabelDistance, nodeSLDistance);
        computePositions(positions, slHolder, PortSide.WEST, edgeEdgeDistance, edgeLabelDistance, nodeSLDistance);
        
        return positions;
    }

    /**
     * Initializes the entry for each routing slot on the north or south sides with the maximum height of labels it
     * houses.
     */
    private void initializeWithMaxLabelHeight(final double[][] positions, final SelfLoopHolder slHolder,
            final PortSide portSide) {
        
        assert portSide == PortSide.NORTH || portSide == PortSide.SOUTH;
        
        double[] sidePositions = positions[portSide.ordinal()];
        
        for (SelfHyperLoop slLoop : slHolder.getSLHyperLoops()) {
            SelfHyperLoopLabels slLabels = slLoop.getSLLabels();
            if (slLabels != null && slLabels.getSide() == portSide) {
                int routingSlot = slLoop.getRoutingSlot(portSide);
                sidePositions[routingSlot] = Math.max(sidePositions[routingSlot], slLabels.getSize().y);
            }
        }
    }
    
    /**
     * Computes the positions for each routing slot, taking labels into account.
     */
    private void computePositions(final double[][] positions, final SelfLoopHolder slHolder, final PortSide portSide,
            final double edgeEdgeDistance, final double edgeLabelDistance, final double nodeSelfLoopDistance) {
        
        double currPos = computeBaselinePosition(slHolder, portSide, nodeSelfLoopDistance);
        
        // For northern and western coordinates, we have to subtract from the current position
        double factor = portSide == PortSide.NORTH || portSide == PortSide.WEST ? -1 : 1;
        
        double[] sidePositions = positions[portSide.ordinal()];
        for (int slot = 0; slot < sidePositions.length; slot++) {
            // The slot entry currently contains the height or width of the largest label in that slot
            double largestLabelSize = sidePositions[slot];
            if (largestLabelSize > 0) {
                // Account for label spacing
                largestLabelSize += edgeLabelDistance;
            }
            
            // Place the slot at the current position and advance the position
            sidePositions[slot] = currPos;
            currPos += factor * (largestLabelSize + edgeEdgeDistance);
        }
    }

    /**
     * Based on the margin area, this computes the offset from the node origin to add to escape the area occupied by
     * ports.
     */
    private double computeBaselinePosition(final SelfLoopHolder slHolder, final PortSide portSide,
            final double nodeSelfLoopDistance) {
        
        LNode lNode = slHolder.getLNode();
        LMargin lMargins = lNode.getMargin();
        
        switch (portSide) {
        case NORTH:
            return -lMargins.top - nodeSelfLoopDistance;
        case EAST:
            return lNode.getSize().x + lMargins.right + nodeSelfLoopDistance;
        case SOUTH:
            return lNode.getSize().y + lMargins.bottom + nodeSelfLoopDistance;
        case WEST:
            return -lMargins.left - nodeSelfLoopDistance;
        default:
            assert false;
            return -1;
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // High-Level Bend Point Computation

    /**
     * Computes the bend points necessary to route the given self loop edge orthogonally.
     */
    protected KVectorChain computeOrthogonalBendPoints(final SelfLoopEdge slEdge,
            final EdgeRoutingDirection routingDirection, final double[][] routingSlotPositions) {
        
        KVectorChain bendPoints = new KVectorChain();
        // Make sure that you consider the label size for centering inline.
        addOuterBendPoint(slEdge, slEdge.getSLSource(), routingSlotPositions, bendPoints);
        addCornerBendPoints(slEdge, routingDirection, routingSlotPositions, bendPoints);
        addOuterBendPoint(slEdge, slEdge.getSLTarget(), routingSlotPositions, bendPoints);
        
        return bendPoints;
    }
    
    /**
     * Allows subclasses to turn the given list of bend points into a new list of bend points.
     */
    protected KVectorChain modifyBendPoints(final SelfLoopEdge slEdge, final EdgeRoutingDirection routingDirection,
            final KVectorChain bendPoints) {
        
        return bendPoints;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Actual Bend Point Computation

    private void addOuterBendPoint(final SelfLoopEdge slEdge, final SelfLoopPort slPort,
            final double[][] routingSlotPositions, final KVectorChain bendPoints) {
        
        SelfHyperLoop slLoop = slEdge.getSLHyperLoop();
        LPort lPort = slPort.getLPort();
        PortSide portSide = lPort.getSide();
        
        // We'll start by computing the coordinate of the level we're on
        KVector result = getBaseVector(portSide, slLoop.getRoutingSlot(portSide), routingSlotPositions);
        
        // Now take care of the port anchor
        KVector anchor = lPort.getPosition().clone().add(lPort.getAnchor());
        switch (lPort.getSide()) {
        case NORTH:
            result.x += anchor.x;
            break;
        case SOUTH:
            result.x += anchor.x;
            break;
            
        case EAST:
            result.y += anchor.y;
            break;
        case WEST:
            result.y += anchor.y;
            break;
            
        default:
            assert false;
        }
        
        bendPoints.add(result);
    }

    private void addCornerBendPoints(final SelfLoopEdge slEdge, final EdgeRoutingDirection routingDirection,
            final double[][] routingSlotPositions, final KVectorChain bendPoints) {
        
        // Check if we even need corner bend points
        LPort lSourcePort = slEdge.getSLSource().getLPort();
        LPort lTargetPort = slEdge.getSLTarget().getLPort();
        
        if (lSourcePort.getSide() == lTargetPort.getSide()) {
            return;
        }
        
        SelfHyperLoop slLoop = slEdge.getSLHyperLoop();
        // Use inline labels side and its size to determine the bendpoint offset required to place the label.
        PortSide labelSide = null;
        KVector lSize = null;
        boolean inline = slEdge.isInline();
        if (inline && slLoop.getSLLabels() != null) {
            labelSide = slEdge.getLabelSide();
            lSize = slLoop.getSLLabels().getSize();
        }
        
        // Compute corner points
        PortSide currPortSide = lSourcePort.getSide();
        PortSide nextPortSide = null;
        
        while (currPortSide != lTargetPort.getSide()) {
            // Next port side depends on the direction we're going
            nextPortSide = routingDirection == EdgeRoutingDirection.CLOCKWISE
                    ? currPortSide.right()
                    : currPortSide.left();
            
            // Compute the coordinates contributes by the current and next port sides
            KVector currPortSideComponent = getBaseVector(
                    currPortSide, slLoop.getRoutingSlot(currPortSide), routingSlotPositions);
            KVector nextPortSideComponent = getBaseVector(
                    nextPortSide, slLoop.getRoutingSlot(nextPortSide), routingSlotPositions);
            // If the label is inline, we need to reserve space for each side to accommodate it.
            if (inline && labelSide != null && lSize != null) {
                if (currPortSide == labelSide) {
                    adjustVectorForLabelSide(currPortSideComponent, labelSide, lSize);
                } else if (nextPortSide == labelSide) {
                    adjustVectorForLabelSide(nextPortSideComponent, labelSide, lSize);
                }
            }
            
            // One has its x coordinate set, the other has its y coordinate set -- their sum is our final bend point
            bendPoints.add(currPortSideComponent.add(nextPortSideComponent));
            
            // Advance to next port side
            currPortSide = nextPortSide;
        }
    }
    
    private KVector getBaseVector(final PortSide portSide, final int routingSlot,
            final double[][] routingSlotPositions) {
        
        double position = routingSlotPositions[portSide.ordinal()][routingSlot];
        
        switch (portSide) {
        case NORTH:
        case SOUTH:
            return new KVector(0, position);
            
        case EAST:
        case WEST:
            return new KVector(position, 0);
            
        default:
            assert false;
            return null;
        }
    }
    
    /**
     * Adjusts the given vector to account for the label side and size. This is necessary to ensure that the
     * label is centered on the bend point when it is inline with the edge.
     * 
     * @param portSideComponent The port side component vector that will be adjusted.
     * @param labelSide The side on which the label is placed.
     * @param labelSize The size of the label.
     */
    private void adjustVectorForLabelSide(KVector portSideComponent, PortSide labelSide, KVector labelSize) {
        switch (labelSide) {
        case NORTH:
            portSideComponent.y -= labelSize.y / 2;
            break;
        case SOUTH:
            portSideComponent.y += labelSize.y / 2;
            break;
        case WEST:
            portSideComponent.x -= labelSize.x / 2;
            break;
        case EAST:
            portSideComponent.x += labelSize.x / 2;
            break;
        }
    }

}
