/*******************************************************************************
 * Copyright (c) 2018, 2019 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.alg.layered.intermediate.loops.routing;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import org.eclipse.elk.alg.layered.graph.LLabel;
import org.eclipse.elk.alg.layered.graph.LMargin;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.intermediate.LabelManagementProcessor;
import org.eclipse.elk.alg.layered.intermediate.SelfLoopRouter;
import org.eclipse.elk.alg.layered.intermediate.loops.SelfHyperLoop;
import org.eclipse.elk.alg.layered.intermediate.loops.SelfHyperLoopLabels;
import org.eclipse.elk.alg.layered.intermediate.loops.SelfHyperLoopLabels.Alignment;
import org.eclipse.elk.alg.layered.intermediate.loops.SelfLoopHolder;
import org.eclipse.elk.alg.layered.intermediate.loops.SelfLoopPort;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.layered.options.SelfLoopOrderingStrategy;
import org.eclipse.elk.core.labels.ILabelManager;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.options.PortSide;
import org.eclipse.elk.core.util.IElkProgressMonitor;

/**
 * Determines where self loop labels will be placed and runs label management, if configured. The placement of each
 * label is described by the port side it will end up on as well as its alignment on that side. For labels assigned to
 * the northern and southern sides, the x coordinate will be computed as well. The same goes for the y coordinate of
 * labels assigned to the eastern and western sides. The other coordinate must be computed by the
 * {@link SelfLoopRouter}.
 */
public final class LabelPlacer {

    /**
     * Places labels for the given {@link SelfLoopHolder}.
     */
    public void placeLabels(final SelfLoopHolder slHolder, final ILabelManager labelManager,
            final IElkProgressMonitor monitor) {
        
        assignSideAndAlignment(slHolder);
        
        for (SelfHyperLoop slLoop : slHolder.getSLHyperLoops()) {
            if (slLoop.getSLLabels() != null) {
                if (labelManager != null) {
                    manageLabels(slLoop, labelManager);
                }
                
                computeCoordinates(slLoop);
            }
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Side and Alignment

    private void assignSideAndAlignment(final SelfLoopHolder slHolder) {
        // Most loops are straightforward, but the one-sides northern and southern ones may not be. If the loops are
        // stacked, we can simply center their labels. If they are sequenced, we want to go from the outer loops to the
        // inner loops in pairs, and align their labels away from each other. To do so, we need to collect all of them
        // first.
        List<SelfHyperLoop> northernOneSidedSLLoops = null;
        List<SelfHyperLoop> southernOneSidedSLLoops = null;
        
        SelfLoopOrderingStrategy orderingStrategy = slHolder.getLNode().getProperty(
                LayeredOptions.EDGE_ROUTING_SELF_LOOP_ORDERING);
        if (orderingStrategy == SelfLoopOrderingStrategy.SEQUENCED) {
            northernOneSidedSLLoops = new ArrayList<>();
            southernOneSidedSLLoops = new ArrayList<>();
        }
        
        // Assign sides and alignments; how this works differs for the different kinds of labels
        for (SelfHyperLoop slLoop : slHolder.getSLHyperLoops()) {
            // If this loop doesn't have any labels, don't bother
            SelfHyperLoopLabels slLabels = slLoop.getSLLabels();
            if (slLabels == null) {
                continue;
            }
            
            assert !slLabels.getLLabels().isEmpty();
            
            // How we place labels largely depends on the self loop type
            switch (slLoop.getSelfLoopType()) {
            case ONE_SIDE:
                PortSide loopSide = slLoop.getOccupiedPortSides().iterator().next();
                
                if (orderingStrategy == SelfLoopOrderingStrategy.SEQUENCED && (loopSide == PortSide.NORTH)) {
                    // Collect for deferred processing
                    northernOneSidedSLLoops.add(slLoop);
                } else if (orderingStrategy == SelfLoopOrderingStrategy.SEQUENCED && (loopSide == PortSide.SOUTH)) {
                    // Collect for deferred processing
                    southernOneSidedSLLoops.add(slLoop);
                } else {
                    assignOneSidedSimpleSideAndAlignment(slLoop, loopSide);
                }
                
                break;
                
            case TWO_SIDES_CORNER:
                assignTwoSidesCornerSideAndAlignment(slLoop);
                break;
                
            case TWO_SIDES_OPPOSING:
            case THREE_SIDES:
                assignTwoSidesOpposingAndThreeSidesSideAndAlignment(slLoop);
                break;
                
            case FOUR_SIDES:
                assignFourSidesSideAndAlignment(slLoop);
                break;
                
            default:
                assert false;
            }
        }
        
        // Process deferred loops
        if (northernOneSidedSLLoops != null) {
            if (!northernOneSidedSLLoops.isEmpty()) {
                assignOneSidedSequencedSideAndAlignment(northernOneSidedSLLoops, PortSide.NORTH);
            }
            
            if (!southernOneSidedSLLoops.isEmpty()) {
                assignOneSidedSequencedSideAndAlignment(southernOneSidedSLLoops, PortSide.SOUTH);
            }
        }
    }
    
    /**
     * Side and alignment assignment for those one-sided self loops where the choice is obvious. This is always the
     * case for eastern and western self loops. For northern and southern self loops, this is the case if they are
     * stacked.
     */
    private void assignOneSidedSimpleSideAndAlignment(final SelfHyperLoop slLoop, final PortSide loopSide) {
        // Remove inline edge label property since it is not valid in this case.
        for (LLabel label : slLoop.getSLLabels().getLLabels()) {
            label.setProperty(LayeredOptions.EDGE_LABELS_INLINE, null);
        }
        
        // The alignment depends on the side. For eastern and western loops, we always align to the topmost port.
        // For northern or southern loops, we always center the label.
        switch (loopSide) {
        case EAST:
        case WEST:
            // The alignment will be relative to the topmost port (which must be either the leftmost or the rightmost
            // port)
            SelfLoopPort topmostPort = slLoop.getLeftmostPort();
            if (slLoop.getRightmostPort().getLPort().getPosition().y < topmostPort.getLPort().getPosition().y) {
                topmostPort = slLoop.getRightmostPort();
            }
            
            assignSideAndAlignment(slLoop, loopSide, Alignment.TOP, topmostPort);
            break;
            
        case NORTH:
        case SOUTH:
            assignSideAndAlignment(slLoop, loopSide, Alignment.CENTER, null);
            break;
        
        default:
            assert false;
        }
    }
    
    private void assignOneSidedSequencedSideAndAlignment(final List<SelfHyperLoop> slLoops, final PortSide portSide) {
        assert !slLoops.isEmpty();
        
        // Ensure sensible port IDs
        int id = 0;
        for (LPort lPort : slLoops.get(0).getSLHolder().getLNode().getPorts()) {
            lPort.id = id++;
        }
        
        // We start by sorting our list according to the ID of the leftmost port. For northern loops, this ensures that
        // the list is sorted from left to right. For southern ports, we need to sort descendingly to have the leftmost
        // port end up at list position 0
        if (portSide == PortSide.NORTH) {
            slLoops.sort((slLoop1, slLoop2) -> Integer.compare(
                    slLoop1.getLeftmostPort().getLPort().id, slLoop2.getLeftmostPort().getLPort().id));
        } else {
            slLoops.sort((slLoop1, slLoop2) -> Integer.compare(
                    slLoop2.getLeftmostPort().getLPort().id, slLoop1.getLeftmostPort().getLPort().id));
        }
        
        // Go from outside loops towards inside loops in pairs
        int leftSlLoopIdx = 0;
        int rightSlLoopIdx = slLoops.size() - 1;
        
        while (leftSlLoopIdx < rightSlLoopIdx) {
            SelfHyperLoop leftSlLoop = slLoops.get(leftSlLoopIdx);
            SelfHyperLoop rightSlLoop = slLoops.get(rightSlLoopIdx);
            
            // If the loop is on the northern side, the leftmost port actually is left of the rightmost port. It's
            // flipped for the southern side
            SelfLoopPort leftLoopAlignmentReference = portSide == PortSide.NORTH
                    ? leftSlLoop.getRightmostPort()
                    : leftSlLoop.getLeftmostPort();
            SelfLoopPort rightLoopAlignmentReference = portSide == PortSide.NORTH
                    ? rightSlLoop.getLeftmostPort()
                    : rightSlLoop.getRightmostPort();
            
            assignSideAndAlignment(leftSlLoop, portSide, Alignment.RIGHT, leftLoopAlignmentReference);
            assignSideAndAlignment(rightSlLoop, portSide, Alignment.LEFT, rightLoopAlignmentReference);
            
            leftSlLoopIdx++;
            rightSlLoopIdx--;
        }
        
        // There might be a single loop in the middle
        if (leftSlLoopIdx == rightSlLoopIdx) {
            assignSideAndAlignment(slLoops.get(leftSlLoopIdx), portSide, Alignment.CENTER, null);
        }
    }
    
    private void assignTwoSidesCornerSideAndAlignment(final SelfHyperLoop slLoop) {
        // The side and alignment depends on where the leftmost and rightmost port is. If the leftmost port is north,
        // the loop goes towards the eastern side and we can left-align the label on the nortern side. The other cases
        // are similar.
        PortSide leftmostPortSide = slLoop.getLeftmostPort().getLPort().getSide();
        PortSide rightmostPortSide = slLoop.getRightmostPort().getLPort().getSide();
        
        // Remove inline edge label property since it is not valid in this case.
        for (LLabel label : slLoop.getSLLabels().getLLabels()) {
            label.setProperty(LayeredOptions.EDGE_LABELS_INLINE, null);
        }
        
        if (leftmostPortSide == PortSide.NORTH) {
            assignSideAndAlignment(slLoop, PortSide.NORTH, Alignment.LEFT, slLoop.getLeftmostPort());
        } else if (rightmostPortSide == PortSide.NORTH) {
            assignSideAndAlignment(slLoop, PortSide.NORTH, Alignment.RIGHT, slLoop.getRightmostPort());
        } else if (leftmostPortSide == PortSide.SOUTH) {
            assignSideAndAlignment(slLoop, PortSide.SOUTH, Alignment.RIGHT, slLoop.getLeftmostPort());
        } else if (rightmostPortSide == PortSide.SOUTH) {
            assignSideAndAlignment(slLoop, PortSide.SOUTH, Alignment.LEFT, slLoop.getRightmostPort());
        } else {
            assert false;
        }
    }
    
    private void assignTwoSidesOpposingAndThreeSidesSideAndAlignment(final SelfHyperLoop slLoop) {
        // Where we place our labels depends on the port side not occupied by the self loop. If that is the northern
        // side, for example, the whole southern side is covered and we can center the label there. If it's the western
        // side, on the other hand, we will left-align the label on the northern side.
        Set<PortSide> occupiedSides = slLoop.getOccupiedPortSides();
        // Check whether the hyperloop has inline labels
        boolean hasInlineLabels = false;
        for (LLabel label : slLoop.getSLLabels().getLLabels()) {
            if (label.getProperty(LayeredOptions.EDGE_LABELS_INLINE)) {
                hasInlineLabels = true;
                break;
            }
        }
        
        if (!occupiedSides.contains(PortSide.NORTH)) {
            // This also works for inline edge labels since the label will be "in the middle".
            assignSideAndAlignment(slLoop, PortSide.SOUTH, Alignment.CENTER, null);
        } else if (!occupiedSides.contains(PortSide.SOUTH)) {
            // This also works for inline edge labels since the label will be "in the middle".
            assignSideAndAlignment(slLoop, PortSide.NORTH, Alignment.CENTER, null);
        } else if (!occupiedSides.contains(PortSide.WEST)) {
            // If we have inline labels, we can center them on the eastern side. Otherwise, we left-align them on the
            // northern side, which is the default for edge labels.
            assignSideAndAlignment(slLoop, hasInlineLabels ? PortSide.EAST : PortSide.NORTH,
                    hasInlineLabels ? Alignment.CENTER : Alignment.LEFT,
                    hasInlineLabels ? null : slLoop.getLeftmostPort());
        } else if (!occupiedSides.contains(PortSide.EAST)) {
            // If we have inline labels, we can center them on the western side. Otherwise, we right-align them on the
            // northern side, which is the default for edge labels.
            assignSideAndAlignment(slLoop, hasInlineLabels ? PortSide.WEST : PortSide.NORTH,
                    hasInlineLabels ? Alignment.CENTER : Alignment.RIGHT,
                    hasInlineLabels ? null : slLoop.getRightmostPort());
        } else {
            assert false;
        }
    }
    
    private void assignFourSidesSideAndAlignment(final SelfHyperLoop slLoop) {
        // There are basically just two cases that we need to distinguish: the loop covers the whole northern side, or
        // the loop covers the whole southern side. We always center its label on the side that's completely covered.
        // A side is not completely covered if either the leftmost or the rightmost port is placed on that side, so
        // that's easy to distinguish
        PortSide leftmostPortSide = slLoop.getLeftmostPort().getLPort().getSide();
        PortSide rightmostPortSide = slLoop.getLeftmostPort().getLPort().getSide();
        
        // Remove inline edge label property since it is not valid in this case.
        for (LLabel label : slLoop.getSLLabels().getLLabels()) {
            label.setProperty(LayeredOptions.EDGE_LABELS_INLINE, null);
        }
        
        if (leftmostPortSide == PortSide.NORTH || rightmostPortSide == PortSide.NORTH) {
            assignSideAndAlignment(slLoop, PortSide.SOUTH, Alignment.CENTER, null);
            
        } else {
            assignSideAndAlignment(slLoop, PortSide.NORTH, Alignment.CENTER, null);
        }
    }
    
    private void assignSideAndAlignment(final SelfHyperLoop slLoop, final PortSide side, final Alignment alignment,
            final SelfLoopPort alignmentReference) {
        
        SelfHyperLoopLabels slLabels = slLoop.getSLLabels();
        slLabels.setSide(side);
        slLabels.setAlignment(alignment);
        slLabels.setAlignmentReferenceSLPort(alignmentReference);
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Label Management
    
    private void manageLabels(final SelfHyperLoop slLoop, final ILabelManager labelManager) {
        // How much space we have available depends entirely on our alignment since that also determines whether we're
        // on the north/south or the east/west sides
        SelfHyperLoopLabels slLabels = slLoop.getSLLabels();
        SelfLoopPort alignRef = slLabels.getAlignmentReferenceSLPort();
        
        LNode lNode = slLoop.getSLHolder().getLNode();
        KVector lNodeSize = lNode.getSize();
        LMargin lNodeMargins = slLoop.getSLHolder().getLNode().getMargin();
        
        // targetWidth determines the maximum horizontal space available for label placement.
        // This may be used by LabelManagement to introduce line breaks in the label.
        double targetWidth = 0;
        
        switch (slLabels.getAlignment()) {
        case CENTER:
            // Clip the label to the extent of the node's bounding box
            targetWidth = lNodeMargins.left + lNodeSize.x + lNodeMargins.right;
            break;
            
        case LEFT:
            targetWidth = lNodeSize.x
                    - alignRef.getLPort().getPosition().x - alignRef.getLPort().getAnchor().x
                    + lNodeMargins.right;
            break;
            
        case RIGHT:
            targetWidth = lNodeMargins.left + alignRef.getLPort().getPosition().x + alignRef.getLPort().getAnchor().x;
            break;
            
        case TOP:
            // We have no way of knowing how far to the left / right our self loop will extend on our placement side,
            // so we need a default here
            targetWidth = LabelManagementProcessor.MIN_WIDTH_EDGE_LABELS;
            break;
            
        default:
            assert false;
        }
        
        // Make the label fit the target width or whatever the label manager deems necessary to do.
        slLabels.applyLabelManagement(
                labelManager,
                Math.max(targetWidth, LabelManagementProcessor.MIN_WIDTH_EDGE_LABELS));
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Coordinate Computation

    private void computeCoordinates(final SelfHyperLoop slLoop) {
        // The coordinates really only depend on the alignment since LEFT, CENTER and RIGHT are specific to the north
        // and south sides while TOP is specific to the east and west sides
        SelfHyperLoopLabels slLabels = slLoop.getSLLabels();
        
        SelfLoopPort alignRef = slLabels.getAlignmentReferenceSLPort();
        KVector size = slLabels.getSize();
        KVector pos = slLabels.getPosition();
        
        switch (slLabels.getAlignment()) {
        case CENTER:
            pos.x = (slLoop.getSLHolder().getLNode().getSize().x - size.x) / 2; 
            break;
            
        case LEFT:
            pos.x = alignRef.getLPort().getPosition().x + alignRef.getLPort().getAnchor().x;
            break;
            
        case RIGHT:
            pos.x = alignRef.getLPort().getPosition().x + alignRef.getLPort().getAnchor().x - size.x;
            break;
            
        case TOP:
            pos.y = alignRef.getLPort().getPosition().y + alignRef.getLPort().getAnchor().y;
            break;
            
        default:
            assert false;
        }
        
    }

}