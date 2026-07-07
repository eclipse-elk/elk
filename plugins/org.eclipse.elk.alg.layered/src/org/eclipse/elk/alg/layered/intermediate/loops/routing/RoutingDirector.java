/*******************************************************************************
 * Copyright (c) 2019 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.layered.intermediate.loops.routing;

import java.util.Comparator;
import java.util.List;
import java.util.Set;

import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.intermediate.loops.SelfHyperLoop;
import org.eclipse.elk.alg.layered.intermediate.loops.SelfLoopHolder;
import org.eclipse.elk.alg.layered.intermediate.loops.SelfLoopPort;
import org.eclipse.elk.alg.layered.intermediate.loops.ordering.PortRestorer;
import org.eclipse.elk.core.options.PortSide;

/**
 * Decides for each {@link SelfHyperLoop} which of its {@link SelfLoopPort}s is the leftmost and which is the rightmost
 * one. This effectively decides how the loop is routed around the node. How complex this is to determine depends on the
 * self loop types (see the method documentations for how we do so).
 */
public class RoutingDirector {

    /** A comparator for {@link SelfLoopPort}s based on the {@link LPort} IDs. */
    private static final Comparator<SelfLoopPort> COMPARE_BY_ID =
            (slPort1, slPort2) -> Integer.compare(slPort1.getLPort().id, slPort2.getLPort().id);
            
    /** The penalty an edge incurs for passing a port without any connections. */
    private static final int UNCONNECTED_PORT_PENALTY = 1;
    /** The penalty an edge incurs for passing a port with connections. */
    private static final int CONNECTED_PORT_PENALTY = 3;

    /**
     * The main data structure we use to compute penalties incurred by connecting two {@link SelfLoopPort}s. Each port
     * has an ID which we use to index into this array. The entry then specifies the penalty an edge that goes through
     * all the ports (from the top left port through the given port) would incur. We can compute penalties for edges
     * between ports with IDs {@code idLeft} and {@code idRight} would incur. If {@code idLeft < idRight}, we simply
     * compute {@code portPenalties[idRight - 1] - portPenalties[idLeft]}. If {@code idLeft > idRight}, we wrap around
     * the top left corner, thus resetting the port IDs back to zero. But that case is not hard either.
     */
    private int[] portPenalties = null;

    /**
     * Sets the leftmost and rightmost ports of each {@link SelfHyperLoop} in the given {@link SelfLoopHolder}.
     */
    public void determineLoopRoutes(final SelfLoopHolder slHolder) {
        // Start by giving IDs to all ports according to their order in the port list
        assignPortIds(slHolder.getLNode().getPorts());
        sortHyperLoopPortLists(slHolder);

        // Now assign stuff! Preferrably, assign leftmost and rightmost ports...
        for (SelfHyperLoop slLoop : slHolder.getSLHyperLoops()) {
            switch (slLoop.getSelfLoopType()) {
            case ONE_SIDE:
                determineOneSideLoopRoutes(slLoop);
                break;

            case TWO_SIDES_CORNER:
                determineTwoSideCornerLoopRoutes(slLoop);
                break;

            case TWO_SIDES_OPPOSING:
                determineTwoSideOpposingLoopRoutes(slLoop);
                break;

            case THREE_SIDES:
                determineThreeSideLoopRoutes(slLoop);
                break;

            case FOUR_SIDES:
                determineFourSideLoopRoutes(slLoop);
                break;

            default:
                assert false;
            }
            
            // Setup the loop's set of occupied port sides
            computeOccupiedPortSides(slLoop);
        }
        
        portPenalties = null;
    }

    private void assignPortIds(final List<LPort> lPorts) {
        for (int i = 0; i < lPorts.size(); i++) {
            lPorts.get(i).id = i;
        }
    }

    private void sortHyperLoopPortLists(final SelfLoopHolder slHolder) {
        slHolder.getSLHyperLoops().stream()
            .map(slLoop -> slLoop.getSLPorts())
            .forEach(slPorts -> slPorts.sort(COMPARE_BY_ID));
    }

    private void computeOccupiedPortSides(final SelfHyperLoop slLoop) {
        PortSide currPortSide = slLoop.getLeftmostPort().getLPort().getSide();
        PortSide targetSide = slLoop.getRightmostPort().getLPort().getSide();
        
        while (currPortSide != targetSide) {
            slLoop.getOccupiedPortSides().add(currPortSide);
            currPortSide = currPortSide.right();
        }
        
        slLoop.getOccupiedPortSides().add(currPortSide);
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Route Determination (One Side)

    /**
     * One-sided self loops should always be routed within their one side, so that's easy.
     */
    private void determineOneSideLoopRoutes(final SelfHyperLoop slLoop) {
        PortSide side = slLoop.getSLPorts().get(0).getLPort().getSide();
        assignLeftmostRightmostPorts(slLoop, side, side);
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Route Determination (Two Sides, corner)

    /**
     * Two-sided corner self loops are easy as well since we only want to span that corner.
     */
    private void determineTwoSideCornerLoopRoutes(final SelfHyperLoop slLoop) {
        PortSide[] sides = PortRestorer.sortedTwoSideLoopPortSides(slLoop);
        assignLeftmostRightmostPorts(slLoop, sides[0], sides[1]);
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Route Determination (One Side, opposing)

    /**
     * Opposing side loops have two ways to be routed. Those can be identified by the otherwise uninvolved port side the
     * loop is routed through. The exact routing we choose determines the number of edge crossings.
     */
    private void determineTwoSideOpposingLoopRoutes(final SelfHyperLoop slLoop) {
        PortSide[] sides = slLoop.getSLPortsBySide().keySet().toArray(new PortSide[2]);
        assert sides.length == 2;

        SelfLoopHolder slHolder = slLoop.getSLHolder();

        // We basically play through both possible options and use the one which has a lower penalty
        SelfLoopPort option1LeftmostPort = lowestPortOnSide(slLoop, sides[0]);
        SelfLoopPort option1RightmostPort = highestPortOnSide(slLoop, sides[1]);
        int option1Penalty = computeEdgePenalty(slHolder, option1LeftmostPort, option1RightmostPort);

        SelfLoopPort option2LeftmostPort = lowestPortOnSide(slLoop, sides[1]);
        SelfLoopPort option2RightmostPort = highestPortOnSide(slLoop, sides[0]);
        int option2Penalty = computeEdgePenalty(slHolder, option2LeftmostPort, option2RightmostPort);

        if (option1Penalty <= option2Penalty) {
            slLoop.setLeftmostPort(option1LeftmostPort);
            slLoop.setRightmostPort(option1RightmostPort);
        } else {
            slLoop.setLeftmostPort(option2LeftmostPort);
            slLoop.setRightmostPort(option2RightmostPort);
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Route Determination (Three Sides)

    /**
     * Three-sided self loops only have one possible way of being routed, so that's easy.
     */
    private void determineThreeSideLoopRoutes(final SelfHyperLoop slLoop) {
        PortSide leftmostSide = null;
        PortSide rightmostSide = null;

        // Determine the leftmost and rightmost of the three port sides (we cannot simply sort the port sides since
        // loops can span the WEST-NORTH corner)
        switch (computeMissingPortSide(slLoop)) {
        case NORTH:
            leftmostSide = PortSide.EAST;
            rightmostSide = PortSide.WEST;
            break;

        case EAST:
            leftmostSide = PortSide.SOUTH;
            rightmostSide = PortSide.NORTH;
            break;

        case SOUTH:
            leftmostSide = PortSide.WEST;
            rightmostSide = PortSide.EAST;
            break;

        case WEST:
            leftmostSide = PortSide.NORTH;
            rightmostSide = PortSide.SOUTH;
            break;

        default:
            assert false;
        }

        assignLeftmostRightmostPorts(slLoop, leftmostSide, rightmostSide);
    }

    private PortSide computeMissingPortSide(final SelfHyperLoop slLoop) {
        Set<PortSide> sides = slLoop.getSLPortsBySide().keySet();

        // SUPPRESS CHECKSTYLE NEXT MagicNumber
        assert sides.size() == 3;
        assert !sides.contains(PortSide.UNDEFINED);

        for (PortSide side : PortSide.values()) {
            if (side != PortSide.UNDEFINED && !sides.contains(side)) {
                return side;
            }
        }

        // This cannot happen
        assert false;
        return null;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Route Determination (Four Sides)

    /**
     * For four-sided self loops, we have at least one port on every side. If we look at the list of ports, this means
     * that we want to insert a break between some pair of ports adjacent in that list: one will become the leftmost
     * port, and the other will become the rightmost port. What we want to do is find the pair of ports with the maximum
     * number of ports with external connections between them, since that will minimize the number of edge crossings our
     * routing will produce.
     */
    private void determineFourSideLoopRoutes(final SelfHyperLoop slLoop) {
        // The self loop ports are sorted by ID
        List<SelfLoopPort> sortedSLPorts = slLoop.getSLPorts();
        SelfLoopHolder slHolder = slLoop.getSLHolder();

        // Go through each pair of adjacent ports and find the one which incurs the highest penalty if we drew an edge
        // between them. That's the pair we want to split the loop at. We start with the uppermost port on the western
        // side and the leftmost on the northern side and then compare successive pairs against those two
        SelfLoopPort worstLeftPort = sortedSLPorts.get(sortedSLPorts.size() - 1);
        SelfLoopPort worstRightPort = sortedSLPorts.get(0);
        int worstPenalty = computeEdgePenalty(slHolder, worstLeftPort, worstRightPort);

        for (int rightPortIndex = 1; rightPortIndex < sortedSLPorts.size(); rightPortIndex++) {
            SelfLoopPort currLeftPort = sortedSLPorts.get(rightPortIndex - 1);
            SelfLoopPort currRightPort = sortedSLPorts.get(rightPortIndex);
            int currPenalty = computeEdgePenalty(slHolder, currLeftPort, currRightPort);

            if (currPenalty > worstPenalty) {
                worstLeftPort = currLeftPort;
                worstRightPort = currRightPort;
                worstPenalty = currPenalty;
            }
        }

        // Since we _don't_ want to draw the self loop between the left and right ports, we switch them here
        slLoop.setLeftmostPort(worstRightPort);
        slLoop.setRightmostPort(worstLeftPort);
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Utility Methods

    /**
     * Assigns's the loop's leftmost and rightmost ports by taking the lowest-ID port on the given leftmost side and the
     * highest-ID port on the given rightmost side.
     */
    private void assignLeftmostRightmostPorts(final SelfHyperLoop slLoop, final PortSide leftmostSide,
            final PortSide rightmostSide) {

        slLoop.setLeftmostPort(lowestPortOnSide(slLoop, leftmostSide));
        slLoop.setRightmostPort(highestPortOnSide(slLoop, rightmostSide));
    }

    /**
     * Returns the port with the lowest ID on the given side.
     */
    private SelfLoopPort lowestPortOnSide(final SelfHyperLoop slLoop, final PortSide side) {
        return slLoop.getSLPortsBySide(side).stream().min(COMPARE_BY_ID).get();
    }

    /**
     * Returns the port with the highest ID on the given side.
     */
    private SelfLoopPort highestPortOnSide(final SelfHyperLoop slLoop, final PortSide side) {
        return slLoop.getSLPortsBySide(side).stream().max(COMPARE_BY_ID).get();
    }

    /**
     * Computes the penalty incurred by an edge running from the leftmost port clockwise to the rightmost port.
     */
    private int computeEdgePenalty(final SelfLoopHolder slHolder, final SelfLoopPort leftmostPort,
            final SelfLoopPort rightmostPort) {

        // Compute penalties on demand
        if (portPenalties == null) {
            computePenalties(slHolder);
        }

        int portCount = slHolder.getLNode().getPorts().size();
                
        int leftmostPortId = leftmostPort.getLPort().id;
        int rightmostPortId = rightmostPort.getLPort().id;
        int leftOfRightmostPortId = rightmostPortId - 1;
        
        // We can have two cases here: either, leftmostPortId < rightmostPortId. Then leftOfRightmostPortId will be > 0
        // and we're happy. If leftmostPortId > rightmostPortId, it can happen that rightmostPortId == 0. Then we need
        // to adjust leftOfRightmostPortId to be a valid index
        if (leftOfRightmostPortId < 0) {
            leftOfRightmostPortId = portCount - 1;
        }

        if (leftmostPortId <= leftOfRightmostPortId) {
            // This can be computed directly
            return portPenalties[leftOfRightmostPortId] - portPenalties[leftmostPortId];
        
        } else {
            // Our edge from the leftmost to the rightmost port would pass the top left corner, where indices reset to
            // zero. We need to split the result into a sum of partial results (from left to corner, from corner to
            // right)
            return portPenalties[portCount - 1] - portPenalties[leftmostPortId] +  portPenalties[leftOfRightmostPortId];
        }
    }
    
    /**
     * Computes accumulated port penalties and stores them in {@link #portPenalties}. Connected and unconnected ports
     * incur different penalties. We could think about having self loop ports not incur a penalty since their edges are
     * not guaranteed to cross other self loop edges, but I leave that as future work.
     */
    private void computePenalties(final SelfLoopHolder slHolder) {
        // Simply iterate over the ports and add to the penalty sum
        List<LPort> ports = slHolder.getLNode().getPorts();
        portPenalties = new int[ports.size()];
        int penaltySum = 0;
        
        for (int i = 0; i < ports.size(); i++) {
            LPort currPort = ports.get(i);
            if (currPort.getIncomingEdges().isEmpty() && currPort.getOutgoingEdges().isEmpty()) {
                penaltySum += UNCONNECTED_PORT_PENALTY;
            } else {
                penaltySum += CONNECTED_PORT_PENALTY;
            }
            
            portPenalties[i] = penaltySum;
        }
    }

}
