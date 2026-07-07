/*******************************************************************************
 * Copyright (c) 2019, 2020 Kiel University and others.
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
package org.eclipse.elk.alg.layered.intermediate.loops.routing;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Queue;
import java.util.Set;
import java.util.stream.Collectors;

import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.intermediate.loops.SelfHyperLoop;
import org.eclipse.elk.alg.layered.intermediate.loops.SelfHyperLoopLabels;
import org.eclipse.elk.alg.layered.intermediate.loops.SelfLoopHolder;
import org.eclipse.elk.alg.layered.intermediate.loops.SelfLoopPort;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.p5edges.orthogonal.HyperEdgeSegment;
import org.eclipse.elk.alg.layered.p5edges.orthogonal.HyperEdgeSegmentDependency;
import org.eclipse.elk.alg.layered.p5edges.orthogonal.OrthogonalRoutingGenerator;
import org.eclipse.elk.core.options.PortSide;

/**
 * Assigns routing slots for each {@link SelfHyperLoop} for each of the port sides it spans. This is done such that it
 * tries to minimize the number of crossings between self loops that ensue. The crossing minimization is based on what
 * the {@link OrthogonalRoutingGenerator} does to order vertical segments of hyperedges.
 */
public class RoutingSlotAssigner {
    
    /** Our segments that constitute our crossing graph. */
    private List<HyperEdgeSegment> hyperEdgeSegments;
    /** Map of segments created for the loops. */
    private Map<SelfHyperLoop, HyperEdgeSegment> slLoopToSegmentMap;
    /** Maps each loop to an array indexed by port IDs and specifying whether the loops runs along the port. */
    private Map<SelfHyperLoop, boolean[]> slLoopActivityOverPorts;

    /**
     * Assigns routing slots to all self loops per port side they span.
     */
    public void assignRoutingSlots(final SelfLoopHolder slHolder) {
        // To be able to check whether labels potentially overlap, we build an overlap matrix to be passed to subsequent
        // phases of the crossing graph
        boolean[][] labelCrossingMatrix = computeLabelCrossingMatrix(slHolder);
        
        // We're using the orthogonal edge router's cycle breaker for this, so create the crossing graph for our loops
        createCrossingGraph(slHolder, labelCrossingMatrix);
        OrthogonalRoutingGenerator.breakNonCriticalCycles(hyperEdgeSegments,
                slHolder.getLNode().getGraph().getProperty(InternalProperties.RANDOM));
        
        // Assign routing slots based on the graph
        doAssignRoutingSlots(slHolder, labelCrossingMatrix);
        
        // Cleanup
        hyperEdgeSegments = null;
        slLoopToSegmentMap = null;
        slLoopActivityOverPorts = null;
    }


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Label Crossing Matrix
    
    /**
     * Computes a two-dimensional array with elements set to {@code true} if the labels with the corresponding IDs
     * overlap.
     */
    private boolean[][] computeLabelCrossingMatrix(final SelfLoopHolder slHolder) {
        // We need to start by giving the labels proper IDs
        int labelID = 0;
        for (SelfHyperLoop slLoop : slHolder.getSLHyperLoops()) {
            if (slLoop.getSLLabels() != null) {
                slLoop.getSLLabels().id = labelID++;
            }
        }
        
        boolean[][] crossingMatrix = new boolean[labelID][labelID];
        
        // Now check for each pair of labels whether or not they overlap
        List<SelfHyperLoop> slLoops = slHolder.getSLHyperLoops();
        for (int sl1Idx = 0; sl1Idx < slLoops.size(); sl1Idx++) {
            SelfHyperLoop slLoop1 = slLoops.get(sl1Idx);
            
            if (slLoop1.getSLLabels() != null) {
                for (int sl2Idx = sl1Idx + 1; sl2Idx < slLoops.size(); sl2Idx++) {
                    SelfHyperLoop slLoop2 = slLoops.get(sl2Idx);
                    
                    if (slLoop2.getSLLabels() != null) {
                        boolean overlap = labelsOverlap(slLoop1, slLoop2);
                        
                        crossingMatrix[slLoop1.getSLLabels().id][slLoop2.getSLLabels().id] = overlap;
                        crossingMatrix[slLoop2.getSLLabels().id][slLoop1.getSLLabels().id] = overlap;
                    }
                }
            }
        }
        
        return crossingMatrix;
    }

    /**
     * Returns {@code true} if the labels of the two loops would overlap if the loops were assigned to the same routing
     * slot.
     */
    private boolean labelsOverlap(final SelfHyperLoop slLoop1, final SelfHyperLoop slLoop2) {
        assert slLoop1 != slLoop2;
        
        SelfHyperLoopLabels slLabels1 = slLoop1.getSLLabels();
        SelfHyperLoopLabels slLabels2 = slLoop2.getSLLabels();
        
        // There won't be overlaps unless both loops have labels
        if (slLabels1 == null || slLabels2 == null) {
            return false;
        }
        
        // The labels must be assigned to the same side, and that side (currently) needs to be either north or south
        if (slLabels1.getSide() != slLabels2.getSide() || slLabels1.getSide() == PortSide.EAST
                || slLabels1.getSide() == PortSide.WEST) {
            
            return false;
        }
        
        // Check if the labels overlap horizontally (remember that at this point, horizontal coordinates of north and
        // south labels have been set by the label placer)
        double start1 = slLabels1.getPosition().x;
        double end1 = start1 + slLabels1.getSize().x;
        double start2 = slLabels2.getPosition().x;
        double end2 = start2 + slLabels2.getSize().x;
        
        return start1 <= end2 && end1 >= start2;
    }
    

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Crossing Graph

    /**
     * Create the crossing graph by creating one {@link HyperEdgeSegment} for each {@link SelfHyperLoop} and computing
     * the proper dependencies between them.
     */
    private void createCrossingGraph(final SelfLoopHolder slHolder, final boolean[][] labelCrossingMatrix) {
        List<SelfHyperLoop> slLoops = slHolder.getSLHyperLoops();
        
        // Create an edge segment for each loop
        hyperEdgeSegments = new ArrayList<>(slLoops.size());
        slLoopToSegmentMap = new HashMap<>();
        
        for (SelfHyperLoop slLoop : slLoops) {
            HyperEdgeSegment segment = new HyperEdgeSegment(null);
            hyperEdgeSegments.add(segment);
            slLoopToSegmentMap.put(slLoop, segment);
        }
        
        // To be able to quickly count crossings later, we remember for each loop whether it's active at a given port
        // ID or not
        slLoopActivityOverPorts = new HashMap<>();
        computeLoopActivity(slHolder);
        
        // For each pair of hyper loops, determine the crossings for one segment routed above the other and vice versa.
        // Note that this will create dependencies between all pairs of self loops, even if they don't share ports. This
        // is intentional. It assigns unique routing slots to each segment, which ensures that there won't be any
        // collisions on sides without ports. The routing space will be compacted at a later step.
        for (int firstIdx = 0; firstIdx < slLoops.size() - 1; firstIdx++) {
            SelfHyperLoop slLoop1 = slHolder.getSLHyperLoops().get(firstIdx);
            for (int secondIdx = firstIdx + 1; secondIdx < slLoops.size(); secondIdx++) {
                createDependencies(slLoop1, slHolder.getSLHyperLoops().get(secondIdx), labelCrossingMatrix);
            }
        }
    }

    /**
     * Initializes {@link #slLoopActivityOverPorts}. Each loop is mapped to an array indexed by port indices indicating
     * whether the loop runs along the given port or not.
     */
    private void computeLoopActivity(final SelfLoopHolder slHolder) {
        List<SelfHyperLoop> slLoops = slHolder.getSLHyperLoops();
        List<LPort> lPorts = slHolder.getLNode().getPorts();
        
        for (SelfHyperLoop slLoop : slLoops) {
            boolean[] loopActivity = new boolean[lPorts.size()];
            slLoopActivityOverPorts.put(slLoop, loopActivity);
            
            // Run from the loop's start port to the end port, possibly wrapping around, and set everything to true
            // along the way
            int lPortIdx = slLoop.getLeftmostPort().getLPort().id - 1;
            int lPortTargetIdx = slLoop.getRightmostPort().getLPort().id;
            
            while (lPortIdx != lPortTargetIdx) {
                lPortIdx = (lPortIdx + 1) % lPorts.size();
                loopActivity[lPortIdx] = true;
            }
        }
    }

    /**
     * Creates the necessary dependencies between the {@link SelfHyperLoop}s with the given indices. The lists must have
     * a one-to-one correspondence between loops and segments, that is, the i-th loop must be represented by the i-th
     * segment. The array says at which ports each loop is active.
     */
    private void createDependencies(final SelfHyperLoop slLoop1, final SelfHyperLoop slLoop2,
            final boolean[][] labelCrossingMatrix) {
        
        // Count the numbers of crossings that would ensue if placing to first segment above the second segment and
        // vice versa
        int firstAboveSecondCrossings = countCrossings(slLoop1, slLoop2);
        int secondAboveFirstCrossings = countCrossings(slLoop2, slLoop1);
        
        // Create dependencies
        HyperEdgeSegment segment1 = slLoopToSegmentMap.get(slLoop1);
        HyperEdgeSegment segment2 = slLoopToSegmentMap.get(slLoop2);
        
        if (firstAboveSecondCrossings < secondAboveFirstCrossings) {
            // The first loop should be above the second loop
            HyperEdgeSegmentDependency.createAndAddRegular(
                    segment1, segment2, secondAboveFirstCrossings - firstAboveSecondCrossings);
            
        } else if (secondAboveFirstCrossings < firstAboveSecondCrossings) {
            // The second loop should be above the first loop
            HyperEdgeSegmentDependency.createAndAddRegular(
                    segment2, segment1, firstAboveSecondCrossings - secondAboveFirstCrossings);
            
        } else if (firstAboveSecondCrossings != 0 || labelsOverlap(slLoop1, slLoop2, labelCrossingMatrix)) {
            // Either both orders cause the same number of crossings (and at least one), or the labels of the two loops
            // overlap and the loops must thus be forced onto different slots
            HyperEdgeSegmentDependency.createAndAddRegular(segment1, segment2, 0);
            HyperEdgeSegmentDependency.createAndAddRegular(segment2, segment1, 0);
        }
    }

    /**
     * Counts the number of crossings that result between a loop routed through a higher and a loop routed through a
     * lower routing slot. The crossings will happen where the upper loop connects to ports within the routing area
     * of the lower loop. That is, for each of the upper loop's ports, if the lower loop is active there, this will
     * cause a crossing.
     */
    private int countCrossings(final SelfHyperLoop slUpperLoop, final SelfHyperLoop slLowerLoop) {
        boolean[] lowerLoopActivity = slLoopActivityOverPorts.get(slLowerLoop);
        int crossings = 0;
        
        for (SelfLoopPort slPort : slUpperLoop.getSLPorts()) {
            if (lowerLoopActivity[slPort.getLPort().id]) {
                crossings++;
            }
        }
        
        return crossings;
    }

    /**
     * Looks at the label crossing matrix to determine whether the labels of the two loops overlap.
     */
    private boolean labelsOverlap(final SelfHyperLoop slLoop1, final SelfHyperLoop slLoop2,
            final boolean[][] labelCrossingMatrix) {
        
        if (slLoop1.getSLLabels() == null || slLoop2.getSLLabels() == null) {
            return false;
        } else {
            return labelCrossingMatrix[slLoop1.getSLLabels().id][slLoop2.getSLLabels().id];
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Slot Assignment

    /**
     * Based on the acyclic crossing graph, assign routing slots to all combinations of self loops and port sides.
     */
    private void doAssignRoutingSlots(final SelfLoopHolder slHolder, final boolean[][] labelCrossingMatrix) {
        // We first compute raw slots
        assignRawRoutingSlotsToSegments();
        
        // We assign raw routing slots, but try to compact the routing slot assignment afterwards (this might not
        // actually be necessary)
        assignRawRoutingSlotsToLoops(slHolder);
        shiftTowardsNode(slHolder, labelCrossingMatrix);
    }

    /**
     * This will assign routing slots to the segments as if the four port sides were completely dependent. However,
     * they are not. If the northern port side has 6 routing slots, that doesn't mean the eastern side will need 6 as
     * well.
     */
    private void assignRawRoutingSlotsToSegments() {
        Queue<HyperEdgeSegment> sinks = new LinkedList<>();
        
        // Fill our queue of sinks, which will ultimately end up being routed nearest to the node. While we go through
        // these, we also reset the in- and out weights to the number of incoming and outgoing dependencies
        for (HyperEdgeSegment segment : hyperEdgeSegments) {
            segment.setInWeight(segment.getIncomingSegmentDependencies().size());
            segment.setOutWeight(segment.getOutgoingSegmentDependencies().size());
            
            if (segment.getOutWeight() == 0) {
                segment.setRoutingSlot(0);
                sinks.add(segment);
            }
        }
        
        // Assign raw routing slots!
        while (!sinks.isEmpty()) {
            HyperEdgeSegment segment = sinks.poll();
            int nextRoutingSlot = segment.getRoutingSlot() + 1;
            
            for (HyperEdgeSegmentDependency inDependency : segment.getIncomingSegmentDependencies()) {
                HyperEdgeSegment sourceSegment = inDependency.getSource();
                sourceSegment.setRoutingSlot(Math.max(sourceSegment.getRoutingSlot(), nextRoutingSlot));
                
                sourceSegment.setOutWeight(sourceSegment.getOutWeight() - 1);
                if (sourceSegment.getOutWeight() == 0) {
                    sinks.add(sourceSegment);
                }
            }
        }
    }

    /**
     * Assigns the computed raw routing slots to each loop on each side the loop occupies.
     */
    private void assignRawRoutingSlotsToLoops(final SelfLoopHolder slHolder) {
        for (SelfHyperLoop slLoop : slHolder.getSLHyperLoops()) {
            int slot = slLoopToSegmentMap.get(slLoop).getRoutingSlot();
            for (PortSide portSide : slLoop.getOccupiedPortSides()) {
                slLoop.setRoutingSlot(portSide, slot);
            }
        }
    }

    /**
     * Moves the self loops towards the node on each of the node's sides to avoid empty routing slots.
     */
    private void shiftTowardsNode(final SelfLoopHolder slHolder, final boolean[][] labelCrossingMatrix) {
        // For each port, this array specifies the next routing slot we can assign a self loop on that side to
        int[] nextFreeRoutingSlotAtPort = new int[slHolder.getLNode().getPorts().size()];
        
        shiftTowardsNodeOnSide(slHolder, PortSide.NORTH, nextFreeRoutingSlotAtPort, labelCrossingMatrix);
        shiftTowardsNodeOnSide(slHolder, PortSide.EAST, nextFreeRoutingSlotAtPort, labelCrossingMatrix);
        shiftTowardsNodeOnSide(slHolder, PortSide.SOUTH, nextFreeRoutingSlotAtPort, labelCrossingMatrix);
        shiftTowardsNodeOnSide(slHolder, PortSide.WEST, nextFreeRoutingSlotAtPort, labelCrossingMatrix);
    }

    private void shiftTowardsNodeOnSide(final SelfLoopHolder slHolder, final PortSide side,
            final int[] nextFreeRoutingSlotAtPort, final boolean[][] labelCrossingMatrix) {
        
        // We will iterate over the self loops that occupy that port side, sorted ascendingly by routing slot
        List<SelfHyperLoop> slLoops = slHolder.getSLHyperLoops().stream()
                .filter(slLoop -> slLoop.getOccupiedPortSides().contains(side))
                .sorted((sl1, sl2) -> Integer.compare(sl1.getRoutingSlot(side), sl2.getRoutingSlot(side)))
                .collect(Collectors.toList());
        
        // Find the indices of the first and last regular port on the port side
        int minLPortIndex = Integer.MAX_VALUE;
        int maxLPortIndex = Integer.MIN_VALUE;
        for (LPort lPort : slHolder.getLNode().getPorts()) {
            if (lPort.getSide() == side) {
                minLPortIndex = Math.min(minLPortIndex, lPort.id);
                maxLPortIndex = Math.max(maxLPortIndex, lPort.id);
            }
        }
        
        if (minLPortIndex == Integer.MAX_VALUE) {
            // There are no ports on this side, so we simply assign the loops to consecutive slots starting with 0. We
            // won't cause label overlaps.
            for (int i = 0; i < slLoops.size(); i++) {
                slLoops.get(i).setRoutingSlot(side, i);
            }
        
        } else {
            int[] slotAssignedToLabel = new int[labelCrossingMatrix.length];
            Arrays.fill(slotAssignedToLabel, -1);
            
            // There are ports on this side. Find the lowest free slot across all ports our loop spans, and ensure that
            // no label our loop label conflicts with is assigned to that slot
            for (SelfHyperLoop slLoop : slLoops) {
                boolean[] activeAtPort = slLoopActivityOverPorts.get(slLoop);
                int lowestAvailableSlot = 0;
                
                for (int portIndex = minLPortIndex; portIndex <= maxLPortIndex; portIndex++) {
                    if (activeAtPort[portIndex]) {
                        lowestAvailableSlot = Math.max(lowestAvailableSlot, nextFreeRoutingSlotAtPort[portIndex]);
                    }
                }
                
                // If we have a label, it could be that we are in conflict with another label placed at the lowest
                // available slot. There might also be labels at higher slots. We thus need to assemble all slots
                // occupied by conflicting labels that have already been assigned to one and find the first
                // conflict-free slot
                if (slLoop.getSLLabels() != null) {
                    int ourLabelIdx = slLoop.getSLLabels().id;
                    Set<Integer> slotsWithLabelConflicts = new HashSet<>();
                    
                    for (int otherLabelIdx = 0; otherLabelIdx < labelCrossingMatrix.length; otherLabelIdx++) {
                        if (labelCrossingMatrix[ourLabelIdx][otherLabelIdx]) {
                            slotsWithLabelConflicts.add(slotAssignedToLabel[otherLabelIdx]);
                        }
                    }
                    
                    // Find the first slot (starting with out lowest available) that does not appear in the set
                    while (slotsWithLabelConflicts.contains(lowestAvailableSlot)) {
                        lowestAvailableSlot++;
                    }
                }
                
                // Assign the loop to that routing slot and update out routing slot array
                slLoop.setRoutingSlot(side, lowestAvailableSlot);
                for (int portIndex = minLPortIndex; portIndex <= maxLPortIndex; portIndex++) {
                    if (activeAtPort[portIndex]) {
                        nextFreeRoutingSlotAtPort[portIndex] = lowestAvailableSlot + 1;
                    }
                }
                
                // If we have a label, update the label's routing slot
                if (slLoop.getSLLabels() != null) {
                    slotAssignedToLabel[slLoop.getSLLabels().id] = lowestAvailableSlot;
                }
            }
        }
    }

}
