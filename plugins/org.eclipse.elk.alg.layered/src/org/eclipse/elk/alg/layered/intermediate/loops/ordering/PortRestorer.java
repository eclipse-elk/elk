/*******************************************************************************
 * Copyright (c) 2019 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.alg.layered.intermediate.loops.ordering;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.EnumSet;
import java.util.List;
import java.util.Set;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.intermediate.loops.SelfHyperLoop;
import org.eclipse.elk.alg.layered.intermediate.loops.SelfLoopHolder;
import org.eclipse.elk.alg.layered.intermediate.loops.SelfLoopPort;
import org.eclipse.elk.alg.layered.intermediate.loops.SelfLoopType;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.layered.options.SelfLoopOrderingStrategy;
import org.eclipse.elk.core.options.PortSide;
import org.eclipse.elk.core.util.IElkProgressMonitor;

import java.util.EnumMap;
import java.util.HashMap;
import java.util.Map;

/**
 * Restores hidden ports assuming fixed side port constraints.
 */
public class PortRestorer {
    
    /** {@link SelfHyperLoop}s indexed by {@link SelfLoopType}. */
    private Map<SelfLoopType, List<SelfHyperLoop>> slLoopsByType;
    /** A table of all of the different areas around the node to place self loop ports. */
    private Map<PortSide, Map<PortSideArea, List<SelfLoopPort>>> targetAreas;

    /**
     * Restores all previously hidden ports at proper locations. This will also reset all hiding-related properties of
     * both the {@link SelfLoopPort}s and the {@link SelfLoopHolder}.
     */
    public void restorePorts(final SelfLoopHolder slHolder, final IElkProgressMonitor monitor) {
        // Initialize all of the port lists we'll be working on before finally merging them into the node's list of
        // ports
        initTargetAreas();
        
        // We distinguish different types of loops depending (mainly) on their number of port sides. They need to be
        // collected first
        slLoopsByType = gatherSelfLoopsByType(slHolder);

        // Now the loops have to be added in a certain order in an attempt to minimize edge crossings. This is not
        // exactly perfect since it doesn't take into account regular eastern / western ports that may be part of self
        // loops as well, but it's at least comparatively easy
        processOneSideLoops(slHolder.getLNode().getProperty(LayeredOptions.EDGE_ROUTING_SELF_LOOP_ORDERING));
        processTwoSideCornerLoops();
        processThreeSideLoops();
        processFourSideLoops();
        processTwoSideOpposingLoops();
        
        // Actually go ahead and add the ports to the real port list
        restorePorts(slHolder);

        // We're not hiding any ports anymore
        targetAreas.values().stream()
            .flatMap(byArea -> byArea.values().stream())
            .flatMap(slPortList -> slPortList.stream())
            .forEach(slPort -> slPort.setHidden(false));
        slHolder.setPortsHidden(false);
        
        // Reset
        slLoopsByType = null;
    }
    
    /**
     * Initializes our table of target areas. Once this method has finished, there is a port list for each combination
     * of port side and port side area in {@link #targetAreas}.
     */
    private void initTargetAreas() {
        targetAreas = new EnumMap<>(PortSide.class);
        for (PortSide side : PortSide.values()) {
            Map<PortSideArea, List<SelfLoopPort>> areaMap = new EnumMap<>(PortSideArea.class);
            for (PortSideArea area : PortSideArea.values()) {
                areaMap.put(area, new ArrayList<>());
            }
            targetAreas.put(side, areaMap);
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Self Loop Gathering

    /**
     * Gathers all of the self loops and distributes them by type.
     */
    private Map<SelfLoopType, List<SelfHyperLoop>> gatherSelfLoopsByType(final SelfLoopHolder slHolder) {
        Map<SelfLoopType, List<SelfHyperLoop>> loops = new HashMap<>();
        slHolder.getSLHyperLoops()
                .forEach(slLoop -> loops.computeIfAbsent(slLoop.getSelfLoopType(), k -> new ArrayList<>())
                        .add(slLoop));
        return loops;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Self Loop Placement (One Side)

    private void processOneSideLoops(final SelfLoopOrderingStrategy ordering) {
        if (ordering == SelfLoopOrderingStrategy.REVERSE_STACKED) {
            Collections.reverse(slLoopsByType.getOrDefault(SelfLoopType.ONE_SIDE, Collections.emptyList()));
        }
        for (SelfHyperLoop slLoop : slLoopsByType.getOrDefault(SelfLoopType.ONE_SIDE, Collections.emptyList())) {
            // Obtain the port side
            PortSide side = slLoop.getSLPorts().get(0).getLPort().getSide();
            
            // We want ports with more outgoing edges to be to the left of ports with more incoming edges, so we need
            // a sorted list of ports
            List<SelfLoopPort> sortedPorts = new ArrayList<>(slLoop.getSLPorts());
            sortedPorts.sort((slPort1, slPort2) -> Integer.compare(slPort1.getSLNetFlow(), slPort2.getSLNetFlow()));
            
            switch (ordering) {
            case SEQUENCED:
                // Simply add all ports according to our list
                addToTargetArea(sortedPorts, side, PortSideArea.MIDDLE, AddMode.APPEND);
                break;
            case REVERSE_STACKED:
            case STACKED:
                // Compute which ports we want to have in the first group and which in the second group
                int splitIndex = computePortListSplitIndex(sortedPorts);
                
                // Prepend the first group to the middle list, and append the second group to that same list
                addToTargetArea(sortedPorts.subList(0, splitIndex), side, PortSideArea.MIDDLE, AddMode.PREPEND);
                addToTargetArea(
                        sortedPorts.subList(splitIndex, sortedPorts.size()), side, PortSideArea.MIDDLE, AddMode.APPEND);
                break;
                
            default:
                assert false;
            }
        }
    }

    /**
     * Given a port list sorted by net flow, this method computes a good split index. The split index splits the list of
     * ports in two halves, with the port at the split index being the start of the second half. That port is either the
     * first port with a positive net flow, the first port with a non-negative net flow, or the port at the center of
     * the list. What is will never be is the first port.
     */
    private int computePortListSplitIndex(final List<SelfLoopPort> sortedPorts) {
        // Find index of the first port with a positive net flow
        int positiveNetFlowIndex = 0;
        for (; positiveNetFlowIndex < sortedPorts.size(); positiveNetFlowIndex++) {
            if (sortedPorts.get(positiveNetFlowIndex).getSLNetFlow() > 0) {
                break;
            }
        }
        
        // If this is neither the first, nor the last port, return its index
        if (positiveNetFlowIndex > 0 && positiveNetFlowIndex < sortedPorts.size() - 1) {
            return positiveNetFlowIndex;
        }
        
        // Find index of the first port with a non-negative net flow
        int nonNegativeNetFlowIndex = 0;
        for (; nonNegativeNetFlowIndex < sortedPorts.size(); nonNegativeNetFlowIndex++) {
            if (sortedPorts.get(nonNegativeNetFlowIndex).getSLNetFlow() > 0) {
                break;
            }
        }
        
        // If this is neither the first, nor the last port, return its index
        if (nonNegativeNetFlowIndex > 0 && positiveNetFlowIndex < sortedPorts.size() - 1) {
            return nonNegativeNetFlowIndex;
        }
        
        // We tried our best; simply return the center port's index
        return sortedPorts.size() / 2;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Self Loop Placement (Two Sides)

    private void processTwoSideCornerLoops() {
        for (SelfHyperLoop slLoop : slLoopsByType.getOrDefault(SelfLoopType.TWO_SIDES_CORNER, Collections.emptyList())) {
            // Sort the port sides such that they follow a clockwise order
            PortSide[] sides = sortedTwoSideLoopPortSides(slLoop);
            
            // Add the ports to their target area
            addToTargetArea(slLoop, sides[0], PortSideArea.END, AddMode.PREPEND);
            addToTargetArea(slLoop, sides[1], PortSideArea.START, AddMode.APPEND);
        }
    }

    private void processTwoSideOpposingLoops() {
        for (SelfHyperLoop slLoop : slLoopsByType.getOrDefault(SelfLoopType.TWO_SIDES_OPPOSING, Collections.emptyList())) {
            // Sort the port sides such that they follow a clockwise order
            PortSide[] sides = sortedTwoSideLoopPortSides(slLoop);
            
            // We prepend to the start side's end area, and append to the target side's start area 
            addToTargetArea(slLoop, sides[0], PortSideArea.END, AddMode.PREPEND);
            addToTargetArea(slLoop, sides[1], PortSideArea.START, AddMode.APPEND);
        }
    }
    
    /**
     * Returns the port sides spanned by a two-side self loop in a clockwise order.
     */
    public static PortSide[] sortedTwoSideLoopPortSides(final SelfHyperLoop slLoop) {
        // Sort the port sides such that they follow a clockwise order
        PortSide[] sides = slLoop.getSLPortsBySide().keySet().toArray(new PortSide[2]);
        Arrays.sort(sides);
        
        assert sides.length == 2;
        
        // NORTH and WEST are, of course, the exception and need to be switched...
        if (sides[0] == PortSide.NORTH && sides[1] == PortSide.WEST) {
            sides[0] = PortSide.WEST;
            sides[1] = PortSide.NORTH;
        }
        
        return sides;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Self Loop Placement (Three Sides)
    
    // Different possible constellations of three-sided loops
    private static final PortSide[] NES = { PortSide.NORTH, PortSide.EAST, PortSide.SOUTH };
    private static final PortSide[] ESW = { PortSide.EAST, PortSide.SOUTH, PortSide.WEST };
    private static final PortSide[] SWN = { PortSide.SOUTH, PortSide.WEST, PortSide.NORTH };
    private static final PortSide[] WNE = { PortSide.WEST, PortSide.NORTH, PortSide.EAST };

    private void processThreeSideLoops() {
        for (SelfHyperLoop slLoop : slLoopsByType.getOrDefault(SelfLoopType.THREE_SIDES, Collections.emptyList())) {
            // This array will yield the loop's start, middle, and end sides
            PortSide[] sides = determineLoopConstellation(slLoop);
            
            // Prepend to the start area, append to the other areas
            addToTargetArea(slLoop, sides[0], PortSideArea.END, AddMode.PREPEND);
            addToTargetArea(slLoop, sides[1], PortSideArea.MIDDLE, AddMode.APPEND);
            addToTargetArea(slLoop, sides[2], PortSideArea.START, AddMode.APPEND);
        }
    }
    
    private PortSide[] determineLoopConstellation(final SelfHyperLoop slLoop) {
        Set<PortSide> portSides = slLoop.getSLPortsBySide().keySet();
        
        if (!portSides.contains(PortSide.NORTH)) {
            return ESW;
        } else if (!portSides.contains(PortSide.EAST)) {
            return SWN;
        } else if (!portSides.contains(PortSide.SOUTH)) {
            return WNE;
        } else if (!portSides.contains(PortSide.WEST)) {
            return NES;
        } else {
            assert false;
            return null;
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Self Loop Placement (Four Sides)

    private void processFourSideLoops() {
        for (SelfHyperLoop slLoop : slLoopsByType.getOrDefault(SelfLoopType.FOUR_SIDES, Collections.emptyList())) {
            // Simply append to all port side's middle areas
            // Prepend to the start area, append to the other areas
            for (PortSide side : slLoop.getSLPortsBySide().keySet()) {
                addToTargetArea(slLoop, side, PortSideArea.MIDDLE, AddMode.APPEND);
            }
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Self Loop Placement (Utilities)

    /**
     * Adds ports of the given loop to one of our several target areas.
     * 
     * @param slLoop
     *            the loop whose ports to place.
     * @param portSide
     *            the side of the ports to be placed, which is also the side of the target area we're going to place
     *            them in.
     * @param area
     *            the target area of the side we're going to place them in.
     * @param addMode
     *            whether to append or prepend the ports to the target area's port list.
     */
    private void addToTargetArea(final SelfHyperLoop slLoop, final PortSide portSide, final PortSideArea area,
            final AddMode addMode) {
        
        addToTargetArea(slLoop.getSLPortsBySide(portSide), portSide, area, addMode);
    }
    
    /**
     * Adds a collection of ports to one of our several target areas.
     * 
     * @param ports
     *            the ports to place.
     * @param portSide
     *            the side of the ports to be placed, which is also the side of the target area we're going to place
     *            them in.
     * @param area
     *            the target area of the side we're going to place them in.
     * @param addMode
     *            whether to append or prepend the ports to the target area's port list.
     */
    private void addToTargetArea(final Collection<SelfLoopPort> slPorts, final PortSide portSide,
            final PortSideArea area, final AddMode addMode) {

        // Gather those ports that are currently hidden (if they're not hidden, there's no point restoring them)
        List<SelfLoopPort> hiddenPorts = slPorts.stream()
                .filter(slPort -> slPort.isHidden())
                .collect(Collectors.toList());

        Collections.reverse(hiddenPorts);
        List<SelfLoopPort> targetArea = targetAreas.get(portSide).get(area);
        if (addMode == AddMode.PREPEND) {
            targetArea.addAll(0, hiddenPorts);
        } else {
            targetArea.addAll(hiddenPorts);
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Port Restoring

    private void restorePorts(final SelfLoopHolder slHolder) {
        // We're building up a new port list and replace the old one once we're finished. That's more efficient than
        // continually inserting ports into the middle of the existing list.
        LNode lNode = slHolder.getLNode();
        
        // We'll add the old ports in bursts and always remember where the next burst starts
        List<LPort> oldPortList = new ArrayList<>(lNode.getPorts());
        int nextOldPortIndex = 0;
        
        List<LPort> newPortList = lNode.getPorts();
        newPortList.clear();
        
        // Go over the target areas and add them in between the regular ports
        addAll(targetAreas.get(PortSide.NORTH).get(PortSideArea.START), lNode);
        nextOldPortIndex = addAllThat(
                oldPortList,
                nextOldPortIndex,
                (lPort) -> lPort.getSide() == PortSide.NORTH && isNorthSouthPortWithWestOrWestEastConnections(lPort),
                newPortList);
        addAll(targetAreas.get(PortSide.NORTH).get(PortSideArea.MIDDLE), lNode);
        nextOldPortIndex = addAllThat(
                oldPortList,
                nextOldPortIndex,
                (lPort) -> lPort.getSide() == PortSide.NORTH,
                newPortList);
        addAll(targetAreas.get(PortSide.NORTH).get(PortSideArea.END), lNode);
        
        addAll(targetAreas.get(PortSide.EAST).get(PortSideArea.START), lNode);
        addAll(targetAreas.get(PortSide.EAST).get(PortSideArea.MIDDLE), lNode);
        nextOldPortIndex = addAllThat(
                oldPortList,
                nextOldPortIndex,
                (lPort) -> lPort.getSide() == PortSide.EAST,
                newPortList);
        addAll(targetAreas.get(PortSide.EAST).get(PortSideArea.END), lNode);
        
        addAll(targetAreas.get(PortSide.SOUTH).get(PortSideArea.START), lNode);
        nextOldPortIndex = addAllThat(
                oldPortList,
                nextOldPortIndex,
                (lPort) -> lPort.getSide() == PortSide.SOUTH && isNorthSouthPortWithEastConnections(lPort),
                newPortList);
        addAll(targetAreas.get(PortSide.SOUTH).get(PortSideArea.MIDDLE), lNode);
        nextOldPortIndex = addAllThat(
                oldPortList,
                nextOldPortIndex,
                (lPort) -> lPort.getSide() == PortSide.SOUTH,
                newPortList);
        addAll(targetAreas.get(PortSide.SOUTH).get(PortSideArea.END), lNode);
        
        addAll(targetAreas.get(PortSide.WEST).get(PortSideArea.START), lNode);
        nextOldPortIndex = addAllThat(
                oldPortList,
                nextOldPortIndex,
                (lPort) -> lPort.getSide() == PortSide.WEST,
                newPortList);
        addAll(targetAreas.get(PortSide.WEST).get(PortSideArea.MIDDLE), lNode);
        addAll(targetAreas.get(PortSide.WEST).get(PortSideArea.END), lNode);
        
        assert newPortList.size() >= oldPortList.size() && newPortList.size() >= slHolder.getSLPortMap().size();
        assert !slHolder.arePortsHidden() || newPortList.size() > oldPortList.size();
    }
    
    /**
     * Adds the {@link LPort}s of all of the {@link SelfLoopPort}s to the given node. It does so by calling
     * {@link LPort#setNode(LNode)}, which will automatically add them to the node's port list.
     */
    private void addAll(final List<SelfLoopPort> slPorts, final LNode lNode) {
        slPorts.stream()
            .map(slPort -> slPort.getLPort())
            .forEach(lPort -> lPort.setNode(lNode));
    }
    
    /**
     * Adds all ports from the list to the target list that satisfy the given predicate. This method starts at a given
     * index and does its thing until it finds the first port that does not satisfy the predicate anymore. That port's
     * index is then returned and can be used as the next start index.
     * 
     * @param lPorts
     *            list of ports to add from.
     * @param fromIndex
     *            the index to start from looking for ports that satisfy our predicate.
     * @param condition
     *            the condition ports to be added must satisfy.
     * @param target
     *            the target list to add the ports to.
     * @return index of the first port that does not satisfy the predicate.
     */
    private int addAllThat(final List<LPort> lPorts, final int fromIndex, final Predicate<LPort> condition,
            final List<LPort> target) {
        
        for (int i = fromIndex; i < lPorts.size(); i++) {
            // If this port is valid, add it to the list
            LPort lPort = lPorts.get(i);
            
            if (condition.test(lPort)) {
                target.add(lPort);
            } else {
                return i;
            }
        }
        
        return lPorts.size();
    }
    
    /**
     * Checks whether the given port is a north / south port that has connections to the west, with optional additional
     * connections to the east.
     */
    private boolean isNorthSouthPortWithWestOrWestEastConnections(final LPort lPort) {
        Set<PortSide> connections = northSouthPortConnectionSides(lPort);
        boolean eastConnections = connections.contains(PortSide.EAST);
        boolean westConnections = connections.contains(PortSide.WEST);
        
        return westConnections || (westConnections && eastConnections);
    }
    
    /**
     * Checks whether the given port is a north / south port that only has connections to the east.
     */
    private boolean isNorthSouthPortWithEastConnections(final LPort lPort) {
        Set<PortSide> connections = northSouthPortConnectionSides(lPort);
        return connections.contains(PortSide.EAST);
    }
    
    /**
     * Gathers the sides to which the given north / south port has connections.
     */
    private Set<PortSide> northSouthPortConnectionSides(final LPort lPort) {
        Set<PortSide> connectionSides = EnumSet.noneOf(PortSide.class);
        
        LNode portDummy = lPort.getProperty(InternalProperties.PORT_DUMMY);
        if (portDummy != null) {
            for (LPort dummyLPort : portDummy.getPorts()) {
                if (dummyLPort.getProperty(InternalProperties.ORIGIN) == lPort) {
                    // This dummy port was indeed created for our original port
                    if (dummyLPort.getConnectedEdges().iterator().hasNext()) {
                        // The port has incident edges!
                        connectionSides.add(dummyLPort.getSide());
                    }
                }
            }
        }
        
        return connectionSides;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Enumerations
    
    /**
     * The three different areas of a port side. These are used to place ports of different kinds of self loops as
     * described in the class comment. The areas are interpreted in clockwise order: the beginning of the
     * {@link PortSide#NORTH} is at its left, the beginning of {@link PortSide#EAST} is at its top, etc.
     */
    private static enum PortSideArea {
        /** Beginning of a side. */
        START,
        /** Middle of a side. */
        MIDDLE,
        /** End of a side. */
        END;
    }
    
    /**
     * Whether to prepend or append items to a list.
     */
    private static enum AddMode {
        /** Prepend to a target list. */
        PREPEND,
        /** Append to a target list. */
        APPEND;
    }

}
