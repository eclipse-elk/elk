/*******************************************************************************
 * Copyright (c) 2010, 2020 Kiel University and others.
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
package org.eclipse.elk.alg.layered.p5edges.orthogonal;

import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.ListIterator;
import java.util.Map;

import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.p5edges.orthogonal.direction.BaseRoutingDirectionStrategy;
import org.eclipse.elk.core.util.Pair;

import com.google.common.collect.Lists;

/**
 * Instances of this class represent the "trunk" of a hyper edge. In left-to-right layouts, this will be the vertical
 * segment of a hyperedge. Such a segment has a list of coordinates for each of its two sides that specifies where
 * incoming or outgoing connections enter and leave the segment, respectively (for horizontal layouts, this will be
 * horizontal edge segments that usually connect to ports). The range of coordinates spanned by a hyper edge segment
 * determines its extent, which can be requested by calling {@link #getStartCoordinate()} and
 * {@link #getEndCoordinate()}.
 * 
 * <p>
 * Whether coordinates refer to the x or y axis depends on the edge routing direction: for horizontal layouts, the edge
 * segment is vertical, which means that the coordinates are y coordinates.
 * </p>
 * 
 * <p>
 * Hyperedge segments can be split by the {@link HyperEdgeSegmentSplitter} if they are part of a cyclic critical
 * dependency (that is, a sequence of edge segments that will cause edge overlaps with each other). Splitting an edge
 * segment will cause involved edges to take a longer detour, but will resolve edge overlaps. The segment keeps two
 * pieces of information about split segments: the new segment introduced to split the segment, and the segment that
 * originally caused the split.
 * </p>
 * 
 * <p>
 * Instances of this class are comparable based on the value of {@link #mark}. {@link #hashCode()} and
 * {@link #equals(Object)} are implemented based on that value as well. That means, of course, that all of those methods
 * only start making sense once {@link #mark} has a meaningful value.
 * </p>
 */
public class HyperEdgeSegment implements Comparable<HyperEdgeSegment> {

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Properties

    /** routing strategy which will ultimately decide how edges will be routed. */
    private final BaseRoutingDirectionStrategy routingStrategy;
    /** ports represented by this hypernode. */
    private final List<LPort> ports = Lists.newArrayList();
    
    /** mark value used for cycle breaking (to be accessed directly). */
    // SUPPRESS CHECKSTYLE NEXT Visibility 
    int mark;
    
    /** the routing slot determines the horizontal distance to the preceding layer. */
    private int routingSlot;
    
    /** start position of this edge segment (in horizontal layouts, this is the topmost y coordinate). */
    private double startPosition = Double.NaN;
    /** end position of this edge segment (in horizontal layouts, this is the bottommost y coordinate). */
    private double endPosition = Double.NaN;
    
    /** sorted list of coordinates where incoming connections enter this segment. */
    private final LinkedList<Double> incomingConnectionCoordinates = Lists.newLinkedList();
    /** sorted list of coordinates where outgoing connections leave this segment. */
    private final LinkedList<Double> outgoingConnectionCoordinates = Lists.newLinkedList();
    
    /** list of outgoing dependencies to other edge segments. */
    private final List<HyperEdgeSegmentDependency> outgoingSegmentDependencies = Lists.newArrayList();
    /** combined weight of all outgoing dependencies. */
    private int outDepWeight;
    /** combined weight of critical outgoing dependencies. */
    private int criticalOutDepWeight;
    /** list of incoming dependencies from other edge segments. */
    private final List<HyperEdgeSegmentDependency> incomingSegmentDependencies = Lists.newArrayList();
    /** combined weight of all incoming dependencies. */
    private int inDepWeight;
    /** combined weight of critical incoming dependencies. */
    private int criticalInDepWeight;

    /** if this segment is the result of a split segment, this is the other segment; otherwise it is {@code null}. */
    private HyperEdgeSegment splitPartner;
    /** the segment that caused this segment to be split, if any (only set on one of the split partners). */
    private HyperEdgeSegment splitBy;

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Initialization

    /**
     * Creates a new instance for the given routing strategy.
     * 
     * @param routingStrategy
     *            the routing strategy. Only required if {@link #addPortPositions(LPort, Map)} will be called to set
     *            this thing up.
     */
    public HyperEdgeSegment(final BaseRoutingDirectionStrategy routingStrategy) {
        this.routingStrategy = routingStrategy;
    }

    /**
     * Adds the positions of the given port and all connected ports.
     *
     * @param port
     *            a port.
     * @param hyperEdgeSegmentMap
     *            map of ports to existing hyperedge segments.
     */
    public void addPortPositions(final LPort port, final Map<LPort, HyperEdgeSegment> hyperEdgeSegmentMap) {
        hyperEdgeSegmentMap.put(port, this);
        ports.add(port);
        double portPos = routingStrategy.getPortPositionOnHyperNode(port);

        // add the new port position to the respective list
        if (port.getSide() == routingStrategy.getSourcePortSide()) {
            insertSorted(incomingConnectionCoordinates, portPos);
        } else {
            insertSorted(outgoingConnectionCoordinates, portPos);
        }
        
        // update start and end coordinates
        recomputeExtent();

        // add connected ports
        for (LPort otherPort : port.getConnectedPorts()) {
            if (!hyperEdgeSegmentMap.containsKey(otherPort)) {
                addPortPositions(otherPort, hyperEdgeSegmentMap);
            }
        }
    }

    private static void insertSorted(final List<Double> list, final double value) {
        ListIterator<Double> listIter = list.listIterator();
        while (listIter.hasNext()) {
            double next = listIter.next().floatValue();
            if (next == value) {
                // an exactly equal value is already present in the list
                return;
            } else if (next > value) {
                listIter.previous();
                break;
            }
        }
        listIter.add(Double.valueOf(value));
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Getters and Setters
    
    /**
     * Returns the ports incident to this segment.
     */
    public List<LPort> getPorts() {
        return ports;
    }

    /**
     * Returns this segment's routing slot.
     */
    public int getRoutingSlot() {
        return routingSlot;
    }

    /**
     * Sets this segment's routing slot.
     */
    public void setRoutingSlot(final int slot) {
        this.routingSlot = slot;
    }

    /**
     * Returns the coordinate where this segment begins.
     */
    public double getStartCoordinate() {
        return startPosition;
    }

    /**
     * Returns the coordinate where this segment ends.
     */
    public double getEndCoordinate() {
        return endPosition;
    }

    /**
     * Returns the (sorted) list of coordinates where incoming connections enter this segment.
     */
    public LinkedList<Double> getIncomingConnectionCoordinates() {
        return incomingConnectionCoordinates;
    }

    /**
     * Returns the (sorted) list of coordinates where outgoing connections leave this segment.
     */
    public LinkedList<Double> getOutgoingConnectionCoordinates() {
        return outgoingConnectionCoordinates;
    }

    /**
     * Return the outgoing dependencies to other hyper edge segments.
     */
    public List<HyperEdgeSegmentDependency> getOutgoingSegmentDependencies() {
        return outgoingSegmentDependencies;
    }

    /**
     * Returns the combined weight of all outgoing dependencies.
     */
    public int getOutWeight() {
        return outDepWeight;
    }

    /**
     * Sets the combined weight of all outgoing dependencies.
     */
    public void setOutWeight(final int outWeight) {
        this.outDepWeight = outWeight;
    }

    /**
     * Returns the combined weight of critical outgoing dependencies.
     */
    public int getCriticalOutWeight() {
        return criticalOutDepWeight;
    }

    /**
     * Sets the combined weight of critical outgoing dependencies.
     */
    public void setCriticalOutWeight(final int outWeight) {
        this.criticalOutDepWeight = outWeight;
    }

    /**
     * Return the incoming dependencies from other hyper edge segments..
     */
    public List<HyperEdgeSegmentDependency> getIncomingSegmentDependencies() {
        return incomingSegmentDependencies;
    }

    /**
     * Returns the weight of incoming dependencies.
     */
    public int getInWeight() {
        return inDepWeight;
    }

    /**
     * Sets the weight of incoming dependencies.
     */
    public void setInWeight(final int inWeight) {
        this.inDepWeight = inWeight;
    }

    /**
     * Returns the combined weight of critical incoming dependencies.
     */
    public int getCriticalInWeight() {
        return criticalInDepWeight;
    }

    /**
     * Sets the combined weight of critical incoming dependencies.
     */
    public void setCriticalInWeight(final int inWeight) {
        this.criticalInDepWeight = inWeight;
    }
    
    /**
     * Returns the split partner, that is, the other segment involved in splitting a segment into two, or {@code null}
     * if this segment was not split. The actual edge routing code will have to know about split partners to route
     * edges accordingly.
     */
    public HyperEdgeSegment getSplitPartner() {
        return splitPartner;
    }
    
    /**
     * Sets the split partner, that is, the other segment involved in splitting a segment into two.
     */
    public void setSplitPartner(final HyperEdgeSegment splitPartner) {
        this.splitPartner = splitPartner;
    }
    
    /**
     * Returns the segment that caused this one to be split, if any. This is only set on one of the split partners.
     */
    public HyperEdgeSegment getSplitBy() {
        return splitBy;
    }
    
    /**
     * Sets the segment that caused this one to be split, if any.
     */
    public void setSplitBy(final HyperEdgeSegment splitBy) {
        this.splitBy = splitBy;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Utilities
    
    /**
     * Convernience method which returns the size of this segment, which is its end coordinate minus its start
     * coordinate.
     */
    public double getLength() {
        return getEndCoordinate() - getStartCoordinate();
    }
    
    /**
     * Checks whether this segment connects two or more ports.
     */
    public boolean representsHyperedge() {
        return getIncomingConnectionCoordinates().size() + getOutgoingConnectionCoordinates().size() > 2;
    }
    
    /**
     * Checks whether this segment was introduced while splitting another segment.
     */
    public boolean isDummy() {
        return splitPartner != null && splitBy == null;
    }
    
    /**
     * Recomputes the start and end coordinate based on incoming and outgoing connection coordinates.
     */
    public void recomputeExtent() {
        startPosition = Double.NaN;
        endPosition = Double.NaN;
        
        recomputeExtent(incomingConnectionCoordinates);
        recomputeExtent(outgoingConnectionCoordinates);
    }
    
    private void recomputeExtent(final LinkedList<Double> positions) {
        // this code assumes that the positions are sorted ascendingly
        if (!positions.isEmpty()) {
            // set new start position
            if (Double.isNaN(startPosition)) {
                startPosition = positions.getFirst();
            } else {
                startPosition = Math.min(startPosition, positions.getFirst());
            }
            
            // set new end position
            if (Double.isNaN(endPosition)) {
                endPosition = positions.getLast();
            } else {
                endPosition = Math.max(endPosition, positions.getLast());
            }
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Splitting
    
    /**
     * Simulates what would happen during a split. The returned pair contains two new {@link HyperEdgeSegment}s: the
     * first simulates what would happen to this instance, and the second simulates the split partner. They are set up
     * with the appropriate incoming and outgoing connection coordinates, except for the linking segment, which this
     * method does not know about. The split partner property of each segment points to the other segment.
     */
    public Pair<HyperEdgeSegment, HyperEdgeSegment> simulateSplit() {
        HyperEdgeSegment newSplit = new HyperEdgeSegment(routingStrategy);
        HyperEdgeSegment newSplitPartner = new HyperEdgeSegment(routingStrategy);
        
        newSplit.incomingConnectionCoordinates.addAll(incomingConnectionCoordinates);
        newSplit.splitBy = splitBy;
        newSplit.splitPartner = newSplitPartner;
        newSplit.recomputeExtent();
        
        newSplitPartner.outgoingConnectionCoordinates.addAll(outgoingConnectionCoordinates);
        newSplitPartner.splitPartner = newSplit;
        newSplitPartner.recomputeExtent();
        
        return Pair.of(newSplit, newSplitPartner);
    }
    
    /**
     * Splits this segment into two and returns the new segment. The segments will be linked at the given position. This
     * segment will retain all incoming connection coordinates, but all of its outgoing connection coordinates will move
     * over to the new segment, to be replaced by the link between the two. This will completely clear all dependencies
     * since they may have become obsolete.
     * 
     * @param splitPosition
     *            position to split the two segments.
     * @return the new segment for convenience, although the two will be linked as split partners.
     */
    public HyperEdgeSegment splitAt(final double splitPosition) {
        splitPartner = new HyperEdgeSegment(routingStrategy);
        splitPartner.setSplitPartner(this);
        
        // Move all target positions over to the new segment
        splitPartner.outgoingConnectionCoordinates.addAll(outgoingConnectionCoordinates);
        this.outgoingConnectionCoordinates.clear();
        
        // Link the two
        this.outgoingConnectionCoordinates.add(splitPosition);
        splitPartner.incomingConnectionCoordinates.add(splitPosition);
        
        // Recompute their outer coordinates
        this.recomputeExtent();
        splitPartner.recomputeExtent();
        
        // Clear dependencies so they can be regenerated later. We could try to be smart about updating them, but that
        // would be more complicated code, so this will do just fine
        while (!incomingSegmentDependencies.isEmpty()) {
            incomingSegmentDependencies.get(0).remove();
        }
        
        while (!outgoingSegmentDependencies.isEmpty()) {
            outgoingSegmentDependencies.get(0).remove();
        }
        
        return splitPartner;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Comparable and Object Overrides

    @Override
    public int compareTo(final HyperEdgeSegment other) {
        return this.mark - other.mark;
    }

    @Override
    public boolean equals(final Object object) {
        if (object instanceof HyperEdgeSegment) {
            HyperEdgeSegment other = (HyperEdgeSegment) object;
            return this.mark == other.mark;
        }
        return false;
    }

    @Override
    public int hashCode() {
        return mark;
    }

    @Override
    public String toString() {
        StringBuilder builder = new StringBuilder("{");
        Iterator<LPort> portIter = ports.iterator();
        while (portIter.hasNext()) {
            LPort port = portIter.next();
            String name = port.getNode().getDesignation();
            if (name == null) {
                name = "n" + port.getNode().getIndex();
            }
            builder.append(name);
            if (portIter.hasNext()) {
                builder.append(',');
            }
        }
        builder.append('}');
        return builder.toString();
    }

}
