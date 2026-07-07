/*******************************************************************************
 * Copyright (c) 2010, 2020 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.layered.p5edges.orthogonal.direction;

import java.util.Set;

import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.layered.p5edges.orthogonal.HyperEdgeSegment;
import org.eclipse.elk.alg.layered.p5edges.orthogonal.OrthogonalRoutingGenerator;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.math.KVectorChain;
import org.eclipse.elk.core.options.PortSide;

import com.google.common.collect.Sets;

/**
 * A routing direction strategy adapts the {@link OrthogonalRoutingGenerator} to different routing directions. Commonly,
 * edges are routed from a left source layer to a right target layer. However, with northern and southern external
 * ports, this can be different. Routing direction strategies abstract from the actual direction.
 * 
 * <p>Use {@link #forRoutingDirection(RoutingDirection)} to obtain a routing direction strategy for a given
 * {@link RoutingDirection}.</p>
 */
public abstract class BaseRoutingDirectionStrategy {

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Properties

    /** set of already created junction points, to avoid multiple points at the same position. */
    private final Set<KVector> createdJunctionPoints = Sets.newHashSet();
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Creation
    
    /**
     * Returns an implementation suitable for the given routing direction.
     * 
     * @param direction
     *            the direction to return a routing direction strategy for.
     * @return a suitable strategy.
     */
    public static BaseRoutingDirectionStrategy forRoutingDirection(final RoutingDirection direction) {
        switch (direction) {
        case WEST_TO_EAST:
            return new WestToEastRoutingStrategy();
        case NORTH_TO_SOUTH:
            return new NorthToSouthRoutingStrategy();
        case SOUTH_TO_NORTH:
            return new SouthToNorthRoutingStrategy();
        default:
            throw new IllegalArgumentException();
        }
    }
    

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Accessors

    /**
     * Add a junction point to the given edge if necessary. It is necessary to add a junction point if the bend point is
     * not at one of the two end positions of the hyperedge segment.
     *
     * @param edge
     *            an edge
     * @param segment
     *            the corresponding hyperedge segment
     * @param pos
     *            the bend point position
     * @param vertical
     *            {@code true} if the connecting segment is vertical, {@code false} if it is horizontal
     */
    protected void addJunctionPointIfNecessary(final LEdge edge, final HyperEdgeSegment segment, final KVector pos,
            final boolean vertical) {

        double p = vertical ? pos.y : pos.x;
        
        // If we already have this junction point, don't bother
        if (createdJunctionPoints.contains(pos)) {
            return;
        }
        
        // Whether the point lies somewhere inside the edge segment (without boundaries)
        boolean pointInsideEdgeSegment = p > segment.getStartCoordinate() && p < segment.getEndCoordinate();
        
        // Check if the point lies somewhere at the segment's boundary
        boolean pointAtSegmentBoundary = false;
        if (!segment.getIncomingConnectionCoordinates().isEmpty()
                && !segment.getOutgoingConnectionCoordinates().isEmpty()) {
            
            // Is the bend point at the start and joins another edge at the same position?
            pointAtSegmentBoundary |=
                    Math.abs(p - segment.getIncomingConnectionCoordinates().getFirst())
                        < OrthogonalRoutingGenerator.TOLERANCE
                    && Math.abs(p - segment.getOutgoingConnectionCoordinates().getFirst())
                        < OrthogonalRoutingGenerator.TOLERANCE;
            
            // Is the bend point at the end and joins another edge at the same position?
            pointAtSegmentBoundary |= 
                    Math.abs(p - segment.getIncomingConnectionCoordinates().getLast())
                        < OrthogonalRoutingGenerator.TOLERANCE
                    && Math.abs(p - segment.getOutgoingConnectionCoordinates().getLast())
                        < OrthogonalRoutingGenerator.TOLERANCE;
        }
        
        if (pointInsideEdgeSegment || pointAtSegmentBoundary) {
            // create a new junction point for the edge at the bend point's position
            KVectorChain junctionPoints = edge.getProperty(LayeredOptions.JUNCTION_POINTS);
            if (junctionPoints == null) {
                junctionPoints = new KVectorChain();
                edge.setProperty(LayeredOptions.JUNCTION_POINTS, junctionPoints);
            }

            KVector jpoint = new KVector(pos);
            junctionPoints.add(jpoint);
            createdJunctionPoints.add(jpoint);
        }
    }

    /**
     * Returns the set of junction points created so far.
     */
    public Set<KVector> getCreatedJunctionPoints() {
        return createdJunctionPoints;
    }
    
    /**
     * Removes all junction points created so far.
     */
    public void clearCreatedJunctionPoints() {
        createdJunctionPoints.clear();
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // To be implemented

    /**
     * Returns the port's position on a hyper edge axis. In the west-to-east routing case, this would be the port's
     * exact y coordinate.
     *
     * @param port
     *            the port.
     * @return the port's coordinate on the hyper edge axis.
     */
    public abstract double getPortPositionOnHyperNode(LPort port);

    /**
     * Returns the side of ports that should be considered on a source layer. For a west-to-east routing, this would
     * probably be the eastern ports of each western layer.
     *
     * @return the side of ports to be considered in the source layer.
     */
    public abstract PortSide getSourcePortSide();

    /**
     * Returns the side of ports that should be considered on a target layer. For a west-to-east routing, this would
     * probably be the western ports of each eastern layer.
     *
     * @return the side of ports to be considered in the target layer.
     */
    public abstract PortSide getTargetPortSide();

    /**
     * Calculates and assigns bend points for edges incident to the ports belonging to the given hyper edge.
     *
     * @param hyperNode
     *            the hyper edge.
     * @param startPos
     *            the position of the trunk of the first hyper edge between the layers. This position, together with the
     *            current hyper node's rank allows the calculation of the hyper node's trunk's position.
     * @param edgeSpacing
     *            spacing between adjacent edges. This is used to turn an edge segment's routing slot into a proper
     *            position.
     */
    public abstract void calculateBendPoints(HyperEdgeSegment hyperNode, double startPos, double edgeSpacing);

}
