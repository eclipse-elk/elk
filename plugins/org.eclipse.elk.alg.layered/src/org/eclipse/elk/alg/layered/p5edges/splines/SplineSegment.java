/*******************************************************************************
 * Copyright (c) 2017, 2018 Kiel University and others.
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
package org.eclipse.elk.alg.layered.p5edges.splines;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.intermediate.FinalSplineBendpointsCalculator;
import org.eclipse.elk.alg.layered.intermediate.LongEdgeSplitter;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.p5edges.splines.SplineEdgeRouter.Dependency;
import org.eclipse.elk.alg.layered.p5edges.splines.SplineEdgeRouter.SideToProcess;
import org.eclipse.elk.core.math.ElkRectangle;
import org.eclipse.elk.core.options.PortSide;
import org.eclipse.elk.core.util.Pair;

/**
 * Represents a segment of a spline that is to be created. A segment represents one or more {@link LEdge}s between
 * adjacent layers. In the case that more than one edge is represented, the segments represent a part of a hyperedge.
 * 
 * Since a segment represents {@link LEdge}(s) between adjacent layers, it may represent {@link LEdge}(s) that were
 * introduced during the splitting of long edges (see {@link LongEdgeSplitter}). The control points for the final spline
 * are created by the {@link FinalSplineBendpointsCalculator}, which runs after long edges have been joined again.
 * Consequently, not all {@link LEdge}s exist anymore. Any information that is required from these edges to create final
 * control points is stored in an accompanying {@link EdgeInformation} class.
 * 
 * <h3>Note on Hyperedges</h3> Splines only support 1:n directed hyperedges. Also all ports on the n-side must lay on
 * the same SideToProcess, looking from the single port. So in a left to right routing, it is not possible to create a
 * hyperedge from an edge going upwards and an edge going downwards.
 * 
 * Furthermore only hyperedges are considered for vertical alignment of edges. All other edges are assumed to go
 * straight to the target port.
 */
public final class SplineSegment implements Comparable<SplineSegment> {
    
    // For convenient usage, members are public
    // SUPPRESS CHECKSTYLE NEXT 70 VisibilityModifier
    /** Spline segments may be part of hyperedges, in which case the {@link FinalSplineBendpointsCalculator} may try to 
     * handle them more than once. We avoid this by flagging the segment as 'handled' as soon as it has been seen. */
    public boolean handled = false;
    
    // the following variables are used during cycle breaking and topological sorting
    /** A sets of ports, determining the ports left of the hyper-edge. */  
    public final Set<LPort> leftPorts = new HashSet<>();
    /** A sets of ports, determining the ports right of the hyper-edge. */  
    public final Set<LPort> rightPorts = new HashSet<>();
    /** Outgoing dependencies are pointing to hyper-edges that must lay right of this hyper edge. */
    public final List<Dependency> outgoing = new ArrayList<>();
    /** Incoming dependencies are pointing to hyper-edges that must lay left of this hyper edge. */
    public final List<Dependency> incoming = new ArrayList<>();
    /** Used to mark nodes in the cycle breaker. */
    public int mark;
    /** Determines how many elements are depending on this. */
    public int inweight;
    /** Determines how many elements this element is depending on. */
    public int outweight;
    /** the rank determines the horizontal distance to the preceding layer. */
    public int rank;
    
    // variables representing the characteristics of the segment 
    /** A set of all LEdges that are combined in this hyper-edge. */
    public final Set<LEdge> edges = new HashSet<>();
    /** If true, the spline segment has no vertical segment, connecting two ports by a horizontal edge path. */
    public final boolean isStraight;
    /** An imaginary bounding box in which this spline segment should reside. */
    public ElkRectangle boundingBox = new ElkRectangle();
    /** Indicates if this segment lies west of the very first layer (e.g. the case for westwards inverted ports). */
    public boolean isWestOfInitialLayer = false;
    /** An x-coordinate offset to be used when the nodeNodeSpacing exceeds the required spacing between a pair of
     *  layers, in which case it looks nicer to move the vertical segments halfway between the layers. */
    public double xDelta;
    
    /** Source port of an edge represented by this segment. */
    public LPort sourcePort;
    /** Target port of an edge represented by this segment. */
    public LPort targetPort;

    /** Flag indicating if this segment is the initial segment of a spline. */
    public boolean initialSegment = false;
    /** Flag indicating if this segment is the last segment of a spline. */
    public boolean lastSegment = false;
    /** If this segment is the {@link SplineSegment#initialSegment initial segment} of a spline, this field holds
     * the source node of the spline. Otherwise it is {@code null}. */
    public LNode sourceNode;
    /** If this segment is the {@link SplineSegment#lastSegment last segment} of a spline, this field holds
     * the target node of the spline. Otherwise it is {@code null}. */
    public LNode targetNode;
    /** Flag indicating that the control points computed for this segment must be interpreted in inverse order. This is
     * necessary for some segments when a {@link org.eclipse.elk.alg.layered.options.WrappingStrategy WrappingStrategy}
     * is used. */
    public boolean inverseOrder = false;

    // hyperedge-related variables
    /** If this segment represents a hyperedge, this is the top y coordinate of the segment that is used during
     *  dependency graph creation for minimize edge crossings of vertical segments. */
    public double hyperEdgeTopYPos;
    /** If this segment represents a hyperedge, this is the bottom y coordinate of the segment that is used during
     *  dependency graph creation for minimize edge crossings of vertical segments. */
    public double hyperEdgeBottomYPos;
    /** If this segment is non-straight, this is the center control point of the spline route to be created. */ 
    public double centerControlPointY;

    
    /**
     * A segment represents an {@link LEdge} between adjacent layers. Thus, it may represent an {@link LEdge}
     * that was introduced during the splitting of long edges (see {@link LongEdgeSplitter}). The final spline 
     * 
     */
    public class EdgeInformation {
        // SUPPRESS CHECKSTYLE NEXT 12 VisibilityModifier|Javadoc
        public double startY;
        public double endY;
        
        public boolean normalSourceNode;
        public boolean normalTargetNode;
        
        public boolean invertedLeft;
        public boolean invertedRight;
    }
    
    /** Map holding for each {@link LEdge} information that was valid during the edge routing phase. */
    public Map<LEdge, EdgeInformation> edgeInformation = new HashMap<>();
    
    /**
     * Constructor for a 1:n hyper-edge.
     * 
     * @param singlePort
     *            The one and only source port.
     * @param edges
     *            All edges that shall be part of this hyper-edge paired with the side (LEFT or RIGHT) 
     *            they are laying on.
     * @param sourceSide
     *            The side of the source.
     */
    public SplineSegment(final LPort singlePort, final Set<Pair<SideToProcess, LEdge>> edges,
            final SideToProcess sourceSide) {
        
        if (sourceSide == SideToProcess.LEFT) {
            leftPorts.add(singlePort);
        } else {
            rightPorts.add(singlePort);
        }

        double yMinPosOfTarget = Double.POSITIVE_INFINITY;
        double yMaxPosOfTarget = Double.NEGATIVE_INFINITY;
        
        for (final Pair<SideToProcess, LEdge> pair : edges) {
            final SideToProcess side = pair.getFirst();
            final LEdge edge = pair.getSecond();
            
            LPort tgtPort = edge.getSource();
            if (tgtPort.equals(singlePort)) {
                tgtPort = edge.getTarget();
            }
            
            if (side == SideToProcess.LEFT) {
                leftPorts.add(tgtPort);
            } else {
                rightPorts.add(tgtPort);
            }
            
            final double yPosOfTarget = anchorY(tgtPort);
            yMinPosOfTarget = Math.min(yMinPosOfTarget, yPosOfTarget);
            yMaxPosOfTarget = Math.max(yMaxPosOfTarget, yPosOfTarget);
        }
            
        final double yPosOfSingleSide = anchorY(singlePort);

        // set the relevant positions 
        setRelevantPositions(yPosOfSingleSide, yMinPosOfTarget, yMaxPosOfTarget);
        
        for (final Pair<SideToProcess, LEdge> pair : edges) {
            addEdge(pair.getSecond());
        }
        isStraight = false;
    }

    /**
     * Constructor for a hyper-edges consisting of a single edge. This MAY be a straight edge.
     * 
     * @param edge
     *            The straight edge that shall be the one and only element of this hyper-edge.
     * @param sourceSide
     *            On witch Layer (LEFT or RIGHT) lays the source port.
     * @param targetSide
     *            On witch Layer (LEFT or RIGHT) lays the target port.
     */
    public SplineSegment(final LEdge edge, final SideToProcess sourceSide, final SideToProcess targetSide) {
        // adding left and right ports
        if (sourceSide == SideToProcess.LEFT) {
            leftPorts.add(edge.getSource());
        } else {
            rightPorts.add(edge.getSource());
        } 
        if (targetSide == SideToProcess.LEFT) {
            leftPorts.add(edge.getTarget());
        } else {
            rightPorts.add(edge.getTarget());
        }
        
        // adding the edges
        addEdge(edge);
        
        // setting relevant positions
        final double sourceY = anchorY(edge.getSource());
        final double targetY = anchorY(edge.getTarget());
        setRelevantPositions(sourceY, targetY, targetY);

        isStraight = SplineEdgeRouter.isStraight(sourceY, targetY);
    }
    
    private void addEdge(final LEdge edge) {
        edges.add(edge);
        
        EdgeInformation ei = new EdgeInformation();
        edgeInformation.put(edge, ei);
        ei.startY = anchorY(edge.getSource());
        ei.endY = anchorY(edge.getTarget());
        ei.normalSourceNode = SplineEdgeRouter.isNormalNode(edge.getSource().getNode());
        ei.normalTargetNode = SplineEdgeRouter.isNormalNode(edge.getTarget().getNode());
        ei.invertedLeft = edge.getSource().getSide() == PortSide.WEST;
        ei.invertedRight = edge.getTarget().getSide() == PortSide.EAST;
    }
    
    private double anchorY(final LPort p) {
        if (PortSide.SIDES_NORTH_SOUTH.contains(p.getSide())) {
            return p.getProperty(InternalProperties.SPLINE_NS_PORT_Y_COORD);
        } else {
            return p.getAbsoluteAnchor().y;
        }
    }
    
    /**
     * @return {@code true} if this spline edge represents a hyperedge, i.e. it represents multiple {@link LEdge}s.
     */
    public boolean isHyperEdge() {
        return edges.size() > 1;
    }

    @Override
    public int compareTo(final SplineSegment other) {
        return this.mark - other.mark;
    }
    
    //////////////////////////////////////////////////
    // Hyper-Edge Constants
    /**
     * Defines the fraction of the outer y position of a hyper-edge for defining the "point of overlap"
     * of two hyper-edges. 1.0 means the point lays on the outer border of the hyper-edge.
     */
    private static final double HYPEREDGE_POS_OUTER_RATE = 0.9;
    /** See RELEVANT_POS_OUTER_RATE! */
    private static final double HYPEREDGE_POS_MID_RATE = 1 - HYPEREDGE_POS_OUTER_RATE;
    /** Avoiding magic number problems. */
    private static final double ONE_HALF = 0.5;
    
    /**
     * Sets the {@link #centerControlPointY} and the top and bottom hyperedge y coordinates to appropriate values.
     * 
     * @param sourceY The y-value of the source point.
     * @param targetYMin The y-value of lowest target point.
     * @param targetYMax The y-value of the highest target point.
     */
    private void setRelevantPositions(final double sourceY, final double targetYMin, final double targetYMax) {
        boundingBox.y = Math.min(sourceY, targetYMin);
        boundingBox.height = Math.max(sourceY, targetYMax) - boundingBox.y;
        
        if (sourceY < targetYMin) {
            // source lays below all target ports
            centerControlPointY = ONE_HALF * (sourceY + targetYMin);
            hyperEdgeTopYPos = HYPEREDGE_POS_MID_RATE * centerControlPointY 
                               + HYPEREDGE_POS_OUTER_RATE * sourceY;
            hyperEdgeBottomYPos = HYPEREDGE_POS_MID_RATE * centerControlPointY 
                                  + HYPEREDGE_POS_OUTER_RATE * targetYMin;
        } else {
            // source lays above all target ports
            centerControlPointY = ONE_HALF * (sourceY + targetYMax);
            hyperEdgeTopYPos = HYPEREDGE_POS_MID_RATE * centerControlPointY
                               + HYPEREDGE_POS_OUTER_RATE * targetYMax;
            hyperEdgeBottomYPos = HYPEREDGE_POS_MID_RATE * centerControlPointY 
                                  + HYPEREDGE_POS_OUTER_RATE * sourceY;
        }
    }
}
