/*******************************************************************************
 * Copyright (c) 2018 Kiel University and others.
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
package org.eclipse.elk.alg.layered.intermediate;

import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LMargin;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.intermediate.compaction.HorizontalGraphCompactor;
import org.eclipse.elk.alg.layered.options.GraphCompactionStrategy;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.layered.options.SplineRoutingMode;
import org.eclipse.elk.alg.layered.p5edges.splines.NubSpline;
import org.eclipse.elk.alg.layered.p5edges.splines.SplineEdgeRouter;
import org.eclipse.elk.alg.layered.p5edges.splines.SplineSegment;
import org.eclipse.elk.alg.layered.p5edges.splines.SplineSegment.EdgeInformation;
import org.eclipse.elk.alg.layered.p5edges.splines.SplinesMath;
import org.eclipse.elk.core.alg.ILayoutProcessor;
import org.eclipse.elk.core.math.ElkMath;
import org.eclipse.elk.core.math.ElkRectangle;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.math.KVectorChain;
import org.eclipse.elk.core.options.PortSide;
import org.eclipse.elk.core.util.IElkProgressMonitor;

/**
 * Intermediate processor that runs some time after the edge routing phase and turns tentative spline routes 
 * into concrete bezier control points that become the bendpoints of {@link LEdge}s. 
 * 
 * For further details on how the control points are assembled, see the {@code calculateControlPoints*(..)} methods.
 * 
 * <dt>Precondition:</dt>
 *    <dd>a layered graph</dd>
 *    <dd>the {@link SplineEdgeRouter} did his job</dd>
 *    <dd>no dummy nodes</dd>
 *  <dt>Postcondition:</dt>
 *    <dd>the bendpoints of edges can be interpreted as Bezier splines</dd>
 *  <dt>Slots:</dt>
 *    <dd>After phase 5.</dd>
 *  <dt>Same-slot dependencies:</dt><dd>After {@link LongEdgeJoiner}</dd>
 *                                  <dd>After {@link LabelDummyRemover}</dd>
 *                                  <dd>After {@link HorizontalGraphCompactor}</dd>
 *                                  <dd>After {@link NorthSouthPortPostprocessor}</dd>
 */
public class FinalSplineBendpointsCalculator implements ILayoutProcessor<LGraph> {

    private double edgeEdgeSpacing;
    private double edgeNodeSpacing;
    private SplineRoutingMode splineRoutingMode;
    private GraphCompactionStrategy compactionStrategy;
    
    /** Avoiding magic number problems. */
    private static final double ONE_HALF = 0.5;
    /** Defines the gap between the source/target anchor of an edge and the control point that is
     *  inserted to straighten the start/end of the edge-spline (for {@link SplineRoutingMode#CONSERVATIVE}). */
    public static final double NODE_TO_STRAIGHTENING_CP_GAP = 5;
    /** Multiplier used to determine the curvature of certain sloppy splines (for {@link SplineRoutingMode#SLOPPY}). */
    private static final double SLOPPY_CENTER_CP_MULTIPLIER = 0.4;
    
    @Override
    public void process(final LGraph graph, final IElkProgressMonitor progressMonitor) {
        this.edgeEdgeSpacing = graph.getProperty(LayeredOptions.SPACING_EDGE_EDGE_BETWEEN_LAYERS);
        this.edgeNodeSpacing = graph.getProperty(LayeredOptions.SPACING_EDGE_NODE_BETWEEN_LAYERS);
        this.splineRoutingMode = graph.getProperty(LayeredOptions.EDGE_ROUTING_SPLINES_MODE);
        this.compactionStrategy = graph.getProperty(LayeredOptions.COMPACTION_POST_COMPACTION_STRATEGY);
        
        // assign indices to nodes to efficiently query neighbors within the same layer
        indexNodesPerLayer(graph);
        
        // collect all edges that represent the first segment of a spline 
        List<LEdge> startEdges = graph.getLayers().stream()
            .flatMap(l -> l.getNodes().stream())
            .flatMap(n -> StreamSupport.stream(n.getOutgoingEdges().spliterator(), false))
            .filter(e -> !e.isSelfLoop())
            .filter(e -> e.hasProperty(InternalProperties.SPLINE_ROUTE_START))
            .collect(Collectors.toList());

        // first determine the NUB control points 
        for (LEdge e : startEdges) {
            List<SplineSegment> spline = e.getProperty(InternalProperties.SPLINE_ROUTE_START);
            spline.forEach(s -> calculateControlPoints(s));
            e.setProperty(InternalProperties.SPLINE_ROUTE_START, null);
        }
        
        // ... then convert them to bezier splines
        for (LEdge e : startEdges) {
            LEdge survivingEdge = e.getProperty(InternalProperties.SPLINE_SURVIVING_EDGE); // may be null
            List<LEdge> edgeChain = e.getProperty(InternalProperties.SPLINE_EDGE_CHAIN);
            calculateBezierBendPoints(edgeChain, survivingEdge);
            // clear property
            e.setProperty(InternalProperties.SPLINE_EDGE_CHAIN, null);
        }
    }
    
    private void indexNodesPerLayer(final LGraph graph) {
        for (Layer l : graph) {
            int index = 0;
            for (LNode n : l) {
                n.id = index++;
            }
        }
    }
    
    /**
     * The method dispatches the calculation of NUB control points for the passed {@code segment} to the
     * {@code calculateControlPoints*(..)} methods based on the following criteria.
     * <ul>
     * <li>the desired {@link SplineRoutingMode}</li>
     * <li>if the segment is straight</li>
     * <li>if the segment represents an in-layer edge due to inverted ports</li>
     * <li>if the segment represents a hyperedge</li>
     * </ul>
     */
    private void calculateControlPoints(final SplineSegment segment) {
        // with hyperedges it can happen that this method is called multiple times for the same segment
        if (segment.handled) {
            return;
        }
        segment.handled = true;

        for (final LEdge edge : segment.edges) {
            if (segment.isStraight && !segment.isHyperEdge()) {
                calculateControlPointsStraight(segment);
                continue;
            }
            
            // Remember that the edge itself is not necessarily valid at this point 
            //  (it may have been remove by the long edge joiner, for instance)
            EdgeInformation ei = segment.edgeInformation.get(edge);
            // inverted ports are handled in the same way for 
            if (ei.invertedLeft || ei.invertedRight) {
                calculateControlPointsInvertedEdge(edge, segment);
                continue;
            }
            
            // to compute sloppy control points at least one of the 
            //  two nodes connected by the segment must be a 'normal' node,
            //  and in case horizontal compaction is active the predicate implemented 
            //  by the 'segmentAllowsSloppyRouting(..)' method must be true
            // for hyperedges sloppy routing is not possible
            boolean sloppy = splineRoutingMode == SplineRoutingMode.SLOPPY
                    && (ei.normalSourceNode || ei.normalTargetNode)
                    && segmentAllowsSloppyRouting(segment)
                    && !segment.isHyperEdge();
            if (sloppy) {
                calculateControlPointsSloppy(edge, segment);
            } else {
                calculateControlPointsConservative(edge, segment);
            }
        }
        
        if (segment.inverseOrder) {
            segment.edges.forEach(e -> Collections.reverse(e.getBendPoints()));
        }
    }
    
    /**
     * Adds a single control point to a straight segment, halfway between the source and target layer.
     */
    private void calculateControlPointsStraight(final SplineSegment segment) {
        double xStartPos = segment.boundingBox.x;
        double xEndPos = segment.boundingBox.x + segment.boundingBox.width;
        KVector halfway = new KVector(xStartPos + (xEndPos - xStartPos) / 2, segment.centerControlPointY); 
        segment.edges.iterator().next().getBendPoints().add(halfway);
    }
    
    /**
     * If {@code edge} is an inverted edge on the western side, control points ({@code +}) are computed as follows: 
     * <pre>
     *        ___
     *       |   |
     *   +--+|   |
     *   |   |   |
     *   |   |   |
     *   +--+|   |
     *       |___|
     * </pre>
     * An inverted edge on the right side works analog.
     */
    private void calculateControlPointsInvertedEdge(final LEdge edge, final SplineSegment containingSegment) {
        final double startXPos = containingSegment.boundingBox.x;
        final double endXPos = containingSegment.boundingBox.x + containingSegment.boundingBox.width;
        
        final EdgeInformation ei = containingSegment.edgeInformation.get(edge);
        final double ySourceAnchor = ei.startY;
        final double yTargetAnchor = ei.endY;
        
        // compute the desired control points
        final KVector sourceStraightCP;
        final KVector targetStraightCP;
        if (ei.invertedLeft) {
            sourceStraightCP = new KVector(endXPos, ySourceAnchor);
        } else {
            sourceStraightCP = new KVector(startXPos, ySourceAnchor);
        }
        if (ei.invertedRight) {
            targetStraightCP = new KVector(startXPos, yTargetAnchor);
        } else {
            targetStraightCP = new KVector(endXPos, yTargetAnchor);
        }
    
        // the center position is the same for all edges but depends on sloppiness of the routing
        double centerXPos = startXPos;
        if (!containingSegment.isWestOfInitialLayer) {
            centerXPos += edgeNodeSpacing;
        }
        centerXPos += containingSegment.xDelta + (containingSegment.rank + 0) * edgeEdgeSpacing;
        
        final KVector sourceVerticalCP = new KVector(centerXPos, ySourceAnchor);
        final KVector targetVerticalCP = new KVector(centerXPos, yTargetAnchor);
        
        // add control points to the edge's bendpoints
        edge.getBendPoints().addAll(sourceStraightCP, sourceVerticalCP);
        boolean isHyperedge = containingSegment.edges.size() > 1; 
        if (isHyperedge) {
            // add an additional center control point to assert that the hyperedge segments  
            // share a part of their route
            final KVector center = new KVector(centerXPos, containingSegment.centerControlPointY);                    
            edge.getBendPoints().add(center);
        }
        edge.getBendPoints().addAll(targetVerticalCP, targetStraightCP);
    }
    
    /**
     * For a regular upward pointing edge and the conservative routing style, control points ({@code +}) are computed 
     * as follows:
     * <pre>
     *               ___
     *              |   |
     *          +--+|   |
     *   ___    |   |___|
     *  |   |   |
     *  |   |+--+
     *  |___|
     *  </pre>
     * An downward pointing edge is handled analog. 
     */
    private void calculateControlPointsConservative(final LEdge edge, final SplineSegment containingSegment) {
        final double startXPos = containingSegment.boundingBox.x;
        final double endXPos = containingSegment.boundingBox.x + containingSegment.boundingBox.width;
               
        final EdgeInformation ei = containingSegment.edgeInformation.get(edge);
        final double ySourceAnchor = ei.startY;
        final double yTargetAnchor = ei.endY;
        
        // Calculate bend points to draw inner layer segments straight
        // to prevent intersections with big nodes
        final KVector sourceStraightCP = new KVector(startXPos, ySourceAnchor);
        final KVector targetStraightCP = new KVector(endXPos, yTargetAnchor);

        double centerXPos = startXPos;
        if (!containingSegment.isWestOfInitialLayer) {
            centerXPos += edgeNodeSpacing;
        }
        centerXPos += containingSegment.xDelta + (containingSegment.rank + 0) * edgeEdgeSpacing;
        final KVector sourceVerticalCP = new KVector(centerXPos, ySourceAnchor);
        final KVector targetVerticalCP = new KVector(centerXPos, yTargetAnchor);

        // Traditional four control points (plus an extra center control point for hyperedges) 
        edge.getBendPoints().addAll(sourceStraightCP, sourceVerticalCP);
        boolean isHyperedge = containingSegment.edges.size() > 1; 
        if (isHyperedge) {
            // add an additional center control point to assert that the hyperedge segments share a part of their route
            final KVector center = new KVector(centerXPos, containingSegment.centerControlPointY);                    
            edge.getBendPoints().add(center);
        }
        edge.getBendPoints().addAll(targetVerticalCP, targetStraightCP);
    }
    
    /**
     * As opposed to the conservative strategy to insert control points (see
     * {@link #calculateControlPointsConservative(LEdge, SplineSegment)}), the sloppy version tries to use less control
     * points to make the splines curvier. The method handles only the case in which one of the two involved nodes is a
     * {@link SplineEdgeRouter#isNormalNode(LNode) normal} node.
     * 
     * Let {@code e = (u,v)} be the edge for which control points are to be determined. The conservative strategy adds
     * control points {@code a,b,c,d} as depicted in the figure below. For a normal node {@code u}, it is checked if
     * {@code a,b} can be omitted. This is the case if the path of the resulting spline doesn't intersect a node that
     * lies above or below {@code u}. The same is checked for a regular node {@code v} and the control points
     * {@code c,d}. In case all control points could be omitted, a single control point is added halfway in between the
     * two nodes, where the y-coordinate is computed as implemented in
     * {@link #computeSloppyCenterY(LEdge, double, double)}. Note that if no control points can be omitted, 
     * the same control points are used as with the conservative strategy.
     * 
     * 
     * <pre>
     *                v
     *               ___
     *              |   |
     *    u     c--d|   |
     *   ___    |   |___|
     *  |   |   |
     *  |   |a--b
     *  |___|
     * </pre>
     * 
     * To illustrate the three possible scenarios in which control points are omitted. 
     * <pre>
     *   shortcut src              shortcut tgt            shortcut both
     * </pre>
     * 
     * <pre>
     *               ___                      ___                     ___
     *              |   |                    |   |                   |   |
     *          +--+|   |                .---|   |               .---|   |
     *   ___    |   |___|         ___    |   |___|        ___    +   |___|
     *  |   |   |                |   |   |               |   |   |
     *  |   |---`                |   |+--+               |   |---`
     *  |___|                    |___|                   |___|
     * </pre>
     * 
     */
    private void calculateControlPointsSloppy(final LEdge edge, final SplineSegment containingSegment) {
        EdgeInformation ei = containingSegment.edgeInformation.get(edge);
        assert ei.normalSourceNode || ei.normalTargetNode;
        
        double startXPos = containingSegment.boundingBox.x;
        double endXPos = containingSegment.boundingBox.x + containingSegment.boundingBox.width;
        
        final double ySourceAnchor = ei.startY;
        final double yTargetAnchor = ei.endY;
        final boolean edgePointsDownwards = ySourceAnchor < yTargetAnchor;
    
        // pre-compute a number of coordinates that we might use as control points
        final KVector sourceStraightCP = new KVector(startXPos, ySourceAnchor);
        final KVector targetStraightCP = new KVector(endXPos, yTargetAnchor);
        final double centerXPos = (startXPos + endXPos) / 2;
        final KVector sourceVerticalCP = new KVector(centerXPos, ySourceAnchor);
        final KVector targetVerticalCP = new KVector(centerXPos, yTargetAnchor);
        
        // evaluate if a rather direct curve is possible
        double centerYPos = computeSloppyCenterY(edge, ySourceAnchor, yTargetAnchor);
        KVector v1 = containingSegment.sourcePort.getAbsoluteAnchor();
        KVector v2 = new KVector(centerXPos, centerYPos);
        KVector v3 = containingSegment.targetPort.getAbsoluteAnchor();
        // approx will be an array holding two values, the zeroth of which represents the 
        // center of the curve
        KVector[] approx = ElkMath.approximateBezierSegment(2, v1, v2, v3);
        
        boolean shortCutSource = false;
        final LNode src = containingSegment.sourcePort.getNode();
        // when graph wrapping is activated, it can happen that the 'src' node doesn't exist anymore (the same goes for
        // the node's layer)
        if (src != null && src.getLayer() != null && ei.normalSourceNode) {
            // possible intersections must only be checked if there are nodes 
            // with which the edge could potentially intersect
            final boolean needToCheckSrc = 
                        edgePointsDownwards && src.id < src.getLayer().getNodes().size() - 1
                    || !edgePointsDownwards && src.id > 0;
        
            if (!needToCheckSrc) {
                shortCutSource = true;
            } else {
                // check within src's layer
                if (needToCheckSrc) {
                    int neighborIndex = src.id;
                    if (edgePointsDownwards) {
                        neighborIndex++;
                    } else { 
                        neighborIndex--; 
                    }
                    LNode neighbor = src.getLayer().getNodes().get(neighborIndex);
                    ElkRectangle box = nodeToBoundingBox(neighbor);
                    shortCutSource = !(ElkMath.intersects(box, v1, approx[0]) || ElkMath.contains(box, v1, approx[0]));
                }            
                
            }
        }
        
        boolean shortCutTarget = false;
        final LNode tgt = containingSegment.targetPort.getNode();
        // see comment above
        if (tgt != null && tgt.getLayer() != null && ei.normalTargetNode) {
            final boolean needToCheckTgt = 
                        edgePointsDownwards && tgt.id > 0
                    || !edgePointsDownwards && tgt.id < tgt.getLayer().getNodes().size() - 1;
     
            // tgt's layer
            if (!needToCheckTgt) {
                shortCutTarget = true;
            } else {
                int neighborIndex = tgt.id;
                if (edgePointsDownwards) {
                    neighborIndex--;
                } else {
                    neighborIndex++;
                }
                LNode neighbor = tgt.getLayer().getNodes().get(neighborIndex);
                ElkRectangle box = nodeToBoundingBox(neighbor);
                shortCutTarget = !(ElkMath.intersects(box, approx[0], v3) || ElkMath.contains(box, approx[0], v3));
            }
        }
        
        // now add the control points 
        if (shortCutSource && shortCutTarget) {
            edge.getBendPoints().add(v2);
        }
        if (!shortCutSource) {
            edge.getBendPoints().addAll(sourceStraightCP, sourceVerticalCP);
        }
        if (!shortCutTarget) {
            edge.getBendPoints().addAll(targetVerticalCP, targetStraightCP);
        }
    }
    
    private ElkRectangle nodeToBoundingBox(final LNode node) {
        KVector pos = node.getPosition();
        KVector size = node.getSize();
        LMargin m = node.getMargin();
        return new ElkRectangle(pos.x - m.left, pos.y - m.top, size.x + m.getHorizontal(), size.y + m.getVertical());
    }
    
    private double computeSloppyCenterY(final LEdge edge, final double ySourceAnchor, final double yTargetAnchor) {
        int indegree = 0;
        int outdegree = 0;
        // Count all the incoming and outgoing edges
        //  remember that the edge's source and target may have been set to null by now 
        //  (in case they were connected to a dummy node beforehand), in which case the degree is assumed to be one
        if (edge.getSource() != null) {
            for (LPort port : edge.getTarget().getNode().getPorts()) {
                indegree += port.getIncomingEdges().size();
            }
        } else {
            indegree = 1;
        }
        if (edge.getTarget() != null) {
            for (LPort port : edge.getSource().getNode().getPorts()) {
                outdegree += port.getOutgoingEdges().size();
            }
        } else {
            outdegree = 1;
        }
        final int degreeDiff = (int) Math.signum(outdegree - indegree);
        final double centerYPos = ((yTargetAnchor + ySourceAnchor) / 2)
                + (yTargetAnchor - ySourceAnchor) * (SLOPPY_CENTER_CP_MULTIPLIER * degreeDiff);
        
        return centerYPos;
    }
    
    /**
     * Collects all NUB control points computed for a spline segment and converts them to bezier control points.
     * Note that additional NUB control points <em>may</em> be added before the conversion is performed.   
     */
    private void calculateBezierBendPoints(final List<LEdge> edgeChain, final LEdge survivingEdge) {
        if (edgeChain.isEmpty()) {
            return;
        }
        
        // in this chain we will put all NURBS control points.
        final KVectorChain allCP = new KVectorChain();
        // add the computed bendpoints to the specified edge (default to the first edge in the edge chain)
        final LEdge edge = survivingEdge != null ? survivingEdge : edgeChain.get(0);
        // Process the source end of the edge-chain. 
        final LPort sourcePort = edge.getSource();
        
        // edge must be the first edge of a chain of edges
        if (!SplineEdgeRouter.isQualifiedAsStartingNode(sourcePort.getNode())) {
            throw new IllegalArgumentException(
                    "The target node of the edge must be a normal node " + "or a northSouthPort.");
        }
        
        // add the source as the very first control point.
        allCP.addLast(sourcePort.getAbsoluteAnchor());
        
        // add an additional control point if the source port is a north or south port
        if (PortSide.SIDES_NORTH_SOUTH.contains(sourcePort.getSide())) {
            double y = sourcePort.getProperty(InternalProperties.SPLINE_NS_PORT_Y_COORD);
            KVector northSouthCP = new KVector(sourcePort.getAbsoluteAnchor().x, y);
            allCP.addLast(northSouthCP);
        }
        
        // copy the calculated control points for all spline segments,
        //  possibly adding additional control points halfway between computed ones
        KVector lastCP = null;
        boolean addMidPoint = false;
        Iterator<LEdge> edgeIterator = edgeChain.iterator();
        while (edgeIterator.hasNext()) {
            LEdge currentEdge = edgeIterator.next(); 
            // read the stored bend-points for vertical segments, calculated by calculateNUBSBendPoint.
            final KVectorChain currentBendPoints = currentEdge.getBendPoints();

            if (!currentBendPoints.isEmpty()) {
                // add a CP in the middle of the straight segment between two vertical segments to
                // get a more straight horizontal segment
                if (addMidPoint) {
                    KVector halfway = lastCP.add(currentBendPoints.getFirst()).scale(ONE_HALF);
                    allCP.addLast(halfway);
                    addMidPoint = false;
                } else {
                    addMidPoint = true;
                }
                lastCP = currentBendPoints.getLast().clone();
                allCP.addAll(currentBendPoints);
                currentBendPoints.clear();
            }
        }

        // finalize the spline
        LPort targetPort = edge.getTarget();
        
        // again, add an additional control point if the target port is a north or sout port
        if (PortSide.SIDES_NORTH_SOUTH.contains(targetPort.getSide())) {
            double y = targetPort.getProperty(InternalProperties.SPLINE_NS_PORT_Y_COORD);
            KVector northSouthCP = new KVector(targetPort.getAbsoluteAnchor().x, y);
            allCP.addLast(northSouthCP);
        }
        
        // finish with the target as last control point.
        allCP.addLast(targetPort.getAbsoluteAnchor());

        // insert straightening control points (if desired)
        if (splineRoutingMode == SplineRoutingMode.CONSERVATIVE) {
            // Add a control point for a straight segment at the very start and at the very end of a spline to prevent
            //  the edge from colliding with self-loops or the like inside the margin of the node. 
            //  This also ensures the correct initial direction of the edge
            insertStraighteningControlPoints(allCP, sourcePort, targetPort);
        }

        // convert list of NUB control points to bezier control points
        final NubSpline nubSpline = new NubSpline(true, SplineEdgeRouter.SPLINE_DIMENSION, allCP);
        // ... and set them as bendpoints of the edge
        edge.getBendPoints().addAll(nubSpline.getBezierCP());
    }
    
    private void insertStraighteningControlPoints(final KVectorChain allCPs, final LPort srcPort, final LPort tgtPort) {
        // KVectorChain extends LinkedList, thus the following operations are fast
        // Further note that the vector computations are necessary to address adding a straightening gap 
        //  into either of the for cardinal directions (based on the port side)
        
        // beginning 
        KVector first = allCPs.getFirst();
        KVector second = allCPs.get(1);
        
        KVector v = new KVector(SplinesMath.portSideToDirection(srcPort.getSide()));
        v.scale(NODE_TO_STRAIGHTENING_CP_GAP);
        KVector v2 = second.clone().sub(first);
        KVector straightenBeginning = new KVector(absMin(v.x, v2.x), absMin(v.y, v2.y));  
        straightenBeginning.add(first);

        allCPs.add(1, straightenBeginning);
        
        // ending
        KVector last = allCPs.getLast();
        KVector secondLast = allCPs.get(allCPs.size() - 2);

        v = new KVector(SplinesMath.portSideToDirection(tgtPort.getSide()));
        v.scale(NODE_TO_STRAIGHTENING_CP_GAP);
        v2 = secondLast.clone().sub(last);
        KVector straightenEnding = new KVector(absMin(v.x,  v2.x), absMin(v.y, v2.y));
        straightenEnding.add(last);

        allCPs.add(allCPs.size() - 1, straightenEnding);
    }
    
    private double absMin(final double d1, final double d2) {
        if (Math.abs(d1) < Math.abs(d2)) {
            return d1;
        } else {
            return d2;
        }
    }
    
    private boolean segmentAllowsSloppyRouting(final SplineSegment segment) {
        // only check this if one dimensional compaction is applied
        if (compactionStrategy == GraphCompactionStrategy.NONE) {
            return true;
        }
        
        double startXPos = segment.boundingBox.x;
        double endXPos = segment.boundingBox.x + segment.boundingBox.width;

        if (segment.initialSegment) {
            LNode n = segment.sourceNode;
            double t = segmentNodeDistanceThreshold(n);
            double nodeSegmentDistance = startXPos - (n.getPosition().x + n.getSize().x); 
            if (nodeSegmentDistance > t) {
                return false;
            }
        }
        if (segment.lastSegment) {
            LNode n = segment.targetNode;
            double t = segmentNodeDistanceThreshold(n);
            double nodeSegmentDistance = n.getPosition().x - endXPos;
            if (nodeSegmentDistance > t) {
                return false;
            }
        }

        return true;
    }
    
    /**
     * @return a magically determined threshold value.  
     */
    private double segmentNodeDistanceThreshold(final LNode n) {
        return n.getLayer().getSize().x - n.getSize().x / 2;
    }
}
