/*******************************************************************************
 * Copyright (c) 2010, 2020 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.layered.p5edges.splines;

import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.ListIterator;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import org.eclipse.elk.alg.layered.LayeredPhases;
import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LGraphUtil;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LNode.NodeType;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.intermediate.FinalSplineBendpointsCalculator;
import org.eclipse.elk.alg.layered.intermediate.IntermediateProcessorStrategy;
import org.eclipse.elk.alg.layered.options.GraphProperties;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.layered.options.SplineRoutingMode;
import org.eclipse.elk.alg.layered.p5edges.PolylineEdgeRouter;
import org.eclipse.elk.core.alg.ILayoutPhase;
import org.eclipse.elk.core.alg.LayoutProcessorConfiguration;
import org.eclipse.elk.core.options.PortSide;
import org.eclipse.elk.core.util.IElkProgressMonitor;
import org.eclipse.elk.core.util.Pair;

import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;

/**
 * Implements a way of routing the edges with splines. Uses the dummy nodes as reference points for
 * a spline calculation, but the dummy nodes do not necessarily lay on the edge. They are only approximated.
 * 
 * <h3>Implementation Details.</h3>
 * A convenience class {@link SplineSegment} is used to represent spline segments between adjacent layers.
 * The segments are combined to a contiguous spline later by the {@link FinalSplineBendpointsCalculator}. As such, 
 * <b>no final</b> bend points are calculated by this layout phase implementation. The x-coordinates of the nodes are 
 * fixed, however, usually to <code> max{nodeNodeSpacing, 2 * edgeNodeSpacing + (slots-1) * edgeEdgeSpacing}</code>.
 * At this, a <code>slot</code> is a vertical strip between a pair of layers in which a non-straight edge will be 
 * routed. The more edges in-between a pair of layers, the more slots are required. 
 * If the splines are to be routed {@link SplineRoutingMode#SLOPPY}, the spacing may be increased further 
 * as computed by the {@link #computeSloppySpacing(Layer, double, double, double)} method.
 * 
 * <dl>
 *   <dt>Precondition:</dt>
 *     <dd>the graph has a proper layering with assigned node and port positions</dd>
 *     <dd>the size of each layer is correctly set</dd>
 *   <dt>Postcondition:</dt>
 *     <dd>each node is assigned a horizontal coordinate</dd>
 *     <dd>the width of the whole graph is set</dd>
 *     <dd>the bend points of each edge are <b>not yet set</b></dd>
 * </dl>
 * 
 * @see FinalSplineBendpointsCalculator
 */
public final class SplineEdgeRouter implements ILayoutPhase<LayeredPhases, LGraph> {
    
    // /////////////////////////////////////////////////////////////////////////////
    // Constants
    /** An edge is drawn as a straight line if the y-difference of source/target is lower than this. */
    private static final double MAX_VERTICAL_DIFF_FOR_STRAIGHT = 0.2;
    /** Default dimension of an edge-spline. */
    public static final int SPLINE_DIMENSION = 3;

    //////////////////////////////////////////////////
    // Intermediate processing configurations
    /** baseline processor dependencies, always required. */
    private static final LayoutProcessorConfiguration<LayeredPhases, LGraph> BASELINE_PROCESSING_ADDITIONS =
            LayoutProcessorConfiguration.<LayeredPhases, LGraph>create()
                .addAfter(LayeredPhases.P5_EDGE_ROUTING,
                    IntermediateProcessorStrategy.FINAL_SPLINE_BENDPOINTS_CALCULATOR)
                .addBefore(LayeredPhases.P3_NODE_ORDERING, IntermediateProcessorStrategy.INVERTED_PORT_PROCESSOR);
    
    /** additional processor dependencies for graphs with self-loops. */
    private static final LayoutProcessorConfiguration<LayeredPhases, LGraph> SELF_LOOP_PROCESSING_ADDITIONS =
            LayoutProcessorConfiguration.<LayeredPhases, LGraph>create()
                .addBefore(LayeredPhases.P1_CYCLE_BREAKING, IntermediateProcessorStrategy.SELF_LOOP_PREPROCESSOR)
                .addAfter(LayeredPhases.P5_EDGE_ROUTING, IntermediateProcessorStrategy.SELF_LOOP_POSTPROCESSOR)
                .before(LayeredPhases.P4_NODE_PLACEMENT)
                    .add(IntermediateProcessorStrategy.SELF_LOOP_PORT_RESTORER)
                    .add(IntermediateProcessorStrategy.SELF_LOOP_ROUTER);

    /** additional processor dependencies for graphs with center edge labels. */
    private static final LayoutProcessorConfiguration<LayeredPhases, LGraph> CENTER_EDGE_LABEL_PROCESSING_ADDITIONS =
            LayoutProcessorConfiguration.<LayeredPhases, LGraph>create()
                .addBefore(LayeredPhases.P2_LAYERING, IntermediateProcessorStrategy.LABEL_DUMMY_INSERTER)
                .addBefore(LayeredPhases.P4_NODE_PLACEMENT, IntermediateProcessorStrategy.LABEL_DUMMY_SWITCHER)
                .addBefore(LayeredPhases.P4_NODE_PLACEMENT, IntermediateProcessorStrategy.LABEL_SIDE_SELECTOR)
                .addAfter(LayeredPhases.P5_EDGE_ROUTING, IntermediateProcessorStrategy.LABEL_DUMMY_REMOVER);

    /** additional processor dependencies for graphs with northern / southern non-free ports. */
    private static final LayoutProcessorConfiguration<LayeredPhases, LGraph> NORTH_SOUTH_PORT_PROCESSING_ADDITIONS =
            LayoutProcessorConfiguration.<LayeredPhases, LGraph>create()
                .addBefore(LayeredPhases.P3_NODE_ORDERING, IntermediateProcessorStrategy.NORTH_SOUTH_PORT_PREPROCESSOR)
                .addBefore(LayeredPhases.P5_EDGE_ROUTING, IntermediateProcessorStrategy.NORTH_SOUTH_PORT_POSTPROCESSOR);

    /** additional processor dependencies for graphs with head or tail edge labels. */
    private static final LayoutProcessorConfiguration<LayeredPhases, LGraph> END_EDGE_LABEL_PROCESSING_ADDITIONS =
            LayoutProcessorConfiguration.<LayeredPhases, LGraph>create()
                .addBefore(LayeredPhases.P4_NODE_PLACEMENT, IntermediateProcessorStrategy.LABEL_SIDE_SELECTOR)
                .addBefore(LayeredPhases.P4_NODE_PLACEMENT, IntermediateProcessorStrategy.END_LABEL_PREPROCESSOR)
                .addAfter(LayeredPhases.P5_EDGE_ROUTING, IntermediateProcessorStrategy.END_LABEL_POSTPROCESSOR);
    
    @Override
    public LayoutProcessorConfiguration<LayeredPhases, LGraph> getLayoutProcessorConfiguration(final LGraph graph) {
        // Basic configuration
        final LayoutProcessorConfiguration<LayeredPhases, LGraph> configuration =
                LayoutProcessorConfiguration.<LayeredPhases, LGraph>create();
        
        configuration.addAll(BASELINE_PROCESSING_ADDITIONS);

        final Set<GraphProperties> graphProperties =
                graph.getProperty(InternalProperties.GRAPH_PROPERTIES);

        if (graphProperties.contains(GraphProperties.SELF_LOOPS)) {
            configuration.addAll(SELF_LOOP_PROCESSING_ADDITIONS);
        }

        if (graphProperties.contains(GraphProperties.CENTER_LABELS)) {
            configuration.addAll(CENTER_EDGE_LABEL_PROCESSING_ADDITIONS);
        }

        if (graphProperties.contains(GraphProperties.NORTH_SOUTH_PORTS)) {
            configuration.addAll(NORTH_SOUTH_PORT_PROCESSING_ADDITIONS);
        }

        if (graphProperties.contains(GraphProperties.END_LABELS)) {
            configuration.addAll(END_EDGE_LABEL_PROCESSING_ADDITIONS);
        }
        return configuration;
    }
    
    //////////////////////////////////////////////////
    
    // some variables valid during one iteration (a pair of layers)
    /** current edges of the layers in a current iteration. */
    private final List<LEdge> edgesRemainingLayer = Lists.newArrayList();
    /** current spline segments of the layers in a current iteration. */
    private final List<SplineSegment> splineSegmentsLayer = Lists.newArrayList();
    /** current ports on the left layer involved in a current iteration. */
    private final Set<LPort> leftPortsLayer = Sets.newLinkedHashSet();
    /** current ports on the right layer involved in a current iteration. */
    private final Set<LPort> rightPortsLayer = Sets.newLinkedHashSet();
    /** current self loops of the layers in a current iteration. */
    private final Set<LEdge> selfLoopsLayer = Sets.newLinkedHashSet();
    
    // variables for the whole edge routing process
    private LGraph lGraph;
    /** a collection of all edges that have a normal node as their source. */
    private final List<LEdge> startEdges = Lists.newArrayList();
    /** all created spline segments. */
    private final List<SplineSegment> allSplineSegments = Lists.newArrayList();
    /** Maps {@link LEdge} to their representing segments. */
    private final Map<LEdge, SplineSegment> edgeToSegmentMap = Maps.newHashMap();
    /** A mapping pointing from an edge to it's succeeding edge, together with their connected 
       friends, they form a long-edge. */
    private final Map<LEdge, LEdge> successingEdge = Maps.newHashMap();
    
    //////////////////////////////////////////////////

    @Override
    public void process(final LGraph layeredGraph, final IElkProgressMonitor monitor) {
        monitor.begin("Spline edge routing", 1);
        
        if (layeredGraph.getLayers().isEmpty()) {
            layeredGraph.getSize().x = 0;
            monitor.done();
            return;
        }
        
        // Retrieve some generic values
        final double nodeNodeSpacing = layeredGraph.getProperty(LayeredOptions.SPACING_NODE_NODE_BETWEEN_LAYERS);
        final double edgeNodeSpacing = layeredGraph.getProperty(LayeredOptions.SPACING_EDGE_NODE_BETWEEN_LAYERS);
        final double edgeEdgeSpacing = layeredGraph.getProperty(LayeredOptions.SPACING_EDGE_EDGE_BETWEEN_LAYERS);
        
        // Find out if splines should be routed thoroughly or sloppy
        final SplineRoutingMode mode = layeredGraph.getProperty(LayeredOptions.EDGE_ROUTING_SPLINES_MODE);
        final boolean sloppyRouting = mode == SplineRoutingMode.SLOPPY;
        final double sloppyLayerSpacingFactor =
                layeredGraph.getProperty(LayeredOptions.EDGE_ROUTING_SPLINES_SLOPPY_LAYER_SPACING_FACTOR);
        
        lGraph = layeredGraph;
        startEdges.clear();
        allSplineSegments.clear();
        successingEdge.clear();

        // check if the first and/or last layer are populated with external port dummies
        final Layer firstLayer = layeredGraph.getLayers().get(0);
        final boolean isLeftLayerExternal =
                Iterables.all(firstLayer.getNodes(), PolylineEdgeRouter.PRED_EXTERNAL_WEST_OR_EAST_PORT);
        final Layer lastLayer = layeredGraph.getLayers().get(layeredGraph.getLayers().size() - 1);
        final boolean isRightLayerExternal =
                Iterables.all(lastLayer.getNodes(), PolylineEdgeRouter.PRED_EXTERNAL_WEST_OR_EAST_PORT);

        final Iterator<Layer> layerIterator = layeredGraph.iterator();
        Layer leftLayer = null;
        Layer rightLayer;
        
        // initial x position 
        double xpos = 0.0;
        do {
            rightLayer = layerIterator.hasNext() ? layerIterator.next() : null;
            
            // fresh start for this pair of layers
            clearThenFillMappings(leftLayer, rightLayer);

            // creation of the SplineSegments
            createSegmentsAndComputeRanking();

            // count the number of required slots for vertical segments
            //  (edges to be drawn straight are assigned a rank but must be omitted here)
            final int slotCount = splineSegmentsLayer.stream()
                    .filter(e -> !e.isStraight)
                    .mapToInt(e -> e.rank + 1).max().orElse(0);
            
            // the code below ensures that at least nodeNodeSpacing is preserved between a pair of layers
            //  if this spacing is larger than what would be required to route the vertical segments in-between
            //  a pair of layers, it looks nicer to move the vertical segments halfway between the layers.
            //  The xSegmentDelta variable holds the required offset 
            double xSegmentDelta = 0;
            double rightLayerPosition = xpos;
            boolean isSpecialLeftLayer = leftLayer == null || (isLeftLayerExternal && leftLayer == firstLayer);
            boolean isSpecialRightLayer = rightLayer == null || (isRightLayerExternal && rightLayer == lastLayer);
            
            // compute horizontal positions just as for the OrthogonalEdgeRouter
            if (slotCount > 0) {
                // the space between each pair of edge segments, and between nodes and edges
                double increment = 0;
                if (leftLayer != null) {
                    increment += edgeNodeSpacing;
                }
                increment += (slotCount - 1) * edgeEdgeSpacing;
                if (rightLayer != null) {
                    increment += edgeNodeSpacing;
                }
                
                // sloppy routing may want to reserve more space in-between a pair of layers
                if (sloppyRouting && rightLayer != null) {
                    increment = Math.max(increment, 
                          computeSloppySpacing(rightLayer, edgeEdgeSpacing, nodeNodeSpacing, sloppyLayerSpacingFactor));
                }
                
                // if we are between two layers, make sure their minimal spacing is preserved
                if (increment < nodeNodeSpacing && !isSpecialLeftLayer && !isSpecialRightLayer) {
                    xSegmentDelta = (nodeNodeSpacing - increment) / 2d;
                    increment = nodeNodeSpacing;
                }
                rightLayerPosition += increment;
            } else if (!isSpecialLeftLayer && !isSpecialRightLayer) {
                // If all edges are straight, use the usual spacing 
                rightLayerPosition += nodeNodeSpacing;
            }

            // place right layer's nodes
            if (rightLayer != null) {
                LGraphUtil.placeNodesHorizontally(rightLayer, rightLayerPosition);
            }
            
            // Assign tentative start and end points to the spline segments
            //  they may be modified before final spline coordinates 
            //  are determined by the FinalSplineBendpointsCaluclator
            for (final SplineSegment segment : splineSegmentsLayer) {
                segment.boundingBox.x = xpos;
                segment.boundingBox.width = rightLayerPosition - xpos;
                segment.xDelta = xSegmentDelta;
                segment.isWestOfInitialLayer = leftLayer == null;
            }
            allSplineSegments.addAll(splineSegmentsLayer);

            // proceed to the next layer
            xpos = rightLayerPosition;
            if (rightLayer != null) {
                xpos += rightLayer.getSize().x;
            }
            
            leftLayer = rightLayer;
            isSpecialLeftLayer = isSpecialRightLayer;
        } while (rightLayer != null);
        
        // all layers have been processed, remember the spline paths for
        //  control point calculation to be done by a later intermediate processor
        for (LEdge edge : startEdges) {
            List<LEdge> edgeChain = getEdgeChain(edge);
            edge.setProperty(InternalProperties.SPLINE_EDGE_CHAIN, edgeChain);

            List<SplineSegment> spline = getSplinePath(edge);
            edge.setProperty(InternalProperties.SPLINE_ROUTE_START, spline);
        }
        
        // assign final width of the layering and thus the overall graph
        layeredGraph.getSize().x = xpos;
        
        lGraph = null;
        monitor.done();
    }

    private void createSegmentsAndComputeRanking() {
        // create the hyperEdges having their start port on the left side.
        createSplineSegmentsForHyperEdges(leftPortsLayer, rightPortsLayer, SideToProcess.LEFT, 
                true, edgesRemainingLayer, splineSegmentsLayer);
        createSplineSegmentsForHyperEdges(leftPortsLayer, rightPortsLayer, SideToProcess.LEFT, 
                false, edgesRemainingLayer, splineSegmentsLayer);

        // create the hyperEdges having their start port on the right side.
        createSplineSegmentsForHyperEdges(leftPortsLayer, rightPortsLayer, SideToProcess.RIGHT, 
                true, edgesRemainingLayer, splineSegmentsLayer);
        createSplineSegmentsForHyperEdges(leftPortsLayer, rightPortsLayer, SideToProcess.RIGHT, 
                false, edgesRemainingLayer, splineSegmentsLayer);

        // remaining edges are single edges that cannot be combined with others to a hyper-edge
        createSplineSegments(edgesRemainingLayer, leftPortsLayer, rightPortsLayer, splineSegmentsLayer);
        
        ////////////////////////////////////
        // Creation of the dependencies of the spline segments
        final ListIterator<SplineSegment> sourceIter = splineSegmentsLayer.listIterator();
        while (sourceIter.hasNext()) {
            final SplineSegment hyperEdge1 = sourceIter.next();
            final ListIterator<SplineSegment> targetIter = splineSegmentsLayer.listIterator(sourceIter.nextIndex());
            while (targetIter.hasNext()) {
                final SplineSegment hyperEdge2 = targetIter.next();
                createDependency(hyperEdge1, hyperEdge2);
            }
        }

        ////////////////////////////////////
        // Apply the topological numbering
        // break cycles
        breakCycles(splineSegmentsLayer, lGraph.getProperty(InternalProperties.RANDOM));
        
        // assign ranks to the hyper-nodes
        topologicalNumbering(splineSegmentsLayer);
    }
    
    /**
     * Initially fills the mappings, collection and sets for a pair of layers.
     * 
     * @param leftLayer The current left layer.
     * @param RightLayer The current right layer.
     */
    private void clearThenFillMappings(final Layer leftLayer, final Layer rightLayer) {

        // clear the mappings
        leftPortsLayer.clear();
        rightPortsLayer.clear();
        edgesRemainingLayer.clear();
        splineSegmentsLayer.clear();
        selfLoopsLayer.clear();

        // iterate over all outgoing edges on the left layer.
        if (leftLayer != null) {
            for (final LNode node : leftLayer.getNodes()) { 
                for (final LPort sourcePort : node.getPorts(PortSide.EAST)) {
                    leftPortsLayer.add(sourcePort);
                    
                    for (final LEdge edge : sourcePort.getOutgoingEdges()) {
                        // Self-loops are handled in the right-layer section below.
                        if (edge.isSelfLoop()) {
                            continue;
                        }
                        
                        // Add edge to set of all edges and find it's successor
                        edgesRemainingLayer.add(edge);  
                        findAndAddSuccessor(edge);

                        // Check if edge is a startingEdge
                        if (isQualifiedAsStartingNode(edge.getSource().getNode())) {
                            startEdges.add(edge);
                        }
                        
                        // Check port-side of target port
                        final LPort targetPort = edge.getTarget();
                        final Layer targetLayer = targetPort.getNode().getLayer();
                        if (targetLayer.equals(rightLayer)) {
                            rightPortsLayer.add(targetPort);
                        } else if (targetLayer.equals(leftLayer)) {
                            leftPortsLayer.add(targetPort);
                        } else {
                            // Unhandled situation. Probably there are incoming and outgoing edges on
                            // the same port. This is not supported.
                            edgesRemainingLayer.remove(edge);
                        }
                    }
                }
            }
        }
        
        if (rightLayer != null) {
            for (final LNode node : rightLayer.getNodes()) {
                // handle all self-loops, no matter on witch port-side they are.
                for (final LPort port : node.getPorts()) {
                    for (final LEdge edge : port.getOutgoingEdges()) {
                        if (edge.isSelfLoop()) {
                            selfLoopsLayer.add(edge);
                        }
                    }
                }

                // iterate over all outgoing edges on the right layer
                for (final LPort sourcePort : node.getPorts(PortSide.WEST)) {
                    rightPortsLayer.add(sourcePort);
                    
                    for (final LEdge edge : sourcePort.getOutgoingEdges()) {
                        // self-loops have been handled before
                        if (edge.isSelfLoop()) {
                            continue;
                        }
                        
                        // Add edge to set of all edges and find it's successor
                        edgesRemainingLayer.add(edge);
                        findAndAddSuccessor(edge);

                        // Check if edge is a startingEdge
                        if (isQualifiedAsStartingNode(edge.getSource().getNode())) {
                            startEdges.add(edge);
                        }
                        
                        // Check port-side of target port
                        final LPort targetPort = edge.getTarget();
                        final Layer targetLayer = targetPort.getNode().getLayer();
                        if (targetLayer.equals(rightLayer)) {
                            rightPortsLayer.add(targetPort);
                        } else if (targetLayer.equals(leftLayer)) {
                            leftPortsLayer.add(targetPort);
                        } else {
                            // Unhandled situation. Probably there are incoming and outgoing edges on
                            // the same port. This is not supported.
                            edgesRemainingLayer.remove(edge);
                        }
                    }
                }
            }
        }
    }
    
    /**
     * For sloppy routing the idea is to limit the angle of vertical segments similar to what is done in the
     * {@link PolylineEdgeRouter}. Consequently, the computed spacing depends on the maximum vertical span of any edge
     * between {@code rightLayer} and the preceding layer.
     */
    private double computeSloppySpacing(final Layer rightLayer, final double edgeEdgeSpacing,
            final double nodeNodeSpacing, final double sloppyLayerSpacingFactor) {
        
        double maxVertDiff = 0.0;
        // Iterate over the layer's nodes
        for (LNode node : rightLayer) {
            // Calculate the maximal vertical span of output edges.
            double maxCurrInputYDiff = 0.0;
            for (LEdge incomingEdge : node.getIncomingEdges()) {
                double sourcePos = incomingEdge.getSource().getAbsoluteAnchor().y;
                double targetPos = incomingEdge.getTarget().getAbsoluteAnchor().y;

                maxCurrInputYDiff = Math.max(maxCurrInputYDiff, Math.abs(targetPos - sourcePos));
            }
            maxVertDiff = Math.max(maxVertDiff, maxCurrInputYDiff);
        }

        // Determine where next layer should start based on the maximal vertical span of edges
        // between the two layers
        double layerSpacing =
                sloppyLayerSpacingFactor * Math.min(1.0, edgeEdgeSpacing / nodeNodeSpacing) * maxVertDiff;
        
        return layerSpacing;
    }
    
    /**
     * Finds the predecessor {@link LEdge} of given edge. It is assumed that the source node of an
     * edge with a predecessor only has one incoming edge. Otherwise the first incoming edge of the
     * source node is added as the predecessor.  
     * 
     * @param edge The {@link LEdge} those predecessor to find.
     */
    private void findAndAddSuccessor(final LEdge edge) {
        final LNode targetNode = edge.getTarget().getNode();
        
        // if target node is a normal node there is no successor
        if (isNormalNode(targetNode)) {
            return;
        }
        
        // otherwise take the first outgoing edge of target node
        final Iterator<LEdge> iter = targetNode.getOutgoingEdges().iterator();
        if (iter.hasNext()) {
            successingEdge.put(edge, iter.next());
        }
    }
    
    /**
     * Creates a "one-edge" hyper-edge for each edges in the collection. The hyper-edges are added to the hyperEdges
     * collection.
     * 
     * @param edges
     *            The edges to process.
     * @param leftPorts
     *            The left ports of the current situation.
     * @param rightPorts
     *            The right ports of the current situation.
     * @param hyperEdges
     *            The new hyper-edges will be added to this collection.
     */
    private void createSplineSegments(final List<LEdge> edges, 
            final Set<LPort> leftPorts, final Set<LPort> rightPorts,
            final List<SplineSegment> hyperEdges) {
        
        for (final LEdge edge : edges) {
            final LPort sourcePort = edge.getSource();
            SideToProcess sourceSide;
            
            if (leftPorts.contains(sourcePort)) {
                sourceSide = SideToProcess.LEFT;
            } else if (rightPorts.contains(sourcePort)) {
                sourceSide = SideToProcess.RIGHT;
            } else {
                throw new IllegalArgumentException("Source port must be in one of the port sets.");
            }

            SideToProcess targetSide;
            final LPort targetPort = edge.getTarget();

            if (leftPorts.contains(targetPort)) {
                targetSide = SideToProcess.LEFT;
            } else if (rightPorts.contains(targetPort)) {
                targetSide = SideToProcess.RIGHT;
            } else {
                throw new IllegalArgumentException("Target port must be in one of the port sets.");
            }

            SplineSegment seg = new SplineSegment(edge, sourceSide, targetSide);
            edgeToSegmentMap.put(edge, seg);
            hyperEdges.add(seg);
        }
    }

    /**
     * Creates hyperEdges. The created hyperEdges all have one port on their source side (if reversed
     * is {@code false}) or on their target side (if reversed is {@code true}). Also only hyperEdges 
     * starting in one of the leftPorts (if sideToProcess is {@code LEFT}) or one of the right ports (if 
     * SideToProcess is {@code RIGHT}) will be created.
     * 
     * @param leftPorts The ports on the left side of current between-layer segment.
     * @param rightPorts The ports on the right side of current between-layer segment.
     * @param sideToProcess Either {@code LEFT} or {@code RIGHT}. 
     * @param reversed {@code true}, if hyperEdges for reversed edges shall be created.
     * @param edgesRemaining Only hyperEdges pointing to a set of ports in this collection will be 
     *          created. 
     * @param hyperEdges The collection of hyperEdges that the created edges will be added to.
     */
    private void createSplineSegmentsForHyperEdges(
            final Set<LPort> leftPorts, 
            final Set<LPort> rightPorts, 
            final SideToProcess sideToProcess,
            final boolean reversed,
            final List<LEdge> edgesRemaining,
            final List<SplineSegment> hyperEdges) {

        Set<LPort> portsToProcess = null;
        if (sideToProcess == SideToProcess.LEFT) {
            portsToProcess = leftPorts;  
        } else if (sideToProcess == SideToProcess.RIGHT) {
            portsToProcess = rightPorts;
        } else {
            assert false : "sideToProcess must be either LEFT or RIGHT.";
        }
        
        // Iterate through all ports on the side to process.
        for (final LPort singlePort : portsToProcess) {
            final double singlePortPosition = singlePort.getAbsoluteAnchor().y;
            final Set<Pair<SideToProcess, LEdge>> upEdges = Sets.newHashSet();
            final Set<Pair<SideToProcess, LEdge>> downEdges = Sets.newHashSet();
            
            // Find edges we could construct a hyper-edge from. If the edge is in the 
            // edgesRemaining set, there is no hyper-edge that represents this edge. 
            for (final LEdge edge : singlePort.getConnectedEdges()) {
                if (edge.getProperty(InternalProperties.REVERSED) != reversed) {
                    continue;
                }
                if (edgesRemaining.contains(edge)) {
                    // find the target port
                    LPort targetPort;
                    if (edge.getTarget() == singlePort) {
                        targetPort = edge.getSource();
                    } else {
                        targetPort = edge.getTarget();
                    }
                    
                    // check if this edge should get drawn as a straight edge
                    final double targetPortPosition = targetPort.getAbsoluteAnchor().y;
                    if (isStraight(targetPortPosition, singlePortPosition)) {
                        continue;
                    }

                    // add the edge to the correct set of up/down-edges 
                    if (targetPortPosition < singlePortPosition) {
                        if (leftPorts.contains(targetPort)) {
                            upEdges.add(Pair.of(SideToProcess.LEFT, edge));
                        } else {
                            upEdges.add(Pair.of(SideToProcess.RIGHT, edge));
                        }
                    } else {
                        if (leftPorts.contains(targetPort)) {
                            downEdges.add(Pair.of(SideToProcess.LEFT, edge));
                        } else {
                            downEdges.add(Pair.of(SideToProcess.RIGHT, edge));
                        }
                    }
                }
            }

            // Create some hyper edges. 
            // We are creating only hyper-edges that have more than one real edge.  
            if (upEdges.size() > 1) {
                SplineSegment seg = new SplineSegment(singlePort, upEdges, sideToProcess);
                upEdges.forEach(e -> edgeToSegmentMap.put(e.getSecond(), seg));
                hyperEdges.add(seg);
                for (final Pair<SideToProcess, LEdge> pair : upEdges) {
                    edgesRemaining.remove(pair.getSecond());
                }
            }
            if (downEdges.size() > 1) {
                SplineSegment seg = new SplineSegment(singlePort, downEdges, sideToProcess);
                downEdges.forEach(e -> edgeToSegmentMap.put(e.getSecond(), seg));
                hyperEdges.add(seg);
                for (final Pair<SideToProcess, LEdge> pair : downEdges) {
                    edgesRemaining.remove(pair.getSecond());
                }
            }
        }
    }
    
    /**
     * Calculate the "must lay left of" dependency for two SplineHyperEdges.
     * @param edge0 First hyper-edge to compare.
     * @param edge1 Second hyper-edge to compare.
     */
    private void createDependency(final SplineSegment edge0, final SplineSegment edge1) {
        if (edge0.hyperEdgeTopYPos > edge1.hyperEdgeBottomYPos 
                || edge1.hyperEdgeTopYPos > edge0.hyperEdgeBottomYPos) {
            // the two hyper-edges do not share a vertical segment
            return;
        }
        int edge0Counter = 0;
        int edge1Counter = 0;
        
        for (final LPort port : edge0.rightPorts) {
            if (SplinesMath.isBetween(port.getAbsoluteAnchor().y, 
                    edge1.hyperEdgeTopYPos, edge1.hyperEdgeBottomYPos)) {
                edge0Counter++;
            }
        }
        for (final LPort port : edge0.leftPorts) {
            if (SplinesMath.isBetween(port.getAbsoluteAnchor().y, 
                    edge1.hyperEdgeTopYPos, edge1.hyperEdgeBottomYPos)) {
                edge0Counter--;
            }
        }
        for (final LPort port : edge1.rightPorts) {
            if (SplinesMath.isBetween(port.getAbsoluteAnchor().y, 
                    edge0.hyperEdgeTopYPos, edge0.hyperEdgeBottomYPos)) {
                edge1Counter++;
            }
        }
        for (final LPort port : edge1.leftPorts) {
            if (SplinesMath.isBetween(port.getAbsoluteAnchor().y, 
                    edge0.hyperEdgeTopYPos, edge0.hyperEdgeBottomYPos)) {
                edge1Counter--;
            }
        }
        
        if (edge0Counter < edge1Counter) {
            // edge0 should lay left of edge1
            new Dependency(edge0, edge1, edge1Counter - edge0Counter);
        } else if (edge1Counter < edge0Counter) {
            // edge0 should lay right of edge1
            new Dependency(edge1, edge0, edge0Counter - edge1Counter);
        } else {
            // in either ordering there would be the same number of crossings
            new Dependency(edge1, edge0, 0);
            new Dependency(edge0, edge1, 0);
        }
    }
    
    ///////////////////////////////////////////////////////////////////////////////
    // Cycle Breaking
    
    /**
     * Breaks all cycles in the given hypernode structure by reversing or removing
     * some dependencies. This implementation assumes that the dependencies of zero
     * weight are exactly the two-cycles of the hypernode structure.
     * 
     * @param edges list of hypernodes
     */
    private static void breakCycles(final List<SplineSegment> edges, final Random random) {
        final LinkedList<SplineSegment> sources = Lists.newLinkedList();
        final LinkedList<SplineSegment> sinks = Lists.newLinkedList();
        
        // initialize values for the algorithm
        int nextMark = -1;
        for (final SplineSegment edge : edges) {
            edge.mark = nextMark--;
            int inweight = 0;
            int outweight = 0;
            
            for (final Dependency dependency : edge.outgoing) {
                outweight += dependency.weight;
            }
            
            for (final Dependency dependency : edge.incoming) {
                inweight += dependency.weight;
            }
            
            edge.inweight = inweight;
            edge.outweight = outweight;
            
            if (outweight == 0) {
                sinks.add(edge);
            } else if (inweight == 0) {
                sources.add(edge);
            }
        }
    
        // assign marks to all nodes, ignore dependencies of weight zero
        final Set<SplineSegment> unprocessed = Sets.newLinkedHashSet(edges);
        final int markBase = edges.size();
        int nextLeft = markBase + 1;
        int nextRight = markBase - 1;
        final List<SplineSegment> maxEdges = Lists.newArrayList();

        while (!unprocessed.isEmpty()) {
            while (!sinks.isEmpty()) {
                final SplineSegment sink = sinks.removeFirst();
                unprocessed.remove(sink);
                sink.mark = nextRight--;
                updateNeighbors(sink, sources, sinks);
            }
            
            while (!sources.isEmpty()) {
                final SplineSegment source = sources.removeFirst();
                unprocessed.remove(source);
                source.mark = nextLeft++;
                updateNeighbors(source, sources, sinks);
            }
            
            int maxOutflow = Integer.MIN_VALUE;
            for (final SplineSegment edge : unprocessed) {
                final int outflow = edge.outweight - edge.inweight;
                if (outflow >= maxOutflow) {
                    if (outflow > maxOutflow) {
                        maxEdges.clear();
                        maxOutflow = outflow;
                    }
                    maxEdges.add(edge);
                }
            }
            
            if (!maxEdges.isEmpty()) {
                // if there are multiple SplineHyperEdges with maximal outflow, select one randomly
                final SplineSegment maxEdge = maxEdges.get(random.nextInt(maxEdges.size()));
                unprocessed.remove(maxEdge);
                maxEdge.mark = nextLeft++;
                updateNeighbors(maxEdge, sources, sinks);
                maxEdges.clear();
            }
        }
    
        // shift ranks that are left of the mark base
        final int shiftBase = edges.size() + 1;
        for (final SplineSegment edge : edges) {
            if (edge.mark < markBase) {
                edge.mark += shiftBase;
            }
        }
    
        // process edges that point left: remove those of zero weight, reverse the others
        for (final SplineSegment source : edges) {
            final ListIterator<Dependency> depIter = source.outgoing.listIterator();
            while (depIter.hasNext()) {
                final Dependency dependency = depIter.next();
                final SplineSegment target = dependency.target;
                
                if (source.mark > target.mark) {
                    depIter.remove();
                    target.incoming.remove(dependency);
                    
                    if (dependency.weight > 0) {
                        dependency.source = target;
                        target.outgoing.add(dependency);
                        dependency.target = source;
                        source.incoming.add(dependency);
                    }
                }
            }
        }
    }
    
    /**
     * Updates in-weight and out-weight values of the neighbors of the given node,
     * simulating its removal from the graph. The sources and sinks lists are
     * also updated.
     * 
     * @param edge node for which neighbors are updated
     * @param sources list of sources
     * @param sinks list of sinks
     */
    private static void updateNeighbors(final SplineSegment edge, 
            final List<SplineSegment> sources,
            final List<SplineSegment> sinks) {
        // process following edges
        for (final Dependency dep : edge.outgoing) {
            if (dep.target.mark < 0 && dep.weight > 0) {
                dep.target.inweight -= dep.weight;
                if (dep.target.inweight <= 0 && dep.target.outweight > 0) {
                    sources.add(dep.target);
                }
            }
        }
        
        // process preceding edges
        for (final Dependency dep : edge.incoming) {
            if (dep.source.mark < 0 && dep.weight > 0) {
                dep.source.outweight -= dep.weight;
                if (dep.source.outweight <= 0 && dep.source.inweight > 0) {
                    sinks.add(dep.source);
                }
            }
        }
    }
        
    ///////////////////////////////////////////////////////////////////////////////
    // Topological Ordering
    
    /**
     * Perform a topological numbering of the given SplineHyperEdges.
     * 
     * @param edges list of SplineHyperEdge
     */
    private static void topologicalNumbering(final List<SplineSegment> edges) {
        // determine sources, targets, incoming count and outgoing count; targets are only
        // added to the list if they only connect westward ports (that is, if all their
        // horizontal segments point to the right)
        final List<SplineSegment> sources = Lists.newLinkedList();
        final List<SplineSegment> rightwardTargets = Lists.newLinkedList();
        for (final SplineSegment edge : edges) {
            edge.rank = 0;
            edge.inweight = edge.incoming.size();
            edge.outweight = edge.outgoing.size();
            
            if (edge.inweight == 0) {
                sources.add(edge);
            }
            
            if (edge.outweight == 0 && edge.leftPorts.isEmpty()) {
                rightwardTargets.add(edge);
            }
        }
        
        int maxRank = -1;
        
        // assign ranks using topological numbering
        while (!sources.isEmpty()) {
            final SplineSegment edge = sources.remove(0);
            for (final Dependency dep : edge.outgoing) {
                SplineSegment target = dep.target;
                target.rank = Math.max(target.rank, edge.rank + 1);
                maxRank = Math.max(maxRank, target.rank);
                
                target.inweight--;
                if (target.inweight == 0) {
                    sources.add(target);
                }
            }
        }
        
        /* If we stopped here, hyper nodes that don't have any horizontal segments pointing
         * leftward would be ranked just like every other hyper node. This would move back
         * edges too far away from their target node. To remedy that, we move all hyper nodes
         * with horizontal segments only pointing rightwards as far right as possible.
         */
        if (maxRank > -1) {
            // assign all target nodes with horizontal segments pointing to the right the
            // rightmost rank
            for (final SplineSegment edge : rightwardTargets) {
                edge.rank = maxRank;
            }
            
            // let all other segments with horizontal segments pointing rightwards move as
            // far right as possible
            while (!rightwardTargets.isEmpty()) {
                final SplineSegment edge = rightwardTargets.remove(0);
                
                // The node only has connections to western ports
                for (final Dependency dep : edge.incoming) {
                    SplineSegment source = dep.source;
                    if (!source.leftPorts.isEmpty()) {
                        continue;
                    }
                    
                    source.rank = Math.min(source.rank, edge.rank - 1);
                    
                    source.outweight--;
                    if (source.outweight == 0) {
                        rightwardTargets.add(source);
                    }
                }
            }
        }
    }
    
    ///////////////////////////////////////////////////////////////////////////////
    // Convenience

    /**
     * Checks if the two Y coordinates are resulting in a straight edge.
     * @param firstY First Y coordinate.
     * @param secondY Second Y coordinate.
     * @return True, if the two Y coordinates result in a straight edge. 
     */
    static boolean isStraight(final double firstY, final double secondY) {
        return Math.abs(firstY - secondY) < MAX_VERTICAL_DIFF_FOR_STRAIGHT;
    }
    
    private List<LEdge> getEdgeChain(final LEdge start) {
        List<LEdge> edgeChain = Lists.newArrayList();
        LEdge current = start;
        do {
            edgeChain.add(current);
            current = successingEdge.get(current);
        } while (current != null);
        return edgeChain;
    }
    
    private List<SplineSegment> getSplinePath(final LEdge start) {
        List<SplineSegment> segmentChain = Lists.newArrayList();
        LEdge current = start;
        do {
            SplineSegment segment = edgeToSegmentMap.get(current);
            segment.sourcePort = current.getSource();
            segment.targetPort = current.getTarget();
            segmentChain.add(segment);
            current = successingEdge.get(current);
        } while (current != null);
        
        SplineSegment initialSegment = segmentChain.get(0);
        initialSegment.initialSegment = true;
        initialSegment.sourceNode = initialSegment.edges.iterator().next().getSource().getNode(); 
                
        SplineSegment lastSegment = segmentChain.get(segmentChain.size() - 1);
        lastSegment.lastSegment = true;
        lastSegment.targetNode = lastSegment.edges.iterator().next().getTarget().getNode();
        
        return segmentChain;
    }
    
    /**
     * Additionally to the {@link NodeType#NORMAL} nodes, we deem {@link NodeType#BIG_NODE} 
     * and {@link NodeType#BREAKING_POINT} to be normal nodes.
     * 
     * @param node the node to test
     * @return true if we consider {@code node} to be normal.
     */
    public static boolean isNormalNode(final LNode node) {
        NodeType nt = node.getType();
        return nt == NodeType.NORMAL 
            || nt == NodeType.BREAKING_POINT;
    }

    /**
     * A node {@code n} qualifies as 'starting node' for a 'starting edge' if it is not a long edge dummy or label
     * dummy. That is, any edge connected to {@code n} can potentially be the first segment of a chain of edges (i.e.
     * long edge).
     * 
     * @param node
     *            the node to test
     * @return true if {@code node} qualifies as starting node.
     */
    public static boolean isQualifiedAsStartingNode(final LNode node) {
        NodeType nt = node.getType();
        return nt == NodeType.NORMAL 
            || nt == NodeType.NORTH_SOUTH_PORT 
            || nt == NodeType.EXTERNAL_PORT 
            || nt == NodeType.BREAKING_POINT;
    }
    
    /**
     * A dependency between two {@link SplineSegment}s that represent hyperedges.
     * A dependency pointing from edge A to edge B means that edge A must lay left of edge B to
     * minimize the number of edge crossings.
     * A dependency with the weight of 0 means that the number of edge crossings does not vary, if the
     * position of the two hyper-edges is swapped. BUT the two hyper-edges share a vertical segment,
     * so they may not lay on the same x-coordinate as this would result in an overlapping segment. 
     */
    final class Dependency {
        /** The source of the dependency. (Should lay left) */
        private SplineSegment source;
        /** The target of the dependency. (Should lay right) */
        private SplineSegment target;
        /** The weight of the dependency. */
        private final int weight;
        
        /**
         * Creates a dependency from the given source to the given target.
         * 
         * @param source the dependency source
         * @param target the dependency target
         * @param weight weight of the dependency
         */
        Dependency(final SplineSegment source, final SplineSegment target, final int weight) {
            this.source = source;
            this.target = target;
            this.weight = weight;
            source.outgoing.add(this);
            target.incoming.add(this);
        }

        @Override
        public String toString() {
            return source + " ->(" + weight + ") " + target;
        }
    }
    
    /**
     * Pretty small enumeration used to define which side to process.
     */
    enum SideToProcess {
        /** Process the left side. */
        LEFT,
        /** Process the right side. */
        RIGHT
    }
}
