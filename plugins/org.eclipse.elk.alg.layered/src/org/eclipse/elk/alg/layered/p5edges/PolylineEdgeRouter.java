/*******************************************************************************
 * Copyright (c) 2010, 2020 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.layered.p5edges;

import java.util.ListIterator;
import java.util.Set;

import org.eclipse.elk.alg.layered.LayeredPhases;
import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LGraphUtil;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LNode.NodeType;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.intermediate.IntermediateProcessorStrategy;
import org.eclipse.elk.alg.layered.options.GraphProperties;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.core.alg.ILayoutPhase;
import org.eclipse.elk.core.alg.LayoutProcessorConfiguration;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.math.KVectorChain;
import org.eclipse.elk.core.options.PortSide;
import org.eclipse.elk.core.util.IElkProgressMonitor;

import com.google.common.base.Predicate;
import com.google.common.collect.Iterables;
import com.google.common.collect.Sets;

/**
 * Edge router module that draws edges with non-orthogonal line segments.
 * 
 * <dl>
 *   <dt>Precondition:</dt>
 *      <dd>the graph has a proper layering with assigned node and port positions</dd>
 *      <dd>the size of each layer is correctly set</dd>
 *      <dd>at least one of the nodes connected by an in-layer edge is a dummy node</dd>
 *   <dt>Postcondition:</dt>
 *      <dd>each node is assigned a horizontal coordinate</dd>
 *      <dd>the bend points of each edge are set</dd>
 *      <dd>the width of the whole graph is set</dd>
 * </dl>
 *
 * @author msp
 * @author cds
 */
public final class PolylineEdgeRouter implements ILayoutPhase<LayeredPhases, LGraph> {
    
    /**
     * Predicate that checks whether nodes represent external ports.
     */
    public static final Predicate<LNode> PRED_EXTERNAL_WEST_OR_EAST_PORT = new Predicate<LNode>() {
        public boolean apply(final LNode node) {
            PortSide extPortSide = node.getProperty(InternalProperties.EXT_PORT_SIDE);
            return node.getType() == NodeType.EXTERNAL_PORT
                    && (extPortSide == PortSide.WEST || extPortSide == PortSide.EAST);
        }
    };
    
    /* The basic processing strategy for this phase is empty. Depending on the graph features,
     * dependencies on intermediate processors are added dynamically as follows:
     * 
     * Before phase 1:
     *   - None.
     * 
     * Before phase 2:
     *   - For center edge labels:
     *     - LABEL_DUMMY_INSERTER
     * 
     * Before phase 3:
     *   - Always:
     *     - INVERTED_PORT_PROCESSOR
     *   
     *   - For non-free ports:
     *     - NORTH_SOUTH_PORT_PREPROCESSOR
     *     
     *   - For edge labels:
     *     - LABEL_SIDE_SELECTOR
     *   
     *   - For center edge labels:
     *     - LABEL_DUMMY_SWITCHER
     * 
     * Before phase 4:
     *   - For center edge labels:
     *     - LABEL_SIDE_SELECTOR
     *     
     *   - For end edge labels:
     *     - END_LABEL_PREPROCESSOR
     * 
     * Before phase 5:
     *   - None.
     * 
     * After phase 5:
     *   - For non-free ports:
     *     - NORTH_SOUTH_PORT_POSTPROCESSOR
     *     
     *   - For center edge labels:
     *     - LABEL_DUMMY_REMOVER
     *     
     *   - For end edge labels:
     *     - END_LABEL_POSTPROCESSOR
     */
    
    /** baseline config. */
    private static final LayoutProcessorConfiguration<LayeredPhases, LGraph> BASELINE_PROCESSOR_CONFIGURATION =
        LayoutProcessorConfiguration.<LayeredPhases, LGraph>create()
            .addBefore(LayeredPhases.P3_NODE_ORDERING, IntermediateProcessorStrategy.INVERTED_PORT_PROCESSOR);
    
    /** additional processor dependencies for graphs with northern / southern non-free ports. */
    private static final LayoutProcessorConfiguration<LayeredPhases, LGraph> NORTH_SOUTH_PORT_PROCESSING_ADDITIONS =
        LayoutProcessorConfiguration.<LayeredPhases, LGraph>create()
            .addBefore(LayeredPhases.P3_NODE_ORDERING, IntermediateProcessorStrategy.NORTH_SOUTH_PORT_PREPROCESSOR)
            .addAfter(LayeredPhases.P5_EDGE_ROUTING, IntermediateProcessorStrategy.NORTH_SOUTH_PORT_POSTPROCESSOR);
    
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
    
    /** additional processor dependencies for graphs with head or tail edge labels. */
    private static final LayoutProcessorConfiguration<LayeredPhases, LGraph> END_EDGE_LABEL_PROCESSING_ADDITIONS =
        LayoutProcessorConfiguration.<LayeredPhases, LGraph>create()
            .addBefore(LayeredPhases.P4_NODE_PLACEMENT, IntermediateProcessorStrategy.LABEL_SIDE_SELECTOR)
            .addBefore(LayeredPhases.P4_NODE_PLACEMENT, IntermediateProcessorStrategy.END_LABEL_PREPROCESSOR)
            .addAfter(LayeredPhases.P5_EDGE_ROUTING, IntermediateProcessorStrategy.END_LABEL_POSTPROCESSOR);
    
    
    @Override
    public LayoutProcessorConfiguration<LayeredPhases, LGraph> getLayoutProcessorConfiguration(final LGraph graph) {
        Set<GraphProperties> graphProperties = graph.getProperty(InternalProperties.GRAPH_PROPERTIES);
        
        // Basic configuration
        LayoutProcessorConfiguration<LayeredPhases, LGraph> configuration =
                LayoutProcessorConfiguration.createFrom(BASELINE_PROCESSOR_CONFIGURATION);
        
        // Additional dependencies
        if (graphProperties.contains(GraphProperties.NORTH_SOUTH_PORTS)) {
            configuration.addAll(NORTH_SOUTH_PORT_PROCESSING_ADDITIONS);
        }
        
        if (graphProperties.contains(GraphProperties.SELF_LOOPS)) {
            configuration.addAll(SELF_LOOP_PROCESSING_ADDITIONS);
        }
        
        if (graphProperties.contains(GraphProperties.CENTER_LABELS)) {
            configuration.addAll(CENTER_EDGE_LABEL_PROCESSING_ADDITIONS);
        }
        
        if (graphProperties.contains(GraphProperties.END_LABELS)) {
            configuration.addAll(END_EDGE_LABEL_PROCESSING_ADDITIONS);
        }
        
        return configuration;
    }

    
    /** the minimal vertical difference for creating bend points. */
    private static final double MIN_VERT_DIFF = 1.0;
    /** factor for spacing apart layers between which edges are routed. */
    private static final double LAYER_SPACE_FAC = 0.4;

    /** Set of already created junction points, to avoid multiple points at the same position. */
    private final Set<KVector> createdJunctionPoints = Sets.newHashSet();
    
    
    /* Implementation Note:
     * 
     * The process method works by going through each layer and possibly adding bend points to incoming
     * edges of nodes and to outgoing edges of nodes to avoid edge-node-overlaps. While other edge
     * routing algorithms insert all bend points required to route edges between a pair of layers, this
     * router actually routes any edge in two steps: first while iterating through the source layer, and
     * second while iterating through the target layer. Understanding causes the code below to make a
     * lot more sense.
     */

    @Override
    public void process(final LGraph layeredGraph, final IElkProgressMonitor monitor) {
        monitor.begin("Polyline edge routing", 1);
        
        final double slopedEdgeZoneWidth =
                layeredGraph.getProperty(LayeredOptions.EDGE_ROUTING_POLYLINE_SLOPED_EDGE_ZONE_WIDTH);
        final double nodeSpacing = layeredGraph.getProperty(LayeredOptions.SPACING_NODE_NODE_BETWEEN_LAYERS);
        final double edgeSpacing = layeredGraph.getProperty(LayeredOptions.SPACING_EDGE_EDGE_BETWEEN_LAYERS);
        final double edgeSpaceFac = Math.min(1.0, edgeSpacing / nodeSpacing);
        
        double xpos = 0.0;
        double layerSpacing = 0.0;
        
        // Determine the horizontal spacing required to route west-side in-layer edges of the first layer
        if (!layeredGraph.getLayers().isEmpty()) {
            double yDiff = calculateWestInLayerEdgeYDiff(layeredGraph.getLayers().get(0));
            xpos = LAYER_SPACE_FAC * edgeSpaceFac * yDiff;
        }
        
        // Iterate over the layers
        ListIterator<Layer> layerIter = layeredGraph.getLayers().listIterator();
        while (layerIter.hasNext()) {
            Layer layer = layerIter.next();
            boolean externalLayer = Iterables.all(layer, PRED_EXTERNAL_WEST_OR_EAST_PORT);
            
            // The rightmost layer is not given any node spacing if it's an external port layer
            if (externalLayer && xpos > 0) {
                xpos -= nodeSpacing;
            }
            
            // Set horizontal coordinates for all nodes of the layer
            LGraphUtil.placeNodesHorizontally(layer, xpos);
            
            // While routing edges, we remember the maximum vertical span of any edge between this and
            // the next layer to insert enough space between the layers to keep the edge slopes from
            // becoming too steep
            double maxVertDiff = 0.0;
            
            // Iterate over the layer's nodes
            for (LNode node : layer) {
                // Calculate the maximal vertical span of output edges. In-layer edges will also be
                // routed at this point
                double maxCurrOutputYDiff = 0.0;
                for (LEdge outgoingEdge : node.getOutgoingEdges()) {
                    double sourcePos = outgoingEdge.getSource().getAbsoluteAnchor().y;
                    double targetPos = outgoingEdge.getTarget().getAbsoluteAnchor().y;
                    
                    if (layer == outgoingEdge.getTarget().getNode().getLayer() && !outgoingEdge.isSelfLoop()) {
                        // In-layer edges require an extra bend point to make them look nice
                        processInLayerEdge(outgoingEdge, xpos,
                                LAYER_SPACE_FAC * edgeSpaceFac * Math.abs(sourcePos - targetPos));
                        
                        if (outgoingEdge.getSource().getSide() == PortSide.WEST) {
                            // The spacing required for routing in-layer edges on the west side doesn't
                            // contribute anything to the spacing required between this and the next
                            // layer and was already taken into account previously
                            sourcePos = 0;
                            targetPos = 0;
                        }
                    }
                    
                    maxCurrOutputYDiff = Math.max(maxCurrOutputYDiff, Math.abs(targetPos - sourcePos));
                }
                
                // We currently only handle certain node types. This might change in the future
                switch (node.getType()) {
                case NORMAL:
                case LABEL:
                case LONG_EDGE:
                case NORTH_SOUTH_PORT:
                case BREAKING_POINT:
                    processNode(node, xpos, slopedEdgeZoneWidth);
                    break;
                }
                
                maxVertDiff = Math.max(maxVertDiff, maxCurrOutputYDiff);
            }
            
            // Consider the span of west-side in-layer edges of the next layer to be sure to reserve
            // enough space for routing them during the next iteration
            if (layerIter.hasNext()) {
                double yDiff = calculateWestInLayerEdgeYDiff(layerIter.next());
                maxVertDiff = Math.max(maxVertDiff, yDiff);
                layerIter.previous();
            }
            
            // Determine where next layer should start based on the maximal vertical span of edges
            // between the two layers
            layerSpacing = LAYER_SPACE_FAC * edgeSpaceFac * maxVertDiff;
            if (!externalLayer && layerIter.hasNext()) {
                layerSpacing += nodeSpacing;
            }
            
            xpos += layer.getSize().x + layerSpacing;
        }
        
        createdJunctionPoints.clear();
        
        // Set the graph's horizontal size
        layeredGraph.getSize().x = xpos;
        
        monitor.done();
    }
    
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Actual Edge Routing Code

    /**
     * Inserts bend points for edges incident to this node. The bend points are inserted such that
     * the segments that cross the layer's area are straight or pretty much straight. A bend point
     * is only added if it's necessary: if the bend point to be inserted differs from the edge's end
     * point.
     * 
     * @param node
     *            the node whose incident edges to insert bend points for.
     * @param layerLeftXPos
     *            the x position of the node's layer.
     * @param maxAcceptableXDiff
     *            the maximum space between node and layer boundary before bend points need to be inserted.
     */
    private void processNode(final LNode node, final double layerLeftXPos, final double maxAcceptableXDiff) {
        // The right side of the layer
        final double layerRightXPos = layerLeftXPos + node.getLayer().getSize().x;
        
        for (LPort port : node.getPorts()) {
            KVector absolutePortAnchor = port.getAbsoluteAnchor();
            
            if (node.getType() == NodeType.NORTH_SOUTH_PORT) {
                // North/south ports require special handling (also see #515):
                // The dummy node that represents the north/south port will have been placed at:
                //  - the very left x-coordinate of the currently processed layer for incoming edges 
                //  - the very right x-coordinate for outgoing edges
                // However, the actual x-coordinate of the port itself is very likely to have a different x-coordinate.
                // The code below introduces bendpoints that prevent edge/node overlaps based on a condition that
                // checks that a node's position is "far enough away" from the left-most x-coordinate of the layer.
                // Consequently, we have to use the north/south port's x-coordinate instead of the dummy node's one.
                LPort correspondingPort = (LPort) port.getProperty(InternalProperties.ORIGIN);
                absolutePortAnchor.x = correspondingPort.getAbsoluteAnchor().x;
                // While it makes sense either way, it is important to move the dummy node to the correct location
                // as well. Otherwise the #addBendPoint(...) method won't add the computed bendpoint as it
                // thinks it's superfluous.
                node.getPosition().x = absolutePortAnchor.x;
            }
            
            KVector bendPoint = new KVector(0, absolutePortAnchor.y);
            
            if (port.getSide() == PortSide.EAST) {
                bendPoint.x = layerRightXPos;
            } else if (port.getSide() == PortSide.WEST) {
                bendPoint.x = layerLeftXPos;
            } else {
                // We only know what to do with eastern and western ports
                continue;
            }
            
            // If the port's absolute anchor equals the bend point, we don't want to insert anything
            // (unless the node represents an in-layer dummy)
            double xDistance = Math.abs(absolutePortAnchor.x - bendPoint.x);
            if (xDistance <= maxAcceptableXDiff && !isInLayerDummy(node)) {
                continue;
            }
            
            // Whether to add a junction point or not
            boolean addJunctionPoint =
                    port.getOutgoingEdges().size() + port.getIncomingEdges().size() > 1;
            
            // Iterate over the edges and add bend (and possibly junction) points
            for (LEdge e : port.getConnectedEdges()) {
                LPort otherPort = e.getSource() == port ? e.getTarget() : e.getSource();
                if (Math.abs(otherPort.getAbsoluteAnchor().y - bendPoint.y) > MIN_VERT_DIFF) {
                    // Insert bend point
                    addBendPoint(e, bendPoint, addJunctionPoint, port);
                }
            }

        }
    }

    /**
     * In-layer edges get an extra bend point in addition to the usual source and target bend
     * points. The additional bend point is inserted halfway between the edge's upper and lower end,
     * a little way from the layer's boundary.
     * 
     * @param edge
     *            the in-layer edge to route.
     * @param layerXPos
     *            the layer's x position.
     * @param edgeSpacing
     *            the spacing to respect for the in-layer edge bend points.
     */
    private void processInLayerEdge(final LEdge edge, final double layerXPos, final double edgeSpacing) {
        LPort sourcePort = edge.getSource();
        LPort targetPort = edge.getTarget();
        
        double sourceAnchorY = sourcePort.getAbsoluteAnchor().y;
        double midY = (sourceAnchorY + targetPort.getAbsoluteAnchor().y) / 2.0;
        
        /* This method is called if an outgoing in-layer edge is found before any other edges of the
         * offending node are routed. Thus, if the edge's list of bend points is not empty, any bend
         * point must be at the target port. It follows that our extra bend point can safely be
         * inserted at the start of the bend point list.
         */
        
        KVector bendPoint = null;
        if (sourcePort.getSide() == PortSide.EAST) {
            bendPoint = new KVector(
                    layerXPos + sourcePort.getNode().getLayer().getSize().x + edgeSpacing,
                    midY);
        } else {
            bendPoint = new KVector(layerXPos - edgeSpacing, midY);
        }
        
        edge.getBendPoints().add(0, bendPoint);
    }
    
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Utility Methods
    
    /**
     * Calculates the maximum vertical span of any in-layer edge connecting west-side ports of nodes
     * in the given layer.
     * 
     * @param layer
     *            the layer to iterate over.
     * @return maximum vertical span of west-side in-layer edges.
     */
    private double calculateWestInLayerEdgeYDiff(final Layer layer) {
        double maxYDiff = 0.0;
        
        for (LNode node : layer) {
            for (LEdge outgoingEdge : node.getOutgoingEdges()) {
                if (layer == outgoingEdge.getTarget().getNode().getLayer()
                        && outgoingEdge.getSource().getSide() == PortSide.WEST) {
                    
                    double sourcePos = outgoingEdge.getSource().getAbsoluteAnchor().y;
                    double targetPos = outgoingEdge.getTarget().getAbsoluteAnchor().y;
                    maxYDiff = Math.max(maxYDiff, Math.abs(targetPos - sourcePos));
                }
            }
        }
        
        return maxYDiff;
    }
    
    /**
     * Adds a copy of the given bend point to the given edge at the appropriate place in the list of
     * bend points. The appropriate place is determined by the given port: if it's the source port,
     * the bend point is prepended to the list; otherwise is is appended. The bend point is not
     * added to the list at all if it wouldn't have any visual effect; that is, if the port's anchor
     * equals the bend point.
     * 
     * @param edge
     *            the edge to add the bend point to.
     * @param bendPoint
     *            the bend point to add.
     * @param addJunctionPoint
     *            if {@code true}, a copy of the bend point will be added to the edge's list of
     *            junction points, if necessary.
     * @param currPort
     *            the port the bend point is near.
     */
    private void addBendPoint(final LEdge edge, final KVector bendPoint, final boolean addJunctionPoint,
            final LPort currPort) {
        
        // Only insert the bend point if necessary; for in-layer edges we are extra save and add the bend point in
        // any case
        if ((edge.isInLayerEdge() || !currPort.getAbsoluteAnchor().equals(bendPoint)) && !edge.isSelfLoop()) {
            if (edge.getSource() == currPort) {
                edge.getBendPoints().add(0, new KVector(bendPoint));
            } else {
                edge.getBendPoints().add(new KVector(bendPoint));
            }
            
            if (addJunctionPoint && !createdJunctionPoints.contains(bendPoint)) {
                // create a new junction point for the edge at the bend point's position
                KVectorChain junctionPoints = edge.getProperty(LayeredOptions.JUNCTION_POINTS);
                if (junctionPoints == null) {
                    junctionPoints = new KVectorChain();
                    edge.setProperty(LayeredOptions.JUNCTION_POINTS, junctionPoints);
                }
                
                KVector jpoint = new KVector(bendPoint);
                junctionPoints.add(jpoint);
                createdJunctionPoints.add(jpoint);
            }
        }
    }
    
    /**
     * At this point a node is considered to be an in-layer dummy 
     * if it is of type {@link NodeType#LONG_EDGE} and has an incident edge that is 
     * an in-layer edge. 
     */
    private boolean isInLayerDummy(final LNode node) {
        if (node.getType() == NodeType.LONG_EDGE) {
            for (LEdge e : node.getConnectedEdges()) {
                if (e.isInLayerEdge()) {
                    return true;
                }
            }
        }
        return false;
    }
}
