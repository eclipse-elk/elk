/*******************************************************************************
 * Copyright (c) 2010, 2015 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.alg.layered.intermediate;

import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.Set;
import java.util.ArrayList;
import java.util.LinkedHashSet;

import org.eclipse.elk.alg.common.nodespacing.NodeDimensionCalculation;
import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LGraphAdapters;
import org.eclipse.elk.alg.layered.graph.LGraphUtil;
import org.eclipse.elk.alg.layered.graph.LLabel;
import org.eclipse.elk.alg.layered.graph.LMargin;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LNode.NodeType;
import org.eclipse.elk.alg.layered.graph.LPadding;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.layered.p5edges.orthogonal.OrthogonalRoutingGenerator;
import org.eclipse.elk.alg.layered.p5edges.orthogonal.direction.RoutingDirection;
import org.eclipse.elk.core.alg.ILayoutProcessor;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.math.KVectorChain;
import org.eclipse.elk.core.options.PortConstraints;
import org.eclipse.elk.core.options.PortLabelPlacement;
import org.eclipse.elk.core.options.PortSide;
import org.eclipse.elk.core.options.SizeConstraint;
import org.eclipse.elk.core.util.IElkProgressMonitor;


/**
 * This processor does the job of routing edges connected to hierarchical ports.
 * 
 * <p>Without going too much into the details, it does so in six distinct steps:</p>
 * <ol>
 *   <li>Dummy nodes created for northern and southern hierarchical ports are restored and
 *     connected to the nodes formerly created in their stead.</li>
 *   <li>The coordinates of these dummy nodes - and thus of the hierarchical ports they
 *     represent - are calculated and applied.</li>
 *   <li>The edges connected to northern and southern hierarchical ports are routed.</li>
 *   <li>The temporary dummy nodes are removed, thereby restoring the original bijective
 *     relationship between external ports and the dummy nodes created to represent them.</li>
 *   <li>Due to the necessity of having to route these edges, additional height may be
 *     required, which in turn may invalidate the y coordinates of eastern and western
 *     hierarchical port dummy nodes. Those are corrected in this step.</li>
 *   </li>That may in turn have produced slanted edge segments, which are corrected in the
 *     final step.</li>
 * </ol>
 * 
 * <dl>
 *   <dt>Precondition:</dt>
 *     <dd>A layered graph</dd>
 *     <dd>with edge routing finished for edges not incident to external ports</dd>
 *     <dd>long edge dummies are not yet joined.</dd>
 *   <dt>Postcondition:</dt>
 *     <dd>All external port dummy nodes left map onto an actual external port</dd>
 *     <dd>the coordinates of external port dummy nodes specify the coordinates of their respective external port</dd>
 *     <dd>all external port dummy nodes have a size of (0, 0)</dd>
 *     <dd>edges connected to external ports have their bend points set.</dd>
 *   <dt>Slots:</dt>
 *     <dd>After phase 5.</dd>
 *   <dt>Same-slot dependencies:</dt>
 *     <dd>None.</dd>
 * </dl>
 * 
 * @see HierarchicalPortConstraintProcessor
 * @see HierarchicalPortDummySizeProcessor
 * @see HierarchicalPortPositionProcessor
 * @see OrthogonalRoutingGenerator
 * @author cds
 */
public final class HierarchicalPortOrthogonalEdgeRouter implements ILayoutProcessor<LGraph> {
    
    // VARIABLES
    
    /**
     * The amount of space necessary to accomodate northern external port edge routing.
     */
    private double northernExtPortEdgeRoutingHeight;

    @Override
    public void process(final LGraph layeredGraph, final IElkProgressMonitor monitor) {
        monitor.begin("Orthogonally routing hierarchical port edges", 1);
        northernExtPortEdgeRoutingHeight = 0.0;
        
        /* Step 1
         * Restore any north / south port dummies removed by the HierarchicalPortConstraintProcessor
         * and connect them to the dummies created in their stead.
         */
        List<LNode> northSouthDummies = restoreNorthSouthDummies(layeredGraph);
        
        /* Step 2
         * Calculate coordinates for the north / south port dummies. Coordinates for the
         * east / west port dummies have already been calculated prior to this processor's
         * execution. The coordinates are relative to the node's content area, just like
         * normal node coordinates. (the content area is the node size minus padding minus offset)
         */
        setNorthSouthDummyCoordinates(layeredGraph, northSouthDummies);
        
        /* Step 3
         * Orthogonal edge routing.
         */
        routeEdges(monitor, layeredGraph, northSouthDummies);
        
        /* Step 4
         * Removal of the temporarily created north / south port dummies.
         */
        removeTemporaryNorthSouthDummies(layeredGraph);
        
        /* Step 5
         * Finally, the coordinates of east / west hierarchical port dummies have to be corrected
         * and set. The x coordinate must be set, and if north / south port routing resulted
         * in a change of offset or graph size, the y coordinates have to be adjusted if port
         * constraints are at FIXED_RATIO or FIXED_POS. The graph's width may also have to be
         * adjusted.
         */
        fixCoordinates(layeredGraph);
        
        /* Step 6
         * Fixing the dummy coordinates can easily lead to slanted edge segments, which needs
         * to be corrected.
         */
        correctSlantedEdgeSegments(layeredGraph);
        
        monitor.done();
    }
    
    
    ///////////////////////////////////////////////////////////////////////////////
    // STEP 1: RESTORE NORTH / SOUTH DUMMIES
    
    /**
     * Iterates through all layers, restoring hierarchical port dummy nodes along the way. The
     * restored nodes are connected to the temporary dummy nodes created for them (or, as Carsten
     * calls them, the "dummy dummies"). The restored nodes are added to the last layer. (which layer
     * they are added to doesn't make any difference)
     * 
     * @param layeredGraph the layered graph.
     * @return the list of restored external port dummies.
     */
    private List<LNode> restoreNorthSouthDummies(final LGraph layeredGraph) {
        List<LNode> restoredDummies = new ArrayList<>();
        
        if (!layeredGraph.hasProperty(InternalProperties.EXT_PORT_REPLACED_DUMMIES)) {
            return restoredDummies;
        }
        
        // Restore the original external port dummies
        for (LNode dummy : layeredGraph.getProperty(InternalProperties.EXT_PORT_REPLACED_DUMMIES)) {
            restoreDummy(dummy, layeredGraph);
            restoredDummies.add(dummy);
        }
        
        // Looking for hierarchical port dummies that replaced the restored hierarchical port dummy
        for (Layer layer : layeredGraph) {
            for (LNode node : layer) {
                if (node.getType() != NodeType.EXTERNAL_PORT) {
                    // Not a hierarchical port dummy - we're not interested. Move along,
                    // please, there's nothing to see here.
                    continue;
                }
                
                LNode replacedDummy = (LNode) node.getProperty(InternalProperties.EXT_PORT_REPLACED_DUMMY);
                if (replacedDummy != null) {
                    assert replacedDummy.getType() == NodeType.EXTERNAL_PORT;
                    
                    // Restore the origin and connect the node to it
                    connectNodeToDummy(layeredGraph, node, replacedDummy);
                }
            }
        }
        
        // Assign the restored dummies to the graph's last layer (this has nothing to do anymore
        // with where they'll be placed; we just need some layer to assign them to, and we simply
        // choose the last one)
        for (LNode dummy : restoredDummies) {
            dummy.setLayer(layeredGraph.getLayers().get(layeredGraph.getLayers().size() - 1));
        }
        
        return restoredDummies;
    }
    
    /**
     * Restores the given dummy.
     */
    private void restoreDummy(final LNode dummy, final LGraph graph) {
        // Depending on the hierarchical port's side, we set the port side of the dummy's ports
        // to be able to route properly (northern dummies must have a southern port)
        PortSide portSide = dummy.getProperty(InternalProperties.EXT_PORT_SIDE);
        LPort dummyPort = dummy.getPorts().get(0);
        
        if (portSide == PortSide.NORTH) {
            dummyPort.setSide(PortSide.SOUTH);
        } else if (portSide == PortSide.SOUTH) {
            dummyPort.setSide(PortSide.NORTH);
        }
        
        // Since the dummy node was hidden from the algorithm, its port labels are not placed properly and its margins
        // are not set acoordingly. That needs to be fixed (if port labels are to be taken into consideration)
        if (graph.getProperty(LayeredOptions.NODE_SIZE_CONSTRAINTS).contains(SizeConstraint.PORT_LABELS)) {
            // The ElkGraphImporter has set the relevant spacings on the dummy node
            double portLabelSpacingHorizontal = dummy.getProperty(LayeredOptions.SPACING_LABEL_PORT_HORIZONTAL);
            double portLabelSpacingVertical = dummy.getProperty(LayeredOptions.SPACING_LABEL_PORT_VERTICAL);
            double labelLabelSpacing = dummy.getProperty(LayeredOptions.SPACING_LABEL_LABEL);
            
            Set<PortLabelPlacement> portLabelPlacement = graph.getProperty(LayeredOptions.PORT_LABELS_PLACEMENT);
            if (portLabelPlacement.contains(PortLabelPlacement.INSIDE)) {
                double currentY = portLabelSpacingVertical;
                double xCenterRelativeToPort = dummy.getSize().x / 2 - dummyPort.getPosition().x;
                
                for (LLabel label : dummyPort.getLabels()) {
                    label.getPosition().y = currentY;
                    label.getPosition().x = xCenterRelativeToPort - label.getSize().x / 2;
                    
                    currentY += label.getSize().y + labelLabelSpacing;
                }
                
            } else if (portLabelPlacement.contains(PortLabelPlacement.OUTSIDE)) {
                // Port labels have a vertical size of 0, but we need to set their x coordinate
                for (LLabel label : dummyPort.getLabels()) {
                    label.getPosition().x = portLabelSpacingHorizontal + dummy.getSize().x - dummyPort.getPosition().x;
                }
            }
            
            // Calculate margins
            NodeDimensionCalculation.getNodeMarginCalculator(LGraphAdapters.adapt(graph, false))
                .processNode(LGraphAdapters.adapt(dummy, false));
        }
    }
    
    /**
     * Adds a port to the given node and connects that to the given dummy node.
     * 
     * @param node the node to connect to the dummy.
     * @param dummy the external port dummy to connect the node to.
     */
    private void connectNodeToDummy(final LGraph layeredGraph, final LNode node, final LNode dummy) {
        // First, add a port to the node. The port side depends on the node's hierarchical port side
        LPort outPort = new LPort();
        outPort.setNode(node);
        
        PortSide extPortSide = node.getProperty(InternalProperties.EXT_PORT_SIDE);
        outPort.setSide(extPortSide);
        
        // Find the dummy node's port
        LPort inPort = dummy.getPorts().get(0);
        
        // Connect the two nodes
        LEdge edge = new LEdge();
        edge.setSource(outPort);
        edge.setTarget(inPort);
    }
    
    
    ///////////////////////////////////////////////////////////////////////////////
    // STEP 2: SET NORTH / SOUTH DUMMY COORDINATES
    
    /**
     * Set coordinates for northern and southern external port dummy nodes.
     * 
     * @param layeredGraph the layered graph.
     * @param northSouthDummies dummy nodes whose position to set.
     */
    private void setNorthSouthDummyCoordinates(final LGraph layeredGraph, final List<LNode> northSouthDummies) {
        PortConstraints constraints = layeredGraph.getProperty(LayeredOptions.PORT_CONSTRAINTS);
        KVector graphSize = layeredGraph.getSize();
        LPadding graphPadding = layeredGraph.getPadding();
        double graphWidth = graphSize.x + graphPadding.left + graphPadding.right;
        double northY = 0 - graphPadding.top - layeredGraph.getOffset().y;
        double southY = graphSize.y + graphPadding.top + graphPadding.bottom - layeredGraph.getOffset().y;
        
        // Lists of northern and southern external port dummies
        List<LNode> northernDummies = new ArrayList<>();
        List<LNode> southernDummies = new ArrayList<>();
        
        for (LNode dummy : northSouthDummies) {
            // Set x coordinate
            switch (constraints) {
            case FREE:
            case FIXED_SIDE:
            case FIXED_ORDER:
                calculateNorthSouthDummyPositions(dummy);
                break;
                
            case FIXED_RATIO:
                applyNorthSouthDummyRatio(dummy, graphWidth);
                dummy.borderToContentAreaCoordinates(true, false);
                break;
            
            case FIXED_POS:
                applyNorthSouthDummyPosition(dummy);
                dummy.borderToContentAreaCoordinates(true, false);
                
                // Ensure that the graph is wide enough to hold the port
                graphSize.x = Math.max(
                        graphSize.x,
                        dummy.getPosition().x + dummy.getSize().x / 2.0);
                break;
            }
            
            // Set y coordinates and add the dummy to its respective list
            switch (dummy.getProperty(InternalProperties.EXT_PORT_SIDE)) {
            case NORTH:
                dummy.getPosition().y = northY;
                northernDummies.add((dummy));
                break;
            
            case SOUTH:
                dummy.getPosition().y = southY;
                southernDummies.add(dummy);
                break;
            }
        }
        
        // Check for correct ordering of nodes in the FIXED_ORDER case and for dummy nodes that
        // have been put on top of one another
        switch (constraints) {
        case FREE:
        case FIXED_SIDE:
            ensureUniquePositions(northernDummies, layeredGraph);
            ensureUniquePositions(southernDummies, layeredGraph);
            break;
            
        case FIXED_ORDER:
            restoreProperOrder(northernDummies, layeredGraph);
            restoreProperOrder(southernDummies, layeredGraph);
            break;
        }
    }

    /**
     * Calculates the positions of northern and southern dummy nodes. The position is based
     * on the nodes the dummy nodes are connected to, which is taken as the position of the
     * port's center.
     * 
     * @param dummy the northern or southern external port dummy node to calculate the position for. 
     */
    private void calculateNorthSouthDummyPositions(final LNode dummy) {
        // We use a simple algorithm that simply adds the horizontal positions of all
        // connected ports and divides the sum by the number of connected ports
        
        // First, get the dummy's port (it has only one)
        LPort dummyInPort = dummy.getPorts().get(0);
        
        // It might be that the node is not actually connected to anything. In that case, we assign a default position
        if (dummyInPort.getDegree() == 0) {
            dummy.getPosition().x = 0;
        } else {
            // Now, iterate over all connected ports, adding their horizontal position
            double posSum = 0.0;
            
            for (LPort connectedPort : dummyInPort.getConnectedPorts()) {
                posSum += connectedPort.getNode().getPosition().x + connectedPort.getPosition().x
                        + connectedPort.getAnchor().x;
            }
            
            // Assign the dummy's x coordinate
            KVector anchor = dummy.getProperty(LayeredOptions.PORT_ANCHOR);
            double offset = anchor == null ? 0 : anchor.x;
            dummy.getPosition().x = posSum / dummyInPort.getDegree() - offset;
        }
    }
    
    /**
     * Sets the dummy's x coordinate to respect the ratio defined for its original port.
     * 
     * @param dummy the dummy.
     * @param width the graph width.
     */
    private void applyNorthSouthDummyRatio(final LNode dummy, final double width) {
        KVector anchor = dummy.getProperty(LayeredOptions.PORT_ANCHOR);
        double offset = anchor == null ? 0 : anchor.x;
        dummy.getPosition().x = width * dummy.getProperty(InternalProperties.PORT_RATIO_OR_POSITION)
                - offset;
    }
    
    /**
     * Sets the dummy's x coordinate to its original port's x coordinate.
     * 
     * @param dummy the dummy.
     */
    private void applyNorthSouthDummyPosition(final LNode dummy) {
        KVector anchor = dummy.getProperty(LayeredOptions.PORT_ANCHOR);
        double offset = anchor == null ? 0 : anchor.x;
        dummy.getPosition().x = dummy.getProperty(InternalProperties.PORT_RATIO_OR_POSITION) - offset;
    }
    
    /**
     * Ensures that no two dummy nodes in the given list are assigned the same x coordinate. This
     * method must only be called if port constraints are set to {@code FREE} or {@code FIXED_SIDE}
     * as it may not preserve the original order of dummy nodes.
     * 
     * @param dummies list of dummy nodes.
     * @param graph the layered graph.
     */
    private void ensureUniquePositions(final List<LNode> dummies, final LGraph graph) {
        if (dummies.isEmpty()) {
            return;
        }
        
        // Turn the list into an array of dummy nodes and sort that by their x coordinate
        LNode[] dummyArray = LGraphUtil.toNodeArray(dummies);
        
        Arrays.sort(dummyArray, new Comparator<LNode>() {
            public int compare(final LNode a, final LNode b) {
                return Double.compare(a.getPosition().x, b.getPosition().x);
            }
        });
        
        assignAscendingCoordinates(dummyArray, graph);
    }

    /**
     * Checks if the automatically calculated node coordinates violate their fixed order and fixes
     * the coordinates. Calling this method only makes sense if port constraints are set to
     * {@code FIXED_ORDER}.
     * 
     * @param dummies list of dummy nodes.
     * @param graph the layered graph.
     */
    private void restoreProperOrder(final List<LNode> dummies, final LGraph graph) {
        if (dummies.isEmpty()) {
            return;
        }
        
        // Turn the list into an array of dummy nodes and sort that by their original x coordinate
        LNode[] dummyArray = LGraphUtil.toNodeArray(dummies);
        
        Arrays.sort(dummyArray, new Comparator<LNode>() {
            public int compare(final LNode a, final LNode b) {
                return Double.compare(a.getProperty(InternalProperties.PORT_RATIO_OR_POSITION),
                        b.getProperty(InternalProperties.PORT_RATIO_OR_POSITION));
            }
        });
        
        assignAscendingCoordinates(dummyArray, graph);
    }
    
    /**
     * Iterates over the given array of dummy nodes, making sure that their x coordinates
     * are strictly ascending. Dummy nodes whose coordinates are in violation of this rule
     * are moved to the right. Once this method is finished, the coordinates of the dummy
     * nodes reflect their order in the array.
     * 
     * @param dummies array of dummy nodes.
     * @param graph the layered graph.
     */
    private void assignAscendingCoordinates(final LNode[] dummies, final LGraph graph) {
        // Find the port-port spacing
        double spacing = graph.getProperty(LayeredOptions.SPACING_PORT_PORT);
        
        // Now, iterate over the array, remembering the last assigned position. If we find a
        // position that is less than or equal to the last position, assign a new position of
        // "lastPosition + edgeSpacing"
        double nextValidCoordinate = dummies[0].getPosition().x
                + dummies[0].getSize().x
                + dummies[0].getMargin().right
                + spacing;
        
        for (int index = 1; index < dummies.length; index++) {
            KVector currentPosition = dummies[index].getPosition();
            KVector currentSize = dummies[index].getSize();
            LMargin currentMargin = dummies[index].getMargin();
            
            // Ensure spacings are adhered to
            double delta = currentPosition.x - currentMargin.left - nextValidCoordinate;
            if (delta < 0) {
                currentPosition.x -= delta;
            }
            
            // Ensure the graph is large enough for this node
            KVector graphSize = graph.getSize();
            graphSize.x = Math.max(graphSize.x, currentPosition.x + currentSize.x);
            
            // Compute next valid coordinate
            nextValidCoordinate = currentPosition.x + currentSize.x + currentMargin.right + spacing;
        }
    }
    
    
    ///////////////////////////////////////////////////////////////////////////////
    // STEP 3: EDGE ROUTING
    
    /**
     * Routes nothern and southern hierarchical port edges and ajusts the graph's height and
     * offsets accordingly.
     * 
     * @param monitor the progress monitor we're using.
     * @param layeredGraph the layered graph.
     * @param northSouthDummies the collection of restored northern and southern port dummies.
     */
    private void routeEdges(final IElkProgressMonitor monitor, final LGraph layeredGraph,
            final Iterable<LNode> northSouthDummies) {
        
        // Prepare south and target layers for northern and southern routing
        Set<LNode> northernSourceLayer = new LinkedHashSet<>();
        Set<LNode> northernTargetLayer = new LinkedHashSet<>();
        Set<LNode> southernSourceLayer = new LinkedHashSet<>();
        Set<LNode> southernTargetLayer = new LinkedHashSet<>();
        
        // Find some routing parameters
        double nodeSpacing = layeredGraph.getProperty(LayeredOptions.SPACING_NODE_NODE).doubleValue();
        double edgeSpacing = layeredGraph.getProperty(LayeredOptions.SPACING_EDGE_EDGE).doubleValue();
        
        // Assemble the northern and southern hierarchical port dummies and the nodes they are
        // connected to
        for (LNode hierarchicalPortDummy : northSouthDummies) {
            PortSide portSide = hierarchicalPortDummy.getProperty(InternalProperties.EXT_PORT_SIDE);
            
            if (portSide == PortSide.NORTH) {
                northernTargetLayer.add(hierarchicalPortDummy);
                
                for (LEdge edge : hierarchicalPortDummy.getIncomingEdges()) {
                    northernSourceLayer.add(edge.getSource().getNode());
                }
            } else if (portSide == PortSide.SOUTH) {
                southernTargetLayer.add(hierarchicalPortDummy);
                
                for (LEdge edge : hierarchicalPortDummy.getIncomingEdges()) {
                    southernSourceLayer.add(edge.getSource().getNode());
                }
            }
        }
        
        // Northern routing
        if (!northernSourceLayer.isEmpty()) {
            // Route the edges using a south-to-north orthogonal edge router
            OrthogonalRoutingGenerator routingGenerator = new OrthogonalRoutingGenerator(
                    RoutingDirection.SOUTH_TO_NORTH, edgeSpacing, "extnorth");
            
            int slots = routingGenerator.routeEdges(
                    monitor,
                    layeredGraph,
                    northernSourceLayer,
                    0,
                    northernTargetLayer,
                    -nodeSpacing - layeredGraph.getOffset().y);
            
            // If anything was routed, adjust the graph's offset and height
            if (slots > 0) {
                northernExtPortEdgeRoutingHeight = nodeSpacing + (slots - 1) * edgeSpacing;

                layeredGraph.getOffset().y += northernExtPortEdgeRoutingHeight;
                layeredGraph.getSize().y += northernExtPortEdgeRoutingHeight;
            }
        }
        
        // Southern routing
        if (!southernSourceLayer.isEmpty()) {
            // Route the edges using a north-to-south orthogonal edge router
            OrthogonalRoutingGenerator routingGenerator = new OrthogonalRoutingGenerator(
                    RoutingDirection.NORTH_TO_SOUTH, edgeSpacing, "extsouth");
            
            int slots = routingGenerator.routeEdges(
                    monitor,
                    layeredGraph,
                    southernSourceLayer,
                    0,
                    southernTargetLayer,
                    layeredGraph.getSize().y + nodeSpacing - layeredGraph.getOffset().y);
            
            // Adjust graph height.
            if (slots > 0) {
                layeredGraph.getSize().y += nodeSpacing + (slots - 1) * edgeSpacing;
            }
        }
    }
    
    
    ///////////////////////////////////////////////////////////////////////////////
    // STEP 4: REMOVE TEMPORARY DUMMIES
    
    /**
     * Removes the temporary hierarchical port dummies, reconnecting their incoming and outgoing
     * edges to the original dummies and setting the appropriate bend points.
     * 
     * @param layeredGraph the layered graph.
     */
    private void removeTemporaryNorthSouthDummies(final LGraph layeredGraph) {
        List<LNode> nodesToRemove = new ArrayList<>();
        
        // Iterate through all layers
        for (Layer layer : layeredGraph) {
            for (LNode node : layer) {
                if (node.getType() != NodeType.EXTERNAL_PORT) {
                    // We're only looking for hierarchical port dummies
                    continue;
                }
                
                if (!node.hasProperty(InternalProperties.EXT_PORT_REPLACED_DUMMY)) {
                    // We're only looking for temporary north / south dummies
                    continue;
                }
                
                // There must be a port where all edges come in, another port where edges go out, and
                // a port with an edge connecting node and origin (that one was added previously by
                // this processor)
                LPort nodeInPort = null;
                LPort nodeOutPort = null;
                LPort nodeOriginPort = null;
                
                for (LPort port : node.getPorts()) {
                    switch (port.getSide()) {
                    case WEST:
                        nodeInPort = port;
                        break;
                    
                    case EAST:
                        nodeOutPort = port;
                        break;
                        
                    default:
                        nodeOriginPort = port;
                    }
                }
                
                // Find the edge connecting this dummy to the original external port dummy that we
                // restored just a while ago
                LEdge nodeToOriginEdge = nodeOriginPort.getOutgoingEdges().get(0);
                
                // Compute bend points for incoming edges
                KVectorChain incomingEdgeBendPoints = new KVectorChain(nodeToOriginEdge.getBendPoints());
                
                KVector firstBendPoint = new KVector(nodeOriginPort.getPosition());
                firstBendPoint.add(node.getPosition());
                incomingEdgeBendPoints.add(0, firstBendPoint);
                
                // Compute bend points for outgoing edges
                KVectorChain outgoingEdgeBendPoints = KVectorChain.reverse(
                        nodeToOriginEdge.getBendPoints());
                
                KVector lastBendPoint = new KVector(nodeOriginPort.getPosition());
                lastBendPoint.add(node.getPosition());
                outgoingEdgeBendPoints.add(lastBendPoint);
                
                // Retrieve the original hierarchical port dummy
                LNode replacedDummy = (LNode) node.getProperty(
                        InternalProperties.EXT_PORT_REPLACED_DUMMY);
                LPort replacedDummyPort = replacedDummy.getPorts().get(0);
                
                // Reroute all the input port's edges
                LEdge[] edges = nodeInPort.getIncomingEdges().toArray(new LEdge[0]);
                
                for (LEdge edge : edges) {
                    edge.setTarget(replacedDummyPort);
                    edge.getBendPoints().addAllAsCopies(
                            edge.getBendPoints().size(), incomingEdgeBendPoints);
                }
                
                // Reroute all the output port's edges
                edges = LGraphUtil.toEdgeArray(nodeOutPort.getOutgoingEdges());
                
                for (LEdge edge : edges) {
                    edge.setSource(replacedDummyPort);
                    edge.getBendPoints().addAllAsCopies(0, outgoingEdgeBendPoints);
                }
                
                // Remove connection between node and original hierarchical port dummy
                nodeToOriginEdge.setSource(null);
                nodeToOriginEdge.setTarget(null);
                
                // Remember the temporary node for removal
                nodesToRemove.add(node);
            }
        }
        
        // Remove nodes
        for (LNode node : nodesToRemove) {
            node.setLayer(null);
        }
    }
    
    
    ///////////////////////////////////////////////////////////////////////////////
    // STEP 5: FIX DUMMY COORDINATES
    
    /**
     * Fixes all hierarchical port dummy coordinates. For east / west external port dummies, this
     * means setting the x coordinate appropriately, and, in case of {@code FIXED_RATIO},
     * checking that the ratio is respected. For north / south hierarchical port dummies, this
     * means setting the y coordinate appropriately.
     * 
     * @param layeredGraph the layered graph.
     */
    private void fixCoordinates(final LGraph layeredGraph) {
        PortConstraints constraints = layeredGraph.getProperty(LayeredOptions.PORT_CONSTRAINTS);
        
        // East port dummies are in the first layer; all other dummies are in the last layer
        List<Layer> layers = layeredGraph.getLayers();
        fixCoordinates(
                layers.get(0),
                constraints,
                layeredGraph);
        fixCoordinates(
                layers.get(layers.size() - 1),
                constraints,
                layeredGraph);
    }
    
    /**
     * Fixes the coordinates of the nodes in the given layer.
     * 
     * @param layer the layer.
     * @param constraints external port constraints.
     * @param graph the graph.
     */
    private void fixCoordinates(final Layer layer, final PortConstraints constraints,
            final LGraph graph) {
        
        // Get some geometric values from the graph
        LPadding padding = graph.getPadding();
        KVector offset = graph.getOffset();
        KVector graphActualSize = graph.getActualSize();
        
        double newActualGraphHeight = graphActualSize.y;
        
        // During the first iteration, EAST and WEST dummy nodes are fixed. This may change the height
        // of the graph, so we're setting y coordinates of NORTH and SOUTH dummies in a second iteration
        for (LNode node : layer) {
            if (node.getType() != NodeType.EXTERNAL_PORT) {
                // We're only looking for hierarchical port dummies
                continue;
            }
            
            PortSide extPortSide = node.getProperty(InternalProperties.EXT_PORT_SIDE);
            KVector extPortSize = node.getProperty(InternalProperties.EXT_PORT_SIZE);
            KVector nodePosition = node.getPosition();
            
            // Set x coordinate
            switch (extPortSide) {
            case EAST:
                nodePosition.x = graph.getSize().x + padding.right - offset.x;
                break;
            
            case WEST:
                nodePosition.x = -offset.x - padding.left;
                break;
            }
            
            // Set y coordinate
            double requiredActualGraphHeight = 0.0;
            
            switch (extPortSide) {
            case EAST:
            case WEST:
                if (constraints == PortConstraints.FIXED_RATIO) {
                    double ratio = node.getProperty(InternalProperties.PORT_RATIO_OR_POSITION);
                    nodePosition.y = graphActualSize.y * ratio
                            - node.getProperty(LayeredOptions.PORT_ANCHOR).y;
                    requiredActualGraphHeight = nodePosition.y + extPortSize.y;
                    node.borderToContentAreaCoordinates(false, true);
                } else if (constraints == PortConstraints.FIXED_POS) {
                    nodePosition.y = node.getProperty(InternalProperties.PORT_RATIO_OR_POSITION)
                            - node.getProperty(LayeredOptions.PORT_ANCHOR).y;
                    requiredActualGraphHeight = nodePosition.y + extPortSize.y;
                    node.borderToContentAreaCoordinates(false, true);
                }
                break;
            }
            
            newActualGraphHeight = Math.max(newActualGraphHeight, requiredActualGraphHeight);
        }
        
        // Make the graph larger, if necessary
        graph.getSize().y += newActualGraphHeight - graphActualSize.y;
        
        // Iterate over NORTH and SOUTH dummies now that the graph's height is fixed
        for (LNode node : layer) {
            if (node.getType() != NodeType.EXTERNAL_PORT) {
                // We're only looking for hierarchical port dummies
                continue;
            }
            
            PortSide extPortSide = node.getProperty(InternalProperties.EXT_PORT_SIDE);
            KVector nodePosition = node.getPosition();
            
            // Set y coordinate
            switch (extPortSide) {
            case NORTH:
                nodePosition.y = -offset.y - padding.top;
                break;
            
            case SOUTH:
                nodePosition.y = graph.getSize().y + padding.bottom - offset.y;
                break;
            }
        }
    }
    
    
    ///////////////////////////////////////////////////////////////////////////////
    // STEP 6: SLANTED EDGE SEGMENT CORRECTION
    
    /**
     * Goes over the eastern and western hierarchical dummy nodes and checks whether their
     * incident edges have slanted segments.
     * 
     * @param layeredGraph the layered graph.
     */
    private void correctSlantedEdgeSegments(final LGraph layeredGraph) {
        // East port dummies are in the first layer; all other dummies are in the last layer
        List<Layer> layers = layeredGraph.getLayers();
        correctSlantedEdgeSegments(layers.get(0));
        correctSlantedEdgeSegments(layers.get(layers.size() - 1));
    }
    
    /**
     * Goes over the eastern and western hierarchical dummy nodes in the given layer and checks
     * whether their incident edges have slanted segments. Note that only the first and last
     * segment needs to be checked.
     * 
     * @param layer the layer.
     */
    private void correctSlantedEdgeSegments(final Layer layer) {
        for (LNode node : layer) {
            if (node.getType() != NodeType.EXTERNAL_PORT) {
                // We're only looking for hierarchical port dummies
                continue;
            }
            
            PortSide extPortSide = node.getProperty(InternalProperties.EXT_PORT_SIDE);
            
            if (extPortSide == PortSide.EAST || extPortSide == PortSide.WEST) {
                for (LEdge edge : node.getConnectedEdges()) {
                    KVectorChain bendPoints = edge.getBendPoints();
                    
                    if (bendPoints.isEmpty()) {
                        // TODO: The edge has no bend points yet, but may still be slanted. Handle that!
                        continue;
                    }
                    
                    // Correct a slanted segment connected to the source port if it belongs to our node
                    LPort sourcePort = edge.getSource();
                    
                    if (sourcePort.getNode() == node) {
                        KVector firstBendPoint = bendPoints.getFirst();
                        firstBendPoint.y = sourcePort.getAbsoluteAnchor().y;
                    }
                    
                    // Correct a slanted segment connected to the target port if it belongs to our node
                    LPort targetPort = edge.getTarget();
                    
                    if (targetPort.getNode() == node) {
                        KVector lastBendPoint = bendPoints.getLast();
                        lastBendPoint.y = targetPort.getAbsoluteAnchor().y;
                    }
                }
            }
        }
    }
    
}
