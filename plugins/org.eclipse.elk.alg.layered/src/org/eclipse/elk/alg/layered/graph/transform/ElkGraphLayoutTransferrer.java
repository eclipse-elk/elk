/*******************************************************************************
 * Copyright (c) 2015, 2020 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.alg.layered.graph.transform;

import java.util.EnumSet;
import java.util.List;
import java.util.Set;
import java.util.ArrayList;

import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LGraphUtil;
import org.eclipse.elk.alg.layered.graph.LLabel;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LPadding;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.intermediate.LabelDummySwitcher;
import org.eclipse.elk.alg.layered.options.GraphProperties;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.layered.options.NodeFlexibility;
import org.eclipse.elk.alg.layered.options.NodePlacementStrategy;
import org.eclipse.elk.core.math.ElkPadding;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.math.KVectorChain;
import org.eclipse.elk.core.options.EdgeRouting;
import org.eclipse.elk.core.options.PortConstraints;
import org.eclipse.elk.core.options.PortLabelPlacement;
import org.eclipse.elk.core.options.SizeConstraint;
import org.eclipse.elk.core.options.SizeOptions;
import org.eclipse.elk.core.util.ElkUtil;
import org.eclipse.elk.graph.ElkEdge;
import org.eclipse.elk.graph.ElkEdgeSection;
import org.eclipse.elk.graph.ElkLabel;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.ElkPort;
import org.eclipse.elk.graph.util.ElkGraphUtil;


/**
 * Implements the graph layout application aspect of {@link ElkGraphTransformer}.
 * 
 * @author cds, sdo
 */
public class ElkGraphLayoutTransferrer {

    /**
     * Applies the layout information contained in the given LGraph to the ElkGraph elements it was
     * created from. All source ElkGraph elements are expected to be accessible through their LGraph
     * counterparts through the {@link InternalProperties#ORIGIN} property.
     * 
     * @param lgraph the LGraph whose layout information to apply.
     */
    public void applyLayout(final LGraph lgraph) {
        Object graphOrigin = lgraph.getProperty(InternalProperties.ORIGIN);
        if (!(graphOrigin instanceof ElkNode)) {
            return;
        }
        
        // The ElkNode that represents this graph in the original ElkGraph
        ElkNode parentElkNode = (ElkNode) graphOrigin;
        
        // The LNode that represents this graph in the upper hierarchy level, if any
        LNode parentLNode = (LNode) lgraph.getParentNode();
        
        // Get the offset to be added to all coordinates
        KVector offset = new KVector(lgraph.getOffset());
        
        // Adjust offset (and with it the positions) by the requested padding
        LPadding lPadding = lgraph.getPadding();
        offset.x += lPadding.left;
        offset.y += lPadding.top;
        
        // Set node padding, if it was computed during layout
        final EnumSet<SizeOptions> sizeOptions = parentElkNode.getProperty(LayeredOptions.NODE_SIZE_OPTIONS);
        if (sizeOptions.contains(SizeOptions.COMPUTE_PADDING)) {
            ElkPadding padding = parentElkNode.getProperty(LayeredOptions.PADDING);
            padding.setBottom(lPadding.bottom);
            padding.setTop(lPadding.top);
            padding.setLeft(lPadding.left);
            padding.setRight(lPadding.right);
        }

        // Along the way, we collect the list of edges to be processed later
        List<LEdge> edgeList = new ArrayList<>();

        // Process the nodes
        for (LNode lnode : lgraph.getLayerlessNodes()) {
            if (representsNode(lnode)) {
                applyNodeLayout(lnode, offset);
                
            } else if (representsExternalPort(lnode) && parentLNode == null) {
                // We have an external port here on the top-most hierarchy level of the current (possibly
                // hierarchical) layout run; set its position
                ElkPort elkport = (ElkPort) lnode.getProperty(InternalProperties.ORIGIN);
                KVector portPosition = LGraphUtil.getExternalPortPosition(
                        lgraph, lnode, elkport.getWidth(), elkport.getHeight());
                elkport.setLocation(portPosition.x, portPosition.y);
            }

            // Collect edges, except if they go into a nested subgraph (those edges need to be
            // processed during one of the recursive calls so that any additional offsets are applied
            // correctly)
            for (LPort port : lnode.getPorts()) {
                port.getOutgoingEdges().stream()
                    .filter(edge -> !LGraphUtil.isDescendant(edge.getTarget().getNode(), lnode))
                    .forEach(edge -> edgeList.add(edge));
            }
        }

        // Collect edges that go from the current graph's representing LNode down into its descendants
        if (parentLNode != null) {
            for (LPort port : parentLNode.getPorts()) {
                port.getOutgoingEdges().stream()
                        .filter(edge -> LGraphUtil.isDescendant(edge.getTarget().getNode(), parentLNode))
                        .forEach(edge -> edgeList.add(edge));
            }
        }
        
        // Iterate through all edges
        EdgeRouting routing = parentElkNode.getProperty(LayeredOptions.EDGE_ROUTING);
        for (LEdge ledge : edgeList) {
            applyEdgeLayout(ledge, routing, offset, lPadding);
        }

        // Setup the parent node
        applyParentNodeLayout(lgraph);
        
        // Process nested subgraphs
        for (LNode lnode : lgraph.getLayerlessNodes()) {
            LGraph nestedGraph = lnode.getNestedGraph();
            if (nestedGraph != null) {
                applyLayout(nestedGraph);
            }
        }
    }

    /**
     * Applies layout information computed for the given node.
     * 
     * @param lnode
     *            the node that has the layout information.
     * @param offset
     *            offset to add to coordinates.
     */
    private void applyNodeLayout(final LNode lnode, final KVector offset) {
        final ElkNode elknode = (ElkNode) lnode.getProperty(InternalProperties.ORIGIN);
        
        
        // Apply the nodeID and layerId that were set on the LGraph on the ElkGraph
        final int nodeID = lnode.getProperty(LayeredOptions.CROSSING_MINIMIZATION_POSITION_ID);
        final int layerID = lnode.getProperty(LayeredOptions.LAYERING_LAYER_ID);
        elknode.setProperty(LayeredOptions.CROSSING_MINIMIZATION_POSITION_ID, nodeID);
        elknode.setProperty(LayeredOptions.LAYERING_LAYER_ID, layerID);
        
        // Set the node position
        elknode.setX(lnode.getPosition().x + offset.x);
        elknode.setY(lnode.getPosition().y + offset.y);
        
        // Set the node size, if necessary
        if (!elknode.getProperty(LayeredOptions.NODE_SIZE_CONSTRAINTS).isEmpty()
                || lnode.getNestedGraph() != null
                || (lnode.getGraph().getProperty(LayeredOptions.NODE_PLACEMENT_STRATEGY) 
                        == NodePlacementStrategy.NETWORK_SIMPLEX 
                    && NodeFlexibility.getNodeFlexibility(lnode).isFlexibleSizeWhereSpacePermits())) {
            
            elknode.setWidth(lnode.getSize().x);
            elknode.setHeight(lnode.getSize().y);
        }

        // Set port positions
        for (LPort lport : lnode.getPorts()) {
            Object origin = lport.getProperty(InternalProperties.ORIGIN);
            if (origin instanceof ElkPort) {
                ElkPort elkport = (ElkPort) origin;
                elkport.setLocation(lport.getPosition().x, lport.getPosition().y);
                elkport.setProperty(LayeredOptions.PORT_SIDE, lport.getSide());
            }
        }
        
        // Set node label positions, if they were not fixed
        // (that is at least one of the node or the label has a node label placement set)
        final boolean nodeHasLabelPlacement = !lnode.getProperty(LayeredOptions.NODE_LABELS_PLACEMENT).isEmpty();
        
        for (LLabel llabel : lnode.getLabels()) {
            if (nodeHasLabelPlacement || !llabel.getProperty(LayeredOptions.NODE_LABELS_PLACEMENT).isEmpty()) {
                ElkLabel elklabel = (ElkLabel) llabel.getProperty(InternalProperties.ORIGIN);
                elklabel.setDimensions(llabel.getSize().x, llabel.getSize().y);
                elklabel.setLocation(llabel.getPosition().x, llabel.getPosition().y);
            }
        }
        
        // Set port label positions, if they were not fixed
        if (!PortLabelPlacement.isFixed(lnode.getProperty(LayeredOptions.PORT_LABELS_PLACEMENT))) {
            for (LPort lport : lnode.getPorts()) {
                for (LLabel llabel : lport.getLabels()) {
                    ElkLabel elklabel = (ElkLabel) llabel.getProperty(InternalProperties.ORIGIN);
                    elklabel.setWidth(llabel.getSize().x);
                    elklabel.setHeight(llabel.getSize().y);
                    elklabel.setLocation(llabel.getPosition().x, llabel.getPosition().y);
                }
            }
        }
    }
    
    /**
     * Applies layout information computed for the given edge.
     * 
     * @param ledge
     *            the edge that has the layout information.
     * @param routing
     *            the kind of routing applied to edges.
     * @param offset
     *            offset to add to coordinates.
     * @param additionalPadding
     *            the additional insets that may have to be taken into account for hierarchical that go
     *            into the bowels of their source node. These are already included in the offset, but
     *            are required separately.
     */
    private void applyEdgeLayout(final LEdge ledge, final EdgeRouting routing, final KVector offset,
            final LPadding additionalPadding) {

        ElkEdge elkedge = (ElkEdge) ledge.getProperty(InternalProperties.ORIGIN);
        
        // Only the orthogonal edge routing algorithm supports self-loops. Thus, leave self-loops
        // untouched if another routing algorithm is selected.
        if (elkedge == null) {
            return;
        }
        
        KVectorChain bendPoints = ledge.getBendPoints();
        
        // The standard offset may need to be modified if the edge needs to end up in a coordinate system of
        // a graph in a higher hierarchy level
        KVector edgeOffset = new KVector(offset);
        edgeOffset.add(calculateHierarchicalOffset(ledge));

        // Adapt the offset value and add the source port position to the vector chain
        KVector sourcePoint;
        if (LGraphUtil.isDescendant(ledge.getTarget().getNode(), ledge.getSource().getNode())) {
            // The external port's anchor position, relative to the node's top left corner
            LPort sourcePort = ledge.getSource();
            sourcePoint = KVector.sum(sourcePort.getPosition(), sourcePort.getAnchor());
            
            // The source point will later have the passed offset added to it, which it doesn't actually
            // need, so we subtract it now (notice that the external port's position was already relative
            // to its node's top left corner, while adding the offset now would mean that it was relative
            // to the top left corner + its insets area)
            sourcePoint.sub(offset);
        } else {
            sourcePoint = ledge.getSource().getAbsoluteAnchor();
        }
        bendPoints.addFirst(sourcePoint);
        
        // Add the target port position to the vector chain, including additional offset
        KVector targetPoint = ledge.getTarget().getAbsoluteAnchor();
        if (ledge.getProperty(InternalProperties.TARGET_OFFSET) != null) {
            targetPoint.add(ledge.getProperty(InternalProperties.TARGET_OFFSET));
        }
        bendPoints.addLast(targetPoint);

        // Translate the bend points by the offset and apply the bend points
        bendPoints.offset(edgeOffset);
        
        // Give the edge a proper edge section to store routing information
        ElkEdgeSection elkedgeSection = ElkGraphUtil.firstEdgeSection(elkedge, true, true);
        elkedgeSection.setIncomingShape(elkedge.getSources().get(0));
        elkedgeSection.setOutgoingShape(elkedge.getTargets().get(0));
        ElkUtil.applyVectorChain(bendPoints, elkedgeSection);

        // Apply layout to labels
        for (LLabel llabel : ledge.getLabels()) {
            ElkLabel elklabel = (ElkLabel) llabel.getProperty(InternalProperties.ORIGIN);
            elklabel.setWidth(llabel.getSize().x);
            elklabel.setHeight(llabel.getSize().y);
            elklabel.setLocation(llabel.getPosition().x + edgeOffset.x,
                                 llabel.getPosition().y + edgeOffset.y);
            
            elklabel.setProperty(
                    LabelDummySwitcher.INCLUDE_LABEL,
                    llabel.getProperty(LabelDummySwitcher.INCLUDE_LABEL));
        }
        
        // Copy junction points
        KVectorChain junctionPoints = ledge.getProperty(LayeredOptions.JUNCTION_POINTS);
        if (junctionPoints != null) {
            junctionPoints.offset(edgeOffset);
            elkedge.setProperty(LayeredOptions.JUNCTION_POINTS, junctionPoints);
        } else {
            elkedge.setProperty(LayeredOptions.JUNCTION_POINTS, null);
        }

        // Mark the edge with information about its routing
        if (routing == EdgeRouting.SPLINES) {
            // SPLINES means that bend points shall be interpreted as control points for splines
            elkedge.setProperty(LayeredOptions.EDGE_ROUTING, EdgeRouting.SPLINES);
        } else {
            // null means that bend points shall be interpreted as bend points
            elkedge.setProperty(LayeredOptions.EDGE_ROUTING, null);
        }
    }
    
    /** A zero vector to avoid having to instantiate new empty vectors for each edge. */
    private static final KVector ZERO_OFFSET = new KVector();
    
    /**
     * If the coordinates of an edge must be relative to a different node than they are in the algorithm, this
     * method returns the correct offset to translate from the algorithm's coordinate system to the necessary
     * target coordinate system.
     * 
     * @return the offset vector, which may simply be zero (but not {@code null}). Must not be modified.
     */
    private KVector calculateHierarchicalOffset(final LEdge ledge) {
        LGraph targetCoordinateSystem = ledge.getProperty(InternalProperties.COORDINATE_SYSTEM_ORIGIN);
        if (targetCoordinateSystem != null) {
            KVector result = new KVector();
            
            // Edges this method is called on are always in the coordinate system of their source
            LGraph currentGraph = ledge.getSource().getNode().getGraph();
            
            while (currentGraph != targetCoordinateSystem) {
                // The current graph should always have an upper level if we have not reached the target graph yet;
                LNode representingNode = currentGraph.getParentNode();
                currentGraph = representingNode.getGraph();
                
                result.add(representingNode.getPosition())
                      .add(currentGraph.getOffset())
                      .add(currentGraph.getPadding().left, currentGraph.getPadding().top);
            }
            
            return result;
        }
        
        // No coordinate system conversion is required
        return ZERO_OFFSET;
    }

    /**
     * Applies layout information computed for the given graph.
     * 
     * @param lgraph
     *            the edge that has the layout information.
     */
    private void applyParentNodeLayout(final LGraph lgraph) {
        ElkNode elknode = (ElkNode) lgraph.getProperty(InternalProperties.ORIGIN);
        
        boolean sizeConstraintsIncludedPortLabels =
                elknode.getProperty(LayeredOptions.NODE_SIZE_CONSTRAINTS).contains(SizeConstraint.PORT_LABELS);

        if (lgraph.getParentNode() == null) {
            Set<GraphProperties> graphProps = lgraph.getProperty(InternalProperties.GRAPH_PROPERTIES);
            KVector actualGraphSize = lgraph.getActualSize();
            
            if (graphProps.contains(GraphProperties.EXTERNAL_PORTS)) {
                // Ports have positions assigned, the graph's size is final
                elknode.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_POS);
                ElkUtil.resizeNode(
                        elknode,
                        actualGraphSize.x,
                        actualGraphSize.y,
                        false,
                        true);
            } else {
                // Ports have not been positioned yet - leave this for next layouter
                if (!elknode.getProperty(LayeredOptions.NODE_SIZE_FIXED_GRAPH_SIZE)) {
                    ElkUtil.resizeNode(
                            elknode,
                            actualGraphSize.x,
                            actualGraphSize.y,
                            true,
                            true);
                }
            }
        }
        
        // Set the size constraints. Thing is, if we always set the size constraints to fixed, we might remove a
        // PORT_LABELS constraints that was formerly present. This in turn will cause ports to be placed without
        // consideration for their labels.
        if (sizeConstraintsIncludedPortLabels) {
            elknode.setProperty(LayeredOptions.NODE_SIZE_CONSTRAINTS, EnumSet.of(SizeConstraint.PORT_LABELS));
        } else {
            elknode.setProperty(LayeredOptions.NODE_SIZE_CONSTRAINTS, SizeConstraint.fixed());
        }
    }
    
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Utility Methods
    
    /**
     * Checks if the given node represents a node in the original ELK graph as well.
     */
    private static boolean representsNode(final LNode lnode) {
        return lnode.getProperty(InternalProperties.ORIGIN) instanceof ElkNode;
    }
    
    /**
     * Checks if the given node represents an external port in the original ELK graph.
     */
    private static boolean representsExternalPort(final LNode lnode) {
        return lnode.getProperty(InternalProperties.ORIGIN) instanceof ElkPort;
    }
    

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Utility Methods for debugging
    
    /**
     * Debug method used to create an ElkGraph from an LGraph to log via the progress monitor.
     * This is done by executing the same steps as in {@link #applyLayout(LGraph)}, but without mutating.
     * 
     * @param lgraph The LGraph
     * @return The ElkNode to log.
     */
    public static ElkNode applyLayoutNonMutating(LGraph lgraph) {
        // The ElkNode that represents this graph in the original ElkGraph
        ElkNode parentElkNode = ElkGraphUtil.createGraph();
        
        // Apply the layout non-mutating
        return applyLayoutNonMutating(lgraph, parentElkNode);
    }
    
    /**
     * Recursive method used to create an ElkGraph from an LGraph to log via the progress monitor.
     * This is done by executing the same steps as in {@link #applyLayout(LGraph)}, but without mutating.
     * 
     * @param lgraph The LGraph
     * @param parentElkNode The ElkNode that represents this graph in the original ElkGraph.
     * @return The ElkNode to log.
     */
    private static ElkNode applyLayoutNonMutating(LGraph lgraph, ElkNode parentElkNode) {
        
        // The LNode that represents this graph in the upper hierarchy level, if any
        LNode parentLNode = (LNode) lgraph.getParentNode();

        
        // Get the offset to be added to all coordinates
        KVector offset = new KVector(lgraph.getOffset().clone());
        
        // Adjust offset (and with it the positions) by the requested padding
        LPadding lPadding = (LPadding) lgraph.getPadding().clone();
        offset.x += lPadding.left;
        offset.y += lPadding.top;
        
        // Set node padding, if it was computed during layout
        final EnumSet<SizeOptions> sizeOptions = parentElkNode.getProperty(LayeredOptions.NODE_SIZE_OPTIONS);
        if (sizeOptions.contains(SizeOptions.COMPUTE_PADDING)) {
            ElkPadding padding = parentElkNode.getProperty(LayeredOptions.PADDING);
            padding.setBottom(lPadding.bottom);
            padding.setTop(lPadding.top);
            padding.setLeft(lPadding.left);
            padding.setRight(lPadding.right);
        }

        // Along the way, we collect the list of edges to be processed later
        List<LEdge> edgeList = new ArrayList<>();

        // Process the nodes
        for (LNode lnode : collectNodes(lgraph)) {
            ElkNode elkNode = null;
            if (representsNode(lnode)) {
                // Create a new node
                elkNode = ElkGraphUtil.createNode(parentElkNode);
                applyNodeLayout(lnode, elkNode, offset);               
                
            } else if (representsExternalPort(lnode) && parentLNode == null) {
                // We have an external port here on the top-most hierarchy level of the current (possibly
                // hierarchical) layout run; set its position
                // Create a new port
                ElkPort elkport = ElkGraphUtil.createPort(parentElkNode);
                KVector portPosition = LGraphUtil.getExternalPortPosition(
                        lgraph, lnode, elkport.getWidth(), elkport.getHeight());
                elkport.setLocation(portPosition.x, portPosition.y);
            }

            // Collect edges, except if they go into a nested subgraph (those edges need to be
            // processed during one of the recursive calls so that any additional offsets are applied
            // correctly)
            for (LPort port : lnode.getPorts()) {
                port.getOutgoingEdges().stream()
                    .filter(edge -> !LGraphUtil.isDescendant(edge.getTarget().getNode(), lnode))
                    .forEach(edge -> edgeList.add(edge));
            }
            
            // Transform child nodes.
            LGraph nestedGraph = lnode.getNestedGraph();
            if (nestedGraph != null) {
                applyLayoutNonMutating(nestedGraph, elkNode);
                parentElkNode.getChildren().add(elkNode);
            }
        }

        // Collect edges that go from the current graph's representing LNode down into its descendants
        if (parentLNode != null) {
            for (LPort port : parentLNode.getPorts()) {
                port.getOutgoingEdges().stream()
                        .filter(edge -> LGraphUtil.isDescendant(edge.getTarget().getNode(), parentLNode))
                        .forEach(edge -> edgeList.add(edge));
            }
        }
        
        // Iterate through all edges
        EdgeRouting routing = parentElkNode.getProperty(LayeredOptions.EDGE_ROUTING);
        for (LEdge ledge : edgeList) {
            // Create a new edge.
            // This edge will not be connected to any source or target, as it is only used for logging.
            ElkEdge elkedge = ElkGraphUtil.createEdge(parentElkNode);
            applyEdgeLayout(ledge, elkedge, routing, offset, lPadding);
        }

        // Setup the parent node
        applyParentNodeLayout(lgraph, parentElkNode);
        return parentElkNode;
        
    }
    
    /**
     * Collects all nodes from the given LGraph that represent a node in the original ELK graph.
     * 
     * @param lgraph The LGraph to collect nodes from.
     * @return A list of all nodes that represent a node in the original ELK graph.
     */
    private static List<LNode> collectNodes(LGraph lgraph) {
        List<LNode> nodes = new ArrayList<>();
        for (Layer layer : lgraph.getLayers()) {
            for (LNode lNode : layer.getNodes()) {
                if (representsNode(lNode)) {
                    nodes.add(lNode);
                } else if (representsExternalPort(lNode)) {
                    // External ports are not nodes, so we do not add them
                } else {
                    // This is a node that does not represent a node, e.g. a dummy node
                    // or a node that is not yet placed.
                    continue;
                }
            }
        }
        for (LNode lNode : lgraph.getLayerlessNodes()) {
            if (representsNode(lNode)) {
                nodes.add(lNode);
            } else if (representsExternalPort(lNode)) {
                // External ports are not nodes, so we do not add them
            } else {
                // This is a node that does not represent a node, e.g. a dummy node
                // or a node that is not yet placed.
                continue;
            }
        }
        return nodes;
    }
    
    /**
     * Applies layout information computed for the given node.
     * 
     * @param lnode
     *            the node that has the layout information.
     * @param elknode
     *            the ElkNode to apply the layout information to.
     * @param offset
     *            offset to add to coordinates.
     */
    private static void applyNodeLayout(final LNode lnode, final ElkNode elknode, final KVector offset) {        
        
        // Apply the nodeID and layerId that were set on the LGraph on the ElkGraph
        final int nodeID = lnode.getProperty(LayeredOptions.CROSSING_MINIMIZATION_POSITION_ID);
        final int layerID = lnode.getProperty(LayeredOptions.LAYERING_LAYER_ID);
        elknode.setProperty(LayeredOptions.CROSSING_MINIMIZATION_POSITION_ID, nodeID);
        elknode.setProperty(LayeredOptions.LAYERING_LAYER_ID, layerID);
        
        // Set the node position
        elknode.setX(lnode.getPosition().x + offset.x);
        elknode.setY(lnode.getPosition().y + offset.y);
        
        // Set the node size, if necessary
        if (!elknode.getProperty(LayeredOptions.NODE_SIZE_CONSTRAINTS).isEmpty()
                || lnode.getNestedGraph() != null
                || (lnode.getGraph().getProperty(LayeredOptions.NODE_PLACEMENT_STRATEGY) 
                        == NodePlacementStrategy.NETWORK_SIMPLEX 
                    && NodeFlexibility.getNodeFlexibility(lnode).isFlexibleSizeWhereSpacePermits())) {
            
            elknode.setWidth(lnode.getSize().x);
            elknode.setHeight(lnode.getSize().y);
        }

        // Set port positions
        for (LPort lport : lnode.getPorts()) {
            // Create a new port
            ElkPort elkport = ElkGraphUtil.createPort(elknode);
            elkport.setLocation(lport.getPosition().x, lport.getPosition().y);
            elkport.setProperty(LayeredOptions.PORT_SIDE, lport.getSide());
            
            // Set port label positions, if they were not fixed
            if (!PortLabelPlacement.isFixed(lnode.getProperty(LayeredOptions.PORT_LABELS_PLACEMENT))) {
                for (LLabel llabel : lport.getLabels()) {
                    // Create a new label
                    ElkLabel elklabel = ElkGraphUtil.createLabel(elkport);
                    elklabel.setWidth(llabel.getSize().x);
                    elklabel.setHeight(llabel.getSize().y);
                    elklabel.setLocation(llabel.getPosition().x, llabel.getPosition().y);
                }
            }
        }
        
        // Set node label positions, if they were not fixed
        // (that is at least one of the node or the label has a node label placement set)
        final boolean nodeHasLabelPlacement = !lnode.getProperty(LayeredOptions.NODE_LABELS_PLACEMENT).isEmpty();
        
        for (LLabel llabel : lnode.getLabels()) {
            if (nodeHasLabelPlacement || !llabel.getProperty(LayeredOptions.NODE_LABELS_PLACEMENT).isEmpty()) {
                // Create a new label
                ElkLabel elklabel = ElkGraphUtil.createLabel(elknode);
                elklabel.setDimensions(llabel.getSize().x, llabel.getSize().y);
                elklabel.setLocation(llabel.getPosition().x, llabel.getPosition().y);
            }
        }
    }
    
    /**
     * Applies layout information computed for the given edge.
     * 
     * @param ledge
     *            the edge that has the layout information.
     * @param routing
     *            the kind of routing applied to edges.
     * @param offset
     *            offset to add to coordinates.
     */
    private static void applyEdgeLayout(final LEdge ledge, final ElkEdge elkedge, final EdgeRouting routing, final KVector offset,
            final LPadding additionalPadding) {
        
        // Only the orthogonal edge routing algorithm supports self-loops. Thus, leave self-loops
        // untouched if another routing algorithm is selected.
        if (elkedge == null) {
            return;
        }
        
        KVectorChain bendPoints = (KVectorChain) ledge.getBendPoints().clone();
        
        // The standard offset may need to be modified if the edge needs to end up in a coordinate system of
        // a graph in a higher hierarchy level
        KVector edgeOffset = new KVector(offset);
        edgeOffset = KVector.sum(edgeOffset, calculateHierarchicalOffsetNonMutating(ledge));

        // Adapt the offset value and add the source port position to the vector chain
        KVector sourcePoint;
        if (LGraphUtil.isDescendant(ledge.getTarget().getNode(), ledge.getSource().getNode())) {
            // The external port's anchor position, relative to the node's top left corner
            LPort sourcePort = ledge.getSource();
            sourcePoint = KVector.sum(sourcePort.getPosition(), sourcePort.getAnchor());
            
            // The source point will later have the passed offset added to it, which it doesn't actually
            // need, so we subtract it now (notice that the external port's position was already relative
            // to its node's top left corner, while adding the offset now would mean that it was relative
            // to the top left corner + its insets area)
            KVector.diff(sourcePoint, offset);
        } else {
            sourcePoint = ledge.getSource().getAbsoluteAnchor().clone();
        }
        bendPoints.addFirst(sourcePoint);
        
        // Add the target port position to the vector chain, including additional offset
        KVector targetPoint = ledge.getTarget().getAbsoluteAnchor().clone();
        if (ledge.getProperty(InternalProperties.TARGET_OFFSET) != null) {
            targetPoint.add(ledge.getProperty(InternalProperties.TARGET_OFFSET));
        }
        bendPoints.addLast(targetPoint);

        // Translate the bend points by the offset and apply the bend points
        // Do not mutate the bend points, but create a new one instead
        bendPoints = bendPoints.offsetNonMutating(edgeOffset);
        
        // Give the edge a proper edge section to store routing information
        // Source and target do not matter here.
        ElkEdgeSection elkedgeSection = ElkGraphUtil.firstEdgeSection(elkedge, true, true);
        ElkUtil.applyVectorChain(bendPoints, elkedgeSection);

        // Apply layout to labels
        for (LLabel llabel : ledge.getLabels()) {
            // Create a new label
            ElkLabel elklabel = ElkGraphUtil.createLabel(elkedge);
            elklabel.setWidth(llabel.getSize().x);
            elklabel.setHeight(llabel.getSize().y);
            elklabel.setLocation(llabel.getPosition().x + edgeOffset.x,
                                 llabel.getPosition().y + edgeOffset.y);
            
            elklabel.setProperty(
                    LabelDummySwitcher.INCLUDE_LABEL,
                    llabel.getProperty(LabelDummySwitcher.INCLUDE_LABEL));
        }
        
        // Copy junction points
        KVectorChain junctionPoints = (KVectorChain) ledge.getProperty(LayeredOptions.JUNCTION_POINTS).clone();
        if (junctionPoints != null) {
            junctionPoints = junctionPoints.offsetNonMutating(edgeOffset);
            elkedge.setProperty(LayeredOptions.JUNCTION_POINTS, junctionPoints);
        } else {
            elkedge.setProperty(LayeredOptions.JUNCTION_POINTS, null);
        }

        // Mark the edge with information about its routing
        if (routing == EdgeRouting.SPLINES) {
            // SPLINES means that bend points shall be interpreted as control points for splines
            elkedge.setProperty(LayeredOptions.EDGE_ROUTING, EdgeRouting.SPLINES);
        } else {
            // null means that bend points shall be interpreted as bend points
            elkedge.setProperty(LayeredOptions.EDGE_ROUTING, null);
        }
    }    
    /**
     * If the coordinates of an edge must be relative to a different node than they are in the algorithm, this
     * method returns the correct offset to translate from the algorithm's coordinate system to the necessary
     * target coordinate system.
     * 
     * @return the offset vector, which may simply be zero (but not {@code null}). Must not be modified.
     */
    private static KVector calculateHierarchicalOffsetNonMutating(final LEdge ledge) {
        LGraph targetCoordinateSystem = ledge.getProperty(InternalProperties.COORDINATE_SYSTEM_ORIGIN);
        if (targetCoordinateSystem != null) {
            KVector result = new KVector();
            
            // Edges this method is called on are always in the coordinate system of their source
            LGraph currentGraph = ledge.getSource().getNode().getGraph();
            
            while (currentGraph != targetCoordinateSystem) {
                // The current graph should always have an upper level if we have not reached the target graph yet;
                LNode representingNode = currentGraph.getParentNode();
                currentGraph = representingNode.getGraph();
                
                // Do not mutate the result, but create a new one instead
                result = KVector.sum(result, representingNode.getPosition(), currentGraph.getOffset(),
                        new KVector(currentGraph.getPadding().left, currentGraph.getPadding().top));
            }
            
            return result;
        }
        
        // No coordinate system conversion is required
        return ZERO_OFFSET;
    }
    
    /**
     * Applies layout information computed for the given graph without mutating the original graph.
     * 
     * @param lgraph
     *            the edge that has the layout information.
     */
    private static void applyParentNodeLayout(LGraph lgraph, ElkNode elknode) {
        boolean sizeConstraintsIncludedPortLabels =
                elknode.getProperty(LayeredOptions.NODE_SIZE_CONSTRAINTS).contains(SizeConstraint.PORT_LABELS);

        if (lgraph.getParentNode() == null) {
            Set<GraphProperties> graphProps = lgraph.getProperty(InternalProperties.GRAPH_PROPERTIES);
            KVector actualGraphSize = lgraph.getActualSize();
            
            if (graphProps.contains(GraphProperties.EXTERNAL_PORTS)) {
                // Ports have positions assigned, the graph's size is final
                elknode.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_POS);
                ElkUtil.resizeNode(
                        elknode,
                        actualGraphSize.x,
                        actualGraphSize.y,
                        false,
                        true);
            } else {
                // Ports have not been positioned yet - leave this for next layouter
                if (!elknode.getProperty(LayeredOptions.NODE_SIZE_FIXED_GRAPH_SIZE)) {
                    ElkUtil.resizeNode(
                            elknode,
                            actualGraphSize.x,
                            actualGraphSize.y,
                            true,
                            true);
                }
            }
        }
        
        // Set the size constraints. Thing is, if we always set the size constraints to fixed, we might remove a
        // PORT_LABELS constraints that was formerly present. This in turn will cause ports to be placed without
        // consideration for their labels.
        if (sizeConstraintsIncludedPortLabels) {
            elknode.setProperty(LayeredOptions.NODE_SIZE_CONSTRAINTS, EnumSet.of(SizeConstraint.PORT_LABELS));
        } else {
            elknode.setProperty(LayeredOptions.NODE_SIZE_CONSTRAINTS, SizeConstraint.fixed());
        }
    }

}
