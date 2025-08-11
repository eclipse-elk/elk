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
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Queue;
import java.util.Set;

import org.eclipse.elk.alg.common.nodespacing.NodeLabelAndSizeCalculator;
import org.eclipse.elk.alg.layered.components.ComponentOrderingStrategy;
import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LGraphElement;
import org.eclipse.elk.alg.layered.graph.LGraphUtil;
import org.eclipse.elk.alg.layered.graph.LLabel;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LPadding;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.options.CrossingMinimizationStrategy;
import org.eclipse.elk.alg.layered.options.CycleBreakingStrategy;
import org.eclipse.elk.alg.layered.options.GraphProperties;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.layered.options.LayeredSpacings;
import org.eclipse.elk.alg.layered.options.LayeringStrategy;
import org.eclipse.elk.alg.layered.options.NodePlacementStrategy;
import org.eclipse.elk.alg.layered.options.NodePromotionStrategy;
import org.eclipse.elk.alg.layered.options.OrderingStrategy;
import org.eclipse.elk.alg.layered.options.PortType;
import org.eclipse.elk.core.UnsupportedGraphException;
import org.eclipse.elk.core.labels.LabelManagementOptions;
import org.eclipse.elk.core.math.ElkPadding;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.math.KVectorChain;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.options.Direction;
import org.eclipse.elk.core.options.EdgeLabelPlacement;
import org.eclipse.elk.core.options.HierarchyHandling;
import org.eclipse.elk.core.options.PortConstraints;
import org.eclipse.elk.core.options.PortLabelPlacement;
import org.eclipse.elk.core.options.PortSide;
import org.eclipse.elk.core.options.SizeConstraint;
import org.eclipse.elk.core.util.ElkUtil;
import org.eclipse.elk.core.util.adapters.ElkGraphAdapters;
import org.eclipse.elk.core.util.adapters.GraphAdapters.GraphAdapter;
import org.eclipse.elk.core.util.adapters.GraphAdapters.NodeAdapter;
import org.eclipse.elk.graph.ElkConnectableShape;
import org.eclipse.elk.graph.ElkEdge;
import org.eclipse.elk.graph.ElkEdgeSection;
import org.eclipse.elk.graph.ElkGraphElement;
import org.eclipse.elk.graph.ElkLabel;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.ElkPort;
import org.eclipse.elk.graph.util.ElkGraphUtil;
import org.eclipse.emf.ecore.util.EcoreUtil;

import com.google.common.base.Strings;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

/**
 * Implements the graph import aspect of {@link ElkGraphTransformer}.
 */
class ElkGraphImporter {
    
    /** map between ElkGraph nodes / ports and the LGraph nodes / ports created for them. */
    private final Map<ElkGraphElement, LGraphElement> nodeAndPortMap = Maps.newHashMap();
    

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Import Entry Points

    /**
     * Imports the given graph.
     * 
     * @param elkgraph
     *            the graph to import.
     * @return the transformed graph.
     */
    public LGraph importGraph(final ElkNode elkgraph) {
        // Create the layered graph
        final LGraph topLevelGraph = createLGraph(elkgraph);
        
        List<ElkPort> ports = elkgraph.getPorts();
        // Assign defined port sides to all external ports 
        elkgraph.getPorts().stream().forEach(elkport -> ensureDefinedPortSide(topLevelGraph, elkport));
        
        // Transform the external ports, if any
        Set<GraphProperties> graphProperties = topLevelGraph.getProperty(InternalProperties.GRAPH_PROPERTIES);
        checkExternalPorts(elkgraph, graphProperties);
        if (graphProperties.contains(GraphProperties.EXTERNAL_PORTS)) {
            for (ElkPort elkport : elkgraph.getPorts()) {
                transformExternalPort(elkgraph, topLevelGraph, elkport);
            }
        }
        
        // Calculate the graph's minimum size
        if (shouldCalculateMinimumGraphSize(elkgraph)) {
            calculateMinimumGraphSize(elkgraph, topLevelGraph);
        }
        
        // Remember things
        if (topLevelGraph.getProperty(LayeredOptions.PARTITIONING_ACTIVATE)) {
            graphProperties.add(GraphProperties.PARTITIONS);
        }
        
        // Apply a spacing configuration based on a base value (if it has been requested)
        //  Note that the computed spacing values are set on the lgraph and not the elkgraph to avoid polluting 
        //  the input graph. If the spacing values were set on the input graph, a second layout run of the same 
        //  input graph - with a different base value - would yield an unexpected result as the computed spacing 
        //  values of the first layout run would be used (explicitly set spacing values are not overwritten).
        if (topLevelGraph.hasProperty(LayeredOptions.SPACING_BASE_VALUE)) {
            LayeredSpacings.withBaseValue(topLevelGraph.getProperty(LayeredOptions.SPACING_BASE_VALUE))
                    .apply(topLevelGraph);
        }

        // Import the graph either with or without multiple nested levels of hierarchy
        if (elkgraph.getProperty(LayeredOptions.HIERARCHY_HANDLING) == HierarchyHandling.INCLUDE_CHILDREN) {
            importHierarchicalGraph(elkgraph, topLevelGraph);
        } else {
            importFlatGraph(elkgraph, topLevelGraph);
        }
        
        return topLevelGraph;
    }
    
    /**
     * Ensures that the given port has a defined port side.
     */
    private void ensureDefinedPortSide(final LGraph lgraph, final ElkPort elkport) {
        Direction layoutDirection = lgraph.getProperty(LayeredOptions.DIRECTION);
        PortSide portSide = elkport.getProperty(LayeredOptions.PORT_SIDE);
        PortConstraints portConstraints = lgraph.getProperty(LayeredOptions.PORT_CONSTRAINTS);
        
        if (!portConstraints.isSideFixed()) {
            // We are free to assign ports to sides, so the port side will depend on the layout direction and the
            // port's net flow
            int netFlow = calculateNetFlow(elkport);
            
            if (netFlow > 0) {
                portSide = PortSide.fromDirection(layoutDirection);
            } else {
                portSide = PortSide.fromDirection(layoutDirection).opposed();
            }
            
        } else {
            // We are not free to assign port sides. If none is set, try inferring it from the port's position
            if (portSide == PortSide.UNDEFINED) {
                portSide = ElkUtil.calcPortSide(elkport, layoutDirection);
                
                // There are cases where ELK may have failed to infer the port side
                if (portSide == PortSide.UNDEFINED) {
                    portSide = PortSide.fromDirection(layoutDirection);
                }
            }
        }
        
        elkport.setProperty(LayeredOptions.PORT_SIDE, portSide);
    }
    
    /**
     * Checks whether {@link #calculateMinimumGraphSize(ElkNode, LGraph)} should be called on the given graph.
     */
    private boolean shouldCalculateMinimumGraphSize(final ElkNode elkgraph) {
        return !elkgraph.getProperty(LayeredOptions.NODE_SIZE_CONSTRAINTS).isEmpty();
    }
    
    /**
     * Asks the label and node size thing to calculate the minimum size necessary for the graph to be large enough for
     * its ports and stuff (if it's not the top level graph).
     * 
     * @param elkgraph
     *            the original ELK graph.
     * @param lgraph
     *            the imported LGraph. Its properties may be updated over the course of this method.
     */
    private void calculateMinimumGraphSize(final ElkNode elkgraph, final LGraph lgraph) {
        // If the graph is on the top level, don't bother
        if (elkgraph.getParent() == null) {
            return;
        }
        
        // If the graph has no size constraints, don't bother either
        EnumSet<SizeConstraint> sizeConstraints = lgraph.getProperty(LayeredOptions.NODE_SIZE_CONSTRAINTS);
        
        // The method should only be called if shouldCalculateMinimumGraphSize(...) returns true
        assert !sizeConstraints.isEmpty();
        
        // Ensure that the port constraints are not UNDEFINED
        if (elkgraph.getProperty(LayeredOptions.PORT_CONSTRAINTS) == PortConstraints.UNDEFINED) {
            elkgraph.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FREE);
        }
        
        // Size constraints are not empty, so calculate the size the node and label placement code thing would like to
        // give the graph
        GraphAdapter<?> graphAdapter = ElkGraphAdapters.adapt(elkgraph.getParent());
        NodeAdapter<?> nodeAdapter = ElkGraphAdapters.adaptSingleNode(elkgraph);
        
        KVector minSize = NodeLabelAndSizeCalculator.process(graphAdapter, nodeAdapter, false, true);
        
        // Apply the minimum size a sa property and make sure the minimum size is respected by ELK Layered by making
        // sure the necessary size constraint exists
        sizeConstraints.add(SizeConstraint.MINIMUM_SIZE);
        
        KVector configuredMinSize = lgraph.getProperty(LayeredOptions.NODE_SIZE_MINIMUM);
        configuredMinSize.x = Math.max(minSize.x, configuredMinSize.x);
        configuredMinSize.y = Math.max(minSize.y, configuredMinSize.y);
    }

    /**
     * Imports the direct children of the given graph.
     * 
     * @param elkgraph
     *            graph to import.
     * @param lgraph
     *            graph to add the imported elements to.
     */
    private void importFlatGraph(final ElkNode elkgraph, final LGraph lgraph) {
        // Transform the node's children, unless we're told not to
        int index = 0;
        HashSet<Integer> cbGroupModelOrders = new HashSet<>();
        for (ElkNode child : elkgraph.getChildren()) {
            if (!child.getProperty(LayeredOptions.NO_LAYOUT)) {
                if (needsModelOrder(child)) {
                    child.setProperty(InternalProperties.MODEL_ORDER, index);
                    index++;
                    if (child.hasProperty(LayeredOptions.CONSIDER_MODEL_ORDER_GROUP_MODEL_ORDER_CYCLE_BREAKING_ID)) {
                        cbGroupModelOrders.add(child.getProperty(LayeredOptions.CONSIDER_MODEL_ORDER_GROUP_MODEL_ORDER_CYCLE_BREAKING_ID));
                    }
                }
                transformNode(child, lgraph);
            }
        }
        // Save the maximum node model order.
        // This is relevant to create graph partitions based on model order and constraints.
        lgraph.setProperty(InternalProperties.MAX_MODEL_ORDER_NODES, index);
        // Save the number of model order groups.
        lgraph.setProperty(InternalProperties.CB_NUM_MODEL_ORDER_GROUPS, cbGroupModelOrders.size());

        // iterate the list of contained edges to preserve the 'input order' of the edges
        // (this is not part of the previous loop since all children must have already been transformed)
        index = 0;
        for (ElkEdge elkedge : elkgraph.getContainedEdges()) {
            if (needsModelOrderBasedOnParent(elkgraph)) {
                elkedge.setProperty(InternalProperties.MODEL_ORDER, index);
                index++;
            }
            ElkNode source = ElkGraphUtil.getSourceNode(elkedge);
            ElkNode target = ElkGraphUtil.getTargetNode(elkedge);
            
            // Is inside self loop processing enabled for this node?
            boolean enableInsideSelfLoops = source.getProperty(LayeredOptions.INSIDE_SELF_LOOPS_ACTIVATE);
            
            // Find out basic information about the edge
            boolean isToBeLaidOut = !elkedge.getProperty(LayeredOptions.NO_LAYOUT);
            boolean isInsideSelfLoop = enableInsideSelfLoops && elkedge.isSelfloop()
                    && elkedge.getProperty(LayeredOptions.INSIDE_SELF_LOOPS_YO);
            boolean connectsSiblings = source.getParent() == elkgraph && source.getParent() == target.getParent();
            boolean connectsToGraph = (source.getParent() == elkgraph && target == elkgraph) 
                    ^ (target.getParent() == elkgraph && source == elkgraph);
            
            // Only transform the edge if we are to layout the edge and if it stays in the current
            // level of hierarchy (which implies that here we don't transform inside self loops)
            if (isToBeLaidOut && !isInsideSelfLoop && (connectsToGraph || connectsSiblings)) {
                transformEdge(elkedge, elkgraph, lgraph);
            } 
        }

        // now collect inside self loops of 'elkgraph'
        if (elkgraph.getParent() != null) {
            for (ElkEdge elkedge : elkgraph.getParent().getContainedEdges()) {
                ElkNode source = ElkGraphUtil.getSourceNode(elkedge);
                if (source == elkgraph && elkedge.isSelfloop()) {
                    boolean isInsideSelfLoop = source.getProperty(LayeredOptions.INSIDE_SELF_LOOPS_ACTIVATE)
                            && elkedge.getProperty(LayeredOptions.INSIDE_SELF_LOOPS_YO);
                    if (isInsideSelfLoop) {
                        transformEdge(elkedge, elkgraph, lgraph);
                    }
                }
            }
        }
        
    }

    /**
     * Imports the graph hierarchy rooted at the given graph.
     * 
     * @param elkgraph
     *            graph to import.
     * @param lgraph
     *            graph to add the direct children of the current hierarchy level to.
     */
    private void importHierarchicalGraph(final ElkNode elkgraph, final LGraph lgraph) {
        final Queue<ElkNode> elkGraphQueue = Lists.newLinkedList();
        
        Direction parentGraphDirection = lgraph.getProperty(LayeredOptions.DIRECTION);

        // Model order index for nodes
        int index = 0;
        HashSet<Integer> cbGroupModelOrders = new HashSet<>();
        // Transform the node's children
        elkGraphQueue.addAll(elkgraph.getChildren());
        while (!elkGraphQueue.isEmpty()) {
            ElkNode elknode = elkGraphQueue.poll();
            
            if (needsModelOrder(elknode)) {
                // Assign a model order to the nodes as they are read
                elknode.setProperty(InternalProperties.MODEL_ORDER, index++);
                if (elknode.hasProperty(LayeredOptions.CONSIDER_MODEL_ORDER_GROUP_MODEL_ORDER_CYCLE_BREAKING_ID)) {
                    cbGroupModelOrders.add(elknode.getProperty(LayeredOptions.CONSIDER_MODEL_ORDER_GROUP_MODEL_ORDER_CYCLE_BREAKING_ID));
                }
            }
            
            // Check if the current node is to be laid out in the first place
            boolean isNodeToBeLaidOut = !elknode.getProperty(LayeredOptions.NO_LAYOUT);
            if (isNodeToBeLaidOut) {
                // Check if there has to be an LGraph for this node (which is the case if it has children or inside
                // self-loops, and if it does not have another layout algorithm configured)
                boolean hasChildren = !elknode.getChildren().isEmpty();
                boolean hasInsideSelfLoops = hasInsideSelfLoops(elknode);
                boolean hasHierarchyHandlingEnabled = elknode.getProperty(LayeredOptions.HIERARCHY_HANDLING)
                        == HierarchyHandling.INCLUDE_CHILDREN;
                boolean usesElkLayered = !elknode.hasProperty(CoreOptions.ALGORITHM)
                        || LayeredOptions.ALGORITHM_ID.endsWith(elknode.getProperty(CoreOptions.ALGORITHM));

                LGraph nestedGraph = null;
                if (usesElkLayered && hasHierarchyHandlingEnabled && (hasChildren || hasInsideSelfLoops)) {
                    nestedGraph = createLGraph(elknode);
                    nestedGraph.setProperty(LayeredOptions.DIRECTION, parentGraphDirection);
                    
                    // Apply a spacing configuration, for details see comment int #importGraph(...)
                    if (nestedGraph.hasProperty(LayeredOptions.SPACING_BASE_VALUE)) {
                        LayeredSpacings.withBaseValue(nestedGraph.getProperty(LayeredOptions.SPACING_BASE_VALUE))
                                .apply(nestedGraph);
                    }
                    
                    // We need to make sure that we make the graph large enough for any ports, node labels, etc.
                    // if the size constraints are not empty
                    if (shouldCalculateMinimumGraphSize(elknode)) {
                        final LGraph finalNestedGraph = nestedGraph;
                        elknode.getPorts().stream()
                            .forEach(elkport -> ensureDefinedPortSide(finalNestedGraph, elkport));
                        calculateMinimumGraphSize(elknode, nestedGraph);
                    }
                }
                
                // Transform da node!!!
                LGraph parentLGraph = lgraph;
                LNode parentLNode = (LNode) nodeAndPortMap.get(elknode.getParent());
                if (parentLNode != null) {
                    parentLGraph = parentLNode.getNestedGraph();
                }
                LNode lnode = transformNode(elknode, parentLGraph);
                
                // Setup hierarchical relationships
                if (nestedGraph != null) {
                    lnode.setNestedGraph(nestedGraph);
                    nestedGraph.setParentNode(lnode);
                    
                    elkGraphQueue.addAll(elknode.getChildren());
                }
            }
        }
        // Save the maximum node model order.
        // This is relevant to create graph partitions based on model order and constraints.
        lgraph.setProperty(InternalProperties.MAX_MODEL_ORDER_NODES, index);
        // Save the number of model order groups.
        lgraph.setProperty(InternalProperties.CB_NUM_MODEL_ORDER_GROUPS, cbGroupModelOrders.size());

        // Model order index for edges. 
        index = 0;
        // Transform the edges
        elkGraphQueue.add(elkgraph);
        while (!elkGraphQueue.isEmpty()) {
            ElkNode elkGraphNode = elkGraphQueue.poll();
            
            for (ElkEdge elkedge : elkGraphNode.getContainedEdges()) {
                // We don't support hyperedges
                checkEdgeValidity(elkedge);
                
                if (needsModelOrderBasedOnParent(elkgraph)) {
                    // Assign a model order to the edges as they are read
                    elkedge.setProperty(InternalProperties.MODEL_ORDER, index++);
                }
                
                ElkNode sourceNode = ElkGraphUtil.connectableShapeToNode(elkedge.getSources().get(0));
                ElkNode targetNode = ElkGraphUtil.connectableShapeToNode(elkedge.getTargets().get(0));
                
                // Don't bother if either the edge or at least one of its end points are excluded from layout
                if (elkedge.getProperty(LayeredOptions.NO_LAYOUT)
                        || sourceNode.getProperty(LayeredOptions.NO_LAYOUT)
                        || targetNode.getProperty(LayeredOptions.NO_LAYOUT)) {
                    continue;
                }
                
                // Check if this edge is an inside self-loop
                boolean isInsideSelfLoop = elkedge.isSelfloop()
                        && sourceNode.getProperty(LayeredOptions.INSIDE_SELF_LOOPS_ACTIVATE)
                        && elkedge.getProperty(LayeredOptions.INSIDE_SELF_LOOPS_YO);
                
                // Find the graph the edge will be placed in. Basically, if the edge is an inside
                // self loop or connects one of its end points to a descendant, the edge will be
                // placed in the graph that represents that end point's insides. Otherwise, it will
                // be placed in the current graph.
                ElkNode parentElkGraph = elkGraphNode;
                
                if (isInsideSelfLoop || ElkGraphUtil.isDescendant(targetNode, sourceNode)) {
                    parentElkGraph = sourceNode;
                } else if (ElkGraphUtil.isDescendant(sourceNode, targetNode)) {
                    parentElkGraph = targetNode;
                }
                
                LGraph parentLGraph = lgraph;
                LNode parentLNode = (LNode) nodeAndPortMap.get(parentElkGraph);
                if (parentLNode != null) {
                    parentLGraph = parentLNode.getNestedGraph();
                }
                
                // Transform the edge, finally...
                LEdge ledge = transformEdge(elkedge, parentElkGraph, parentLGraph);
                
                // Find the graph the edge's coordinates will have to be made relative to during export. This will only
                // do something if the edge containment inside ELK Layered differs from the edge containment in the
                // ELK graph
                ledge.setProperty(InternalProperties.COORDINATE_SYSTEM_ORIGIN,
                        findCoordinateSystemOrigin(elkedge, elkgraph, lgraph));
            }
            
            // We may need to look at edges contained in the current graph node's children as well.
            // this is true unless either the current graph node does not have hierarchy handling
            // enabled, or a child has another layout algorithm configured
            boolean hasHierarchyHandlingEnabled = elkGraphNode.getProperty(LayeredOptions.HIERARCHY_HANDLING)
                    == HierarchyHandling.INCLUDE_CHILDREN;
            if (hasHierarchyHandlingEnabled) {
                for (ElkNode elkChildGraphNode : elkGraphNode.getChildren()) {
                    boolean usesElkLayered =
                            !elkChildGraphNode.hasProperty(CoreOptions.ALGORITHM) || LayeredOptions.ALGORITHM_ID
                                    .endsWith(elkChildGraphNode.getProperty(CoreOptions.ALGORITHM));
                    boolean partOfSameLayoutRun = elkChildGraphNode.getProperty(LayeredOptions.HIERARCHY_HANDLING)
                            == HierarchyHandling.INCLUDE_CHILDREN;
                    
                    if (usesElkLayered && partOfSameLayoutRun) {
                        elkGraphQueue.add(elkChildGraphNode);
                    }
                }
            }
        }
    }
    
    /**
     * Calculates whether a node needs to have a model order set.
     * 
     * @param child The node.
     * @return True, if model order should be set.
     */
    private boolean needsModelOrder(final ElkNode child) {
        ElkNode elkgraph = child.getParent();
        return needsModelOrderBasedOnParent(elkgraph)
                && !child.getProperty(LayeredOptions.CONSIDER_MODEL_ORDER_NO_MODEL_ORDER);
    }
    
    /**
     * Calculates whether the parent graph of an element would need model order to be set.
     * 
     * @param elkgraph The graph.
     * @return True, if model order should be set.
     */
    private boolean needsModelOrderBasedOnParent(final ElkNode elkgraph) {
        
        boolean modelOrderCycleBreaking = elkgraph.getProperty(LayeredOptions.CYCLE_BREAKING_STRATEGY) == CycleBreakingStrategy.MODEL_ORDER
                || elkgraph.getProperty(LayeredOptions.CYCLE_BREAKING_STRATEGY) == CycleBreakingStrategy.BFS_NODE_ORDER
                || elkgraph.getProperty(LayeredOptions.CYCLE_BREAKING_STRATEGY) == CycleBreakingStrategy.DFS_NODE_ORDER
                || elkgraph.getProperty(LayeredOptions.CYCLE_BREAKING_STRATEGY) == CycleBreakingStrategy.GREEDY_MODEL_ORDER
                || elkgraph.getProperty(LayeredOptions.CYCLE_BREAKING_STRATEGY) == CycleBreakingStrategy.SCC_CONNECTIVITY
                || elkgraph.getProperty(LayeredOptions.CYCLE_BREAKING_STRATEGY) == CycleBreakingStrategy.SCC_NODE_TYPE;
        boolean modelOrderLayering = elkgraph.getProperty(LayeredOptions.LAYERING_STRATEGY) == LayeringStrategy.BF_MODEL_ORDER
                || elkgraph.getProperty(LayeredOptions.LAYERING_STRATEGY) == LayeringStrategy.DF_MODEL_ORDER
                || elkgraph.getProperty(LayeredOptions.LAYERING_NODE_PROMOTION_STRATEGY) == NodePromotionStrategy.MODEL_ORDER_LEFT_TO_RIGHT
                || elkgraph.getProperty(LayeredOptions.LAYERING_NODE_PROMOTION_STRATEGY) == NodePromotionStrategy.MODEL_ORDER_RIGHT_TO_LEFT;
        boolean modelOrderCrossingMinimization =
                // Maybe add the explicit strategies here
                elkgraph.getProperty(LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY) != OrderingStrategy.NONE
                || elkgraph.getProperty(LayeredOptions.CROSSING_MINIMIZATION_FORCE_NODE_MODEL_ORDER)
                // Maybe add the explicit strategies here
                || elkgraph.getProperty(LayeredOptions.CONSIDER_MODEL_ORDER_COMPONENTS) != ComponentOrderingStrategy.NONE
                
                || elkgraph.getProperty(LayeredOptions.CONSIDER_MODEL_ORDER_CROSSING_COUNTER_NODE_INFLUENCE) != 0
                || elkgraph.getProperty(LayeredOptions.CONSIDER_MODEL_ORDER_CROSSING_COUNTER_PORT_INFLUENCE) != 0;
        return modelOrderCycleBreaking || modelOrderLayering || modelOrderCrossingMinimization;
    }    
    

    /**
     * Checks if the given node has any inside self loops.
     * 
     * @param elknode the node to check for inside self loops.
     * @return {@code true} if the node has inside self loops, {@code false} otherwise.
     */
    private boolean hasInsideSelfLoops(final ElkNode elknode) {
        if (elknode.getProperty(LayeredOptions.INSIDE_SELF_LOOPS_ACTIVATE)) {
            for (ElkEdge edge : ElkGraphUtil.allOutgoingEdges(elknode)) {
                if (edge.isSelfloop()) {
                    if (edge.getProperty(LayeredOptions.INSIDE_SELF_LOOPS_YO)) {
                        return true;
                    }
                }
            }
        }
        
        return false;
    }
    
    /**
     * Finds the LGraph the edge's coordinates should be relative to when the layout results are applied back. This
     * is only relevant if this differs from the graph the coordinates are relative to inside ELK Layered. In fact,
     * this method only returns something for edges that connect nodes that are not in an anscestor-descendant
     * relationship.
     */
    private LGraph findCoordinateSystemOrigin(final ElkEdge elkedge, final ElkNode topLevelElkGraph,
            final LGraph topLevelLGraph) {
        
        ElkNode source = ElkGraphUtil.connectableShapeToNode(elkedge.getSources().get(0));
        ElkNode target = ElkGraphUtil.connectableShapeToNode(elkedge.getTargets().get(0));
        
        // We're going to rule out one case after the other. If the source and the target are siblings, we're good
        // (this also includes self-loops)
        if (source.getParent() == target.getParent()) {
            return null;
        }
        
        // If the target is a descendant of the source, ELK Layered uses the source's top left corner as the origin
        // of the coordinate system, which matches how ELK graph should be constructed
        if (ElkGraphUtil.isDescendant(target, source)) {
            return null;
        }
        
        // If the source is a descendant of the target, ELK Layered uses the source's parent graph as the origin of
        // the coordinate system, while ELK will expect the first common ancestor (the target) to be the origin.
        // 
        // If source and target have no relationship to each other, ELK Layered again uses the source's parent graph
        // as the origin of the coordinate system, while ELK will expect the first common ancestor to be the origin
        ElkNode origin = elkedge.getContainingNode();
        
        // The origin must always be an ancestor of both, source and parent
        assert source == origin || ElkGraphUtil.isDescendant(source, origin);
        assert target == origin || ElkGraphUtil.isDescendant(target, origin);
        
        // Find the associated LGraph
        if (origin == topLevelElkGraph) {
            return topLevelLGraph;
        } else {
            LNode lnode = (LNode) nodeAndPortMap.get(origin);
            if (lnode != null) {
                // Find the graph that represents the node's insides
                LGraph lgraph = lnode.getNestedGraph();
                if (lgraph != null) {
                    return lgraph;
                }
            }
        }
        
        return null;
    }
    

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Graph Transformation

    /**
     * Create an LGraph from the given node.
     * 
     * @param elkgraph
     *            the parent node from which to create the LGraph
     * @return a new LGraph instance
     */
    private LGraph createLGraph(final ElkNode elkgraph) {
        LGraph lgraph = new LGraph();
        
        // Copy the properties of the KGraph to the layered graph
        lgraph.copyProperties(elkgraph);
        if (lgraph.getProperty(LayeredOptions.DIRECTION) == Direction.UNDEFINED) {
            lgraph.setProperty(LayeredOptions.DIRECTION, LGraphUtil.getDirection(lgraph));
        }
        
        // The root may have a label manager installed
        if (lgraph.getProperty(LabelManagementOptions.LABEL_MANAGER) == null) {
            ElkGraphElement root = (ElkGraphElement) EcoreUtil.getRootContainer(elkgraph);
            lgraph.setProperty(LabelManagementOptions.LABEL_MANAGER,
                    root.getProperty(LabelManagementOptions.LABEL_MANAGER));
        }
        
        // Remember the KGraph parent the LGraph was created from
        lgraph.setProperty(InternalProperties.ORIGIN, elkgraph);

        // Initialize the graph properties discovered during the transformations
        lgraph.setProperty(InternalProperties.GRAPH_PROPERTIES,
                EnumSet.noneOf(GraphProperties.class));
        
        // Adjust the padding to respect inside labels (if the graph has a parent, we need to supply that as well
        // since size information stored there may apply to the current graph node)
        ElkPadding nodeLabelpadding = NodeLabelAndSizeCalculator.computeInsideNodeLabelPadding(
                elkgraph.getParent() == null ? null : ElkGraphAdapters.adapt(elkgraph.getParent()),
                ElkGraphAdapters.adaptSingleNode(elkgraph),
                Direction.RIGHT);
        ElkPadding nodePadding = lgraph.getProperty(LayeredOptions.PADDING);

        // Setup the graph's padding
        LPadding lPadding = lgraph.getPadding();
        lPadding.add(nodePadding);
        lPadding.add(nodeLabelpadding);

        return lgraph;
    }
    

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // External Port Transformation
    
    /**
     * Checks if external ports processing should be active. This is the case if the parent node has
     * ports and at least one of the following conditions is true:
     * <ul>
     *   <li>
     *     Port label placement is set to {@code INSIDE} and at least one of the ports has a label.
     *   </li>
     *   <li>
     *     At least one of the ports has an edge that connects to the insides of the parent node.
     *   </li>
     *   <li>
     *     There is a self-loop that should be routed inside the node.
     *   </li>
     * </ul>
     * 
     * @param elkgraph
     *            a KGraph we want to check for external ports.
     * @param graphProperties
     *            the set of graph properties to store our results in.
     */
    private void checkExternalPorts(final ElkNode elkgraph, final Set<GraphProperties> graphProperties) {
        final boolean enableSelfLoops = elkgraph.getProperty(LayeredOptions.INSIDE_SELF_LOOPS_ACTIVATE);
        final Set<PortLabelPlacement> portLabelPlacement = elkgraph.getProperty(LayeredOptions.PORT_LABELS_PLACEMENT);

        // We're iterating over the ports until we've determined that we have both external ports and
        // hyperedges, or if there are no more ports left
        boolean hasExternalPorts = false;
        boolean hasHyperedges = false;
        
        final Iterator<ElkPort> portIterator = elkgraph.getPorts().iterator();
        while (portIterator.hasNext() && (!hasExternalPorts || !hasHyperedges)) {
            final ElkPort elkport = portIterator.next();
            
            // Find out if there are edges connected to external ports of the graph (this is the case
            // for inside self loops as well as for edges connected to children)
            int externalPortEdges = 0;
            
            for (ElkEdge elkedge : ElkGraphUtil.allIncidentEdges(elkport)) {
                boolean isInsideSelfLoop = enableSelfLoops && elkedge.isSelfloop()
                        && elkedge.getProperty(LayeredOptions.INSIDE_SELF_LOOPS_YO);
                boolean connectsToChild = elkedge.getSources().contains(elkport)
                        ? elkgraph == ElkGraphUtil.connectableShapeToNode(elkedge.getTargets().get(0)).getParent()
                        : elkgraph == ElkGraphUtil.connectableShapeToNode(elkedge.getSources().get(0)).getParent();
                
                if (isInsideSelfLoop || connectsToChild) {
                    externalPortEdges++;
                    if (externalPortEdges > 1) {
                        break;
                    }
                }
            }
            
            // External ports?
            if (externalPortEdges > 0) {
                hasExternalPorts = true;
            } else if (portLabelPlacement.contains(PortLabelPlacement.INSIDE) && elkport.getLabels().size() > 0) {
                hasExternalPorts = true;
            }
            
            // Hyperedges, even?
            if (externalPortEdges > 1) {
                hasHyperedges = true;
            }
        }
        
        // Update graph properties
        if (hasExternalPorts) {
            graphProperties.add(GraphProperties.EXTERNAL_PORTS);
        }
        
        if (hasHyperedges) {
            graphProperties.add(GraphProperties.HYPEREDGES);
        }
    }

    /**
     * Transforms the given external port into a dummy node.
     * 
     * @param elkgraph
     *            the original KGraph
     * @param lgraph
     *            the corresponding layered graph
     * @param elkport
     *            the port to be transformed
     */
    private void transformExternalPort(final ElkNode elkgraph, final LGraph lgraph, final ElkPort elkport) {
        // We need some information about the port
        KVector elkportPosition = new KVector(
                elkport.getX() + elkport.getWidth() / 2.0,
                elkport.getY() + elkport.getHeight() / 2.0);
        int netFlow = calculateNetFlow(elkport);
        PortConstraints portConstraints = elkgraph.getProperty(LayeredOptions.PORT_CONSTRAINTS);
        
        // If we don't have a proper port side, calculate one
        PortSide portSide = elkport.getProperty(LayeredOptions.PORT_SIDE);
        assert portSide != PortSide.UNDEFINED;
        
        // If we don't have a port offset, infer one
        if (!elkport.getAllProperties().containsKey(LayeredOptions.PORT_BORDER_OFFSET)) {
            double portOffset;
            // if port coordinates are (0,0), we default to port offset 0 to make the common case frustration-free
            if (elkport.getX() == 0.0 && elkport.getY() == 0.0) {
                portOffset = 0.0;
            } else {
                portOffset = ElkUtil.calcPortOffset(elkport, portSide);
            }
            elkport.setProperty(LayeredOptions.PORT_BORDER_OFFSET, portOffset);
        }
        
        // Create the external port dummy node
        KVector graphSize = new KVector(elkgraph.getWidth(), elkgraph.getHeight());
        LNode dummy = LGraphUtil.createExternalPortDummy(
                elkport, portConstraints, portSide, netFlow, graphSize,
                elkportPosition, new KVector(elkport.getWidth(), elkport.getHeight()),
                lgraph.getProperty(LayeredOptions.DIRECTION), lgraph);
        dummy.setProperty(InternalProperties.ORIGIN, elkport);
        
        // The dummy only has one port
        LPort dummyPort = dummy.getPorts().get(0);
        dummyPort.setConnectedToExternalNodes(isConnectedToExternalNodes(elkport));
        dummy.setProperty(LayeredOptions.PORT_LABELS_PLACEMENT, PortLabelPlacement.outside());
        
        // If the compound node wants to have its port labels placed on the inside, we need to leave
        // enough space for them by creating an LLabel for the KLabels. If the compound node wants to
        // have its port labels placed on the outside, we still need to leave enough space for them
        // so the port placement does not cause problems on the outside, but we also don't want to waste
        // space inside. Thus, for east and west ports, we reduce the label width to zero, otherwise
        // we reduce the label height to zero
        boolean insidePortLabels =
                elkgraph.getProperty(LayeredOptions.PORT_LABELS_PLACEMENT).contains(PortLabelPlacement.INSIDE);

        // Transform all of the port's labels
        for (ElkLabel elklabel : elkport.getLabels()) {
            if (!elklabel.getProperty(LayeredOptions.NO_LAYOUT) && !Strings.isNullOrEmpty(elklabel.getText())) {
                LLabel llabel = transformLabel(elklabel);
                dummyPort.getLabels().add(llabel);
                
                // If port labels are placed outside, modify the size
                // If the port labels are fixed, we should consider the part that is inside the node and not 0.
                if (!insidePortLabels) {
                    double insidePart = 0;
                    if (PortLabelPlacement.isFixed(elkgraph.getProperty(LayeredOptions.PORT_LABELS_PLACEMENT))) {
                        // We use 0 as port border offset here, as we only want the label part that is
                        // inside the node "after" the port.
                        insidePart = ElkUtil.computeInsidePart(new KVector(elklabel.getX(), elklabel.getY()),
                                new KVector(elklabel.getWidth(), elklabel.getHeight()),
                                new KVector(elkport.getWidth(), elkport.getHeight()), 0, portSide);
                    }
                    switch (portSide) {
                    case EAST:
                    case WEST:
                        llabel.getSize().x = insidePart;
                        break;

                    case NORTH:
                    case SOUTH:
                        llabel.getSize().y = insidePart;
                        break;
                    }
                }
            }
        }
        
        // Remember the relevant spacings that will apply to the labels here. It's not the spacings in the graph, but
        // in the parent
        dummy.setProperty(LayeredOptions.SPACING_LABEL_PORT_HORIZONTAL,
                elkgraph.getParent().getProperty(LayeredOptions.SPACING_LABEL_PORT_HORIZONTAL));
        dummy.setProperty(LayeredOptions.SPACING_LABEL_PORT_VERTICAL,
                elkgraph.getParent().getProperty(LayeredOptions.SPACING_LABEL_PORT_VERTICAL));
        dummy.setProperty(LayeredOptions.SPACING_LABEL_LABEL,
                elkgraph.getParent().getProperty(LayeredOptions.SPACING_LABEL_LABEL));
        
        // Put the external port dummy into our graph and associate it with the original KPort
        lgraph.getLayerlessNodes().add(dummy);
        nodeAndPortMap.put(elkport, dummy);
    }
    
    /**
     * Count how many edges want the port to be an output port of the parent and how many want it to
     * be an input port. An edge coming into the port from the inside votes for the port to be an
     * output port of the parent, as does an edge leaving the port for the outside. The result returned
     * by this method is the so-called net flow as fed into {@code createExternalPort(..)}.
     * 
     * @param elkport
     *            the port to look at.
     * @return the port's net flow.
     */
    private int calculateNetFlow(final ElkPort elkport) {
        final ElkNode elkgraph = elkport.getParent();
        final boolean insideSelfLoopsEnabled = elkgraph.getProperty(LayeredOptions.INSIDE_SELF_LOOPS_ACTIVATE);

        int outputPortVote = 0, inputPortVote = 0;
        
        // Iterate over outgoing edges
        for (ElkEdge outgoingEdge : elkport.getOutgoingEdges()) {
            final boolean isSelfLoop = outgoingEdge.isSelfloop();
            final boolean isInsideSelfLoop = isSelfLoop && insideSelfLoopsEnabled
                    && outgoingEdge.getProperty(LayeredOptions.INSIDE_SELF_LOOPS_YO);
            final ElkNode targetNode = ElkGraphUtil.connectableShapeToNode(outgoingEdge.getTargets().get(0));

            if (isSelfLoop && isInsideSelfLoop) {
                inputPortVote++;
            } else if (isSelfLoop && !isInsideSelfLoop) {
                outputPortVote++;
            } else if (targetNode.getParent() == elkgraph || targetNode  == elkgraph) {
                inputPortVote++;
            } else {
                outputPortVote++;
            }
        }
        
        // Iterate over incoming edges
        for (ElkEdge incomingEdge : elkport.getIncomingEdges()) {
            final boolean isSelfLoop = incomingEdge.isSelfloop();
            final boolean isInsideSelfLoop = isSelfLoop && insideSelfLoopsEnabled
                    && incomingEdge.getProperty(LayeredOptions.INSIDE_SELF_LOOPS_YO);
            final ElkNode sourceNode = ElkGraphUtil.connectableShapeToNode(incomingEdge.getSources().get(0));

            if (isSelfLoop && isInsideSelfLoop) {
                outputPortVote++;
            } else if (isSelfLoop && !isInsideSelfLoop) {
                inputPortVote++;
            } else if (sourceNode.getParent() == elkgraph || sourceNode  == elkgraph) {
                outputPortVote++;
            } else {
                inputPortVote++;
            }
        }
        
        return outputPortVote - inputPortVote;
    }
    
    /**
     * Checks whether the given (external) port has connections to the outside (that is, to non-descendants).
     */
    private boolean isConnectedToExternalNodes(final ElkPort elkport) {
        ElkNode parent = elkport.getParent();
        
        for (ElkEdge outEdge : elkport.getOutgoingEdges()) {
            ElkNode targetNode = ElkGraphUtil.connectableShapeToNode(outEdge.getTargets().get(0));
            
            if (!ElkGraphUtil.isDescendant(targetNode, parent)) {
                return true;
            }
        }
        
        for (ElkEdge inEdge : elkport.getIncomingEdges()) {
            ElkNode sourceNode = ElkGraphUtil.connectableShapeToNode(inEdge.getSources().get(0));
            
            if (!ElkGraphUtil.isDescendant(sourceNode, parent)) {
                return true;
            }
        }
        
        return false;
    }
    

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Node Transformation

    /**
     * Transforms the given node and its contained ports.
     * 
     * @param elknode
     *            the node to transform
     * @param lgraph
     *            the layered graph into which the transformed node is put
     * @return the transformed node
     */
    private LNode transformNode(final ElkNode elknode, final LGraph lgraph) {
        // add a new node to the layered graph, copying its position
        LNode lnode = new LNode(lgraph);
        lnode.copyProperties(elknode);
        lnode.setProperty(InternalProperties.ORIGIN, elknode);
        
        lnode.getSize().x = elknode.getWidth();
        lnode.getSize().y = elknode.getHeight();
        lnode.getPosition().x = elknode.getX();
        lnode.getPosition().y = elknode.getY();
        
        lgraph.getLayerlessNodes().add(lnode);
        nodeAndPortMap.put(elknode, lnode);
        
        // check if the node is a compound node in the original graph
        if (!elknode.getChildren().isEmpty() || elknode.getProperty(LayeredOptions.INSIDE_SELF_LOOPS_ACTIVATE)) {
            lnode.setProperty(InternalProperties.COMPOUND_NODE, true);
        }

        Set<GraphProperties> graphProperties = lgraph.getProperty(InternalProperties.GRAPH_PROPERTIES);
        
        // port constraints and sides cannot be undefined
        PortConstraints portConstraints = lnode.getProperty(LayeredOptions.PORT_CONSTRAINTS);
        if (portConstraints == PortConstraints.UNDEFINED) {
            lnode.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FREE);
        } else if (portConstraints != PortConstraints.FREE) {
            // if the port constraints are not free, set the appropriate graph property
            graphProperties.add(GraphProperties.NON_FREE_PORTS);
        }

        int portModelOrder = 0;
        // transform the ports
        Direction direction = lgraph.getProperty(LayeredOptions.DIRECTION);
        for (ElkPort elkport : elknode.getPorts()) {
            if (needsModelOrder(elknode)) {                
                elkport.setProperty(InternalProperties.MODEL_ORDER, portModelOrder++);
            }
            if (!elkport.getProperty(LayeredOptions.NO_LAYOUT)) {
                transformPort(elkport, lnode, graphProperties, direction, portConstraints);
            }
        }

        // add the node's labels
        for (ElkLabel elklabel : elknode.getLabels()) {
            if (!elklabel.getProperty(LayeredOptions.NO_LAYOUT) && !Strings.isNullOrEmpty(elklabel.getText())) {
                lnode.getLabels().add(transformLabel(elklabel));
            }
        }

        if (lnode.getProperty(LayeredOptions.COMMENT_BOX)) {
            graphProperties.add(GraphProperties.COMMENTS);
        }

        // if we have a hypernode without ports, create a default input and output port
        if (lnode.getProperty(LayeredOptions.HYPERNODE)) {
            graphProperties.add(GraphProperties.HYPERNODES);
            graphProperties.add(GraphProperties.HYPEREDGES);
            lnode.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FREE);
        }
        
        return lnode;
    }
    

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Port Transformation
    
    /**
     * Transforms the given port. The new port will be added to the given node and will be
     * registered with the {@code transformMap}.
     * 
     * @param elkport
     *            the port to transform.
     * @param parentLNode
     *            the node the port should be added to.
     * @param graphProperties
     *            the graph properties of the graph the transformed port will be part of. The graph
     *            properties are modified depending on the port's properties.
     * @param layoutDirection
     *            the layout direction in the graph the port will be part of.
     * @param portConstraints
     *            the port constraints of the port's node.
     * @return the transformed port.
     */
    private LPort transformPort(final ElkPort elkport, final LNode parentLNode,
            final Set<GraphProperties> graphProperties, final Direction layoutDirection,
            final PortConstraints portConstraints) {
        
        // create layered port, copying its position
        LPort lport = new LPort();
        lport.copyProperties(elkport);
        lport.setSide(elkport.getProperty(LayeredOptions.PORT_SIDE));
        lport.setProperty(InternalProperties.ORIGIN, elkport);
        lport.setNode(parentLNode);
        
        KVector portSize = lport.getSize();
        portSize.x = elkport.getWidth();
        portSize.y = elkport.getHeight();
        
        KVector portPos = lport.getPosition();
        portPos.x = elkport.getX();
        portPos.y = elkport.getY();
        
        nodeAndPortMap.put(elkport, lport);
        
        // check if the original port has any outgoing connections to descendants of its node
        boolean connectionsToDescendants = elkport.getOutgoingEdges().stream()
             // All targets of each edge
            .flatMap(edge -> edge.getTargets().stream())
             // Target connectable shapes to nodes
            .map(ElkGraphUtil::connectableShapeToNode)
            // Check if any target is a descendant of the port's parent node
            .anyMatch(targetNode -> ElkGraphUtil.isDescendant(targetNode, elkport.getParent()));
        
        // there could be yet incoming connections from descendants
        if (!connectionsToDescendants) {
            // check if the original port has any incoming connections from descendants of its node
            connectionsToDescendants = elkport.getIncomingEdges().stream()
                 // All sources of each edge
                .flatMap(edge -> edge.getSources().stream())
                 // Source connectable shapes to nodes
                .map(ElkGraphUtil::connectableShapeToNode)
                // Check if any source is a descendant of the port's parent node
                .anyMatch(sourceNode -> ElkGraphUtil.isDescendant(sourceNode, elkport.getParent()));
        }
        
        // if there are still no connections to descendants, there might yet be inside self loops involved
        if (!connectionsToDescendants) {
            // check if the original port has any incoming connections from descendants of its node
            connectionsToDescendants = elkport.getOutgoingEdges().stream()
                    // All targets of each edge
                   .anyMatch(edge -> edge.isSelfloop() && edge.getProperty(LayeredOptions.INSIDE_SELF_LOOPS_YO));
        }
        
        // if we have found connections to / from descendants, mark the port accordingly
        lport.setProperty(InternalProperties.INSIDE_CONNECTIONS, connectionsToDescendants);

        // initialize the port's side, offset, and anchor point
        LGraphUtil.initializePort(lport, portConstraints, layoutDirection,
                elkport.getProperty(LayeredOptions.PORT_ANCHOR));

        // create the port's labels
        for (ElkLabel elklabel : elkport.getLabels()) {
            if (!elklabel.getProperty(LayeredOptions.NO_LAYOUT) && !Strings.isNullOrEmpty(elklabel.getText())) {
                lport.getLabels().add(transformLabel(elklabel));
            }
        }
        
        switch (layoutDirection) {
        case LEFT:
        case RIGHT:
            if (lport.getSide() == PortSide.NORTH || lport.getSide() == PortSide.SOUTH) {
                graphProperties.add(GraphProperties.NORTH_SOUTH_PORTS);
            }
            break;
        case UP:
        case DOWN:
            if (lport.getSide() == PortSide.EAST || lport.getSide() == PortSide.WEST) {
                graphProperties.add(GraphProperties.NORTH_SOUTH_PORTS);
            }
            break;
        }
        
        return lport;
    }
    

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Edge Transformation

    /**
     * Transforms the given edge if it's not a hyperedge. If it is a hyperedge, throws an exception.
     * 
     * @param elkedge the edge to transform
     * @param elkparent the node in the original graph which currently gets transformed into {@code lgraph}
     * @param lgraph the layered graph
     * @return the transformed edge, or {@code null} if it cannot be transformed
     * @throws UnsupportedGraphException if the edge is a hyperedge.
     */
    private LEdge transformEdge(final ElkEdge elkedge, final ElkNode elkparent, final LGraph lgraph) {
        checkEdgeValidity(elkedge);
        
        // Get a few basic information about the edge
        ElkConnectableShape elkSourceShape = elkedge.getSources().get(0);
        ElkConnectableShape elkTargetShape = elkedge.getTargets().get(0);
        ElkNode elkSourceNode = ElkGraphUtil.connectableShapeToNode(elkSourceShape);
        ElkNode elkTargetNode = ElkGraphUtil.connectableShapeToNode(elkTargetShape);
        
        ElkEdgeSection edgeSection = elkedge.getSections().isEmpty() ? null : elkedge.getSections().get(0);
        
        // Find the transformed source and target nodes
        LNode sourceLNode = (LNode) nodeAndPortMap.get(elkSourceNode);
        LNode targetLNode = (LNode) nodeAndPortMap.get(elkTargetNode);
        LPort sourceLPort = null;
        LPort targetLPort = null;

        // Find the transformed source port, if any
        if (elkSourceShape instanceof ElkPort) {
            // If the ElkPort is a regular port, it will map to an LPort; if it's an external port, it
            // will map to an LNode
            LGraphElement sourceElem = nodeAndPortMap.get(elkSourceShape);
            if (sourceElem instanceof LPort) {
                sourceLPort = (LPort) sourceElem;
            } else if (sourceElem instanceof LNode) {
                sourceLNode = (LNode) sourceElem;
                sourceLPort = sourceLNode.getPorts().get(0);
            }
        }

        // Find the transformed target port, if any
        if (elkTargetShape instanceof ElkPort) {
            // If the ElkPort is a regular port, it will map to an LPort; if it's an external port, it
            // will map to an LNode
            LGraphElement targetElem = nodeAndPortMap.get(elkTargetShape);
            if (targetElem instanceof LPort) {
                targetLPort = (LPort) targetElem;
            } else if (targetElem instanceof LNode) {
                targetLNode = (LNode) targetElem;
                targetLPort = targetLNode.getPorts().get(0);
            }
        }
        
        // If either the source or the target of the edge wasn't properly transformed for whatever
        // reason, we back out
        if (sourceLNode == null || targetLNode == null) {
            throw new UnsupportedGraphException("The source or the target of edge " + elkedge + " could not be found. "
                    + "This usually happens when an edge connects a node laid out by ELK Layered to a node in "
                    + "another level of hierarchy laid out by either another instance of ELK Layered or another "
                    + "layout algorithm alltogether. The former can be solved by setting the hierarchyHandling "
                    + "option to INCLUDE_CHILDREN.");
        }
        
        // Create a layered edge
        LEdge ledge = new LEdge();
        ledge.copyProperties(elkedge);
        ledge.setProperty(InternalProperties.ORIGIN, elkedge);
        
        // Clear junction points, since they are recomputed from scratch
        ledge.setProperty(LayeredOptions.JUNCTION_POINTS, null);
        
        // If we have a self-loop, set the appropriate graph property
        Set<GraphProperties> graphProperties = lgraph.getProperty(InternalProperties.GRAPH_PROPERTIES);
        if (sourceLNode == targetLNode) {
            graphProperties.add(GraphProperties.SELF_LOOPS);
        }

        // Create source and target ports if they do not exist yet
        if (sourceLPort == null) {
            PortType portType = PortType.OUTPUT;
            KVector sourcePoint = null;
            
            if (edgeSection != null && sourceLNode.getProperty(LayeredOptions.PORT_CONSTRAINTS).isSideFixed()) {
                sourcePoint = new KVector(edgeSection.getStartX(), edgeSection.getStartY());
                
                // The coordinates need to be relative to us
                ElkUtil.toAbsolute(sourcePoint, elkedge.getContainingNode());
                ElkUtil.toRelative(sourcePoint, elkparent);
                
                // If the edge is hierarchical (in which case it can only be that the target is a descendant of the
                // source), we may need to adjust the coordinates
                if (ElkGraphUtil.isDescendant(elkTargetNode, elkSourceNode)) {
                    // External source port: put it on the west side
                    portType = PortType.INPUT;
                    sourcePoint.add(sourceLNode.getPosition());
                }
            }
            sourceLPort = LGraphUtil.createPort(sourceLNode, sourcePoint, portType, lgraph);
        }
        
        if (targetLPort == null) {
            PortType portType = PortType.INPUT;
            KVector targetPoint = null;
            
            if (edgeSection != null && targetLNode.getProperty(LayeredOptions.PORT_CONSTRAINTS).isSideFixed()) {
                targetPoint = new KVector(edgeSection.getEndX(), edgeSection.getEndY());
                
                // Adjust the coordinates
                // MIGRATE Not sure yet if this really does what we want it to do
                ElkUtil.toAbsolute(targetPoint, elkedge.getContainingNode());
                ElkUtil.toRelative(targetPoint, elkparent);
            }
            
            targetLPort = LGraphUtil.createPort(targetLNode, targetPoint, portType, targetLNode.getGraph());
        }
        
        // Finally set the source and target of the edge
        ledge.setSource(sourceLPort);
        ledge.setTarget(targetLPort);
        
        // If the ports have multiple incoming or outgoing edges, the HYPEREDGE property needs to be set
        if (sourceLPort.getIncomingEdges().size() > 1 || sourceLPort.getOutgoingEdges().size() > 1
                || targetLPort.getIncomingEdges().size() > 1 || targetLPort.getOutgoingEdges().size() > 1) {
            
            graphProperties.add(GraphProperties.HYPEREDGES);
        }

        // Transform the edge's labels
        for (ElkLabel elklabel : elkedge.getLabels()) {
            if (!elklabel.getProperty(LayeredOptions.NO_LAYOUT) && !Strings.isNullOrEmpty(elklabel.getText())) {
                LLabel llabel = transformLabel(elklabel);
                ledge.getLabels().add(llabel);
                
                // Depending on the label placement, we want to set graph properties and make sure the
                // edge label placement is actually properly defined
                switch (llabel.getProperty(LayeredOptions.EDGE_LABELS_PLACEMENT)) {
                case HEAD:
                case TAIL:
                    graphProperties.add(GraphProperties.END_LABELS);
                    break;
                    
                case CENTER:
                    graphProperties.add(GraphProperties.CENTER_LABELS);
                    llabel.setProperty(LayeredOptions.EDGE_LABELS_PLACEMENT, EdgeLabelPlacement.CENTER);
                }
            }
        }
        
        // Copy the original bend points of the edge in case they are required
        CrossingMinimizationStrategy crossMinStrat = lgraph.getProperty(LayeredOptions.CROSSING_MINIMIZATION_STRATEGY);
        NodePlacementStrategy nodePlaceStrat = lgraph.getProperty(LayeredOptions.NODE_PLACEMENT_STRATEGY);
        boolean bendPointsRequired = crossMinStrat == CrossingMinimizationStrategy.INTERACTIVE
                || nodePlaceStrat == NodePlacementStrategy.INTERACTIVE;
        
        if (edgeSection != null && !edgeSection.getBendPoints().isEmpty() && bendPointsRequired) {
            KVectorChain originalBendpoints = ElkUtil.createVectorChain(edgeSection);
            KVectorChain importedBendpoints = new KVectorChain();
            
            // MIGRATE We may have to do some coordinate conversion here
            for (KVector point : originalBendpoints) {
                importedBendpoints.add(new KVector(point));
            }
            ledge.setProperty(InternalProperties.ORIGINAL_BENDPOINTS, importedBendpoints);
        }
        
        return ledge;
    }
    
    /**
     * Checks if the given edge has exactly one source and one parent.
     * 
     * @param edge the edge to check.
     * @throws UnsupportedGraphException if the edge does not meet the criteria.
     */
    private void checkEdgeValidity(final ElkEdge edge) {
        if (edge.getSources().isEmpty()) {
            throw new UnsupportedGraphException("Edges must have a source.");
        } else if (edge.getTargets().isEmpty()) {
            throw new UnsupportedGraphException("Edges must have a target.");
        } else if (edge.isHyperedge()) {
            throw new UnsupportedGraphException("Hyperedges are not supported.");
        }
    }
    

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Label Transformation

    /**
     * Transform the given {@code ElkLabel} into an {@code LLabel}.
     * 
     * @param elklabel the label to transform.
     * @return the created {@code LLabel}.
     */
    private LLabel transformLabel(final ElkLabel elklabel) {
        LLabel newLabel = new LLabel(elklabel.getText());
        
        newLabel.copyProperties(elklabel);
        newLabel.setProperty(InternalProperties.ORIGIN, elklabel);
        
        newLabel.getSize().x = elklabel.getWidth();
        newLabel.getSize().y = elklabel.getHeight();
        newLabel.getPosition().x = elklabel.getX();
        newLabel.getPosition().y = elklabel.getY();
        
        return newLabel;
    }

}
