/*******************************************************************************
 * Copyright (c) 2011, 2015 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.alg.layered.components;

import java.util.Arrays;
import java.util.Collections;
import java.util.EnumSet;
import java.util.List;
import java.util.Set;

import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LGraphUtil;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LNode.NodeType;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.options.GraphProperties;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.core.options.PortConstraints;
import org.eclipse.elk.core.options.PortSide;
import org.eclipse.elk.core.util.Pair;

import com.google.common.collect.Lists;

/**
 * A processor that is able to split an input graph into connected components and to pack those
 * components after layout.
 * 
 * <p>If the graph has no external ports, splitting it into components is straightforward and always
 * works. If on the other hand it does have external ports, splitting the graph into connected
 * components is problematic because the port positions might introduce constraints on the placement
 * of the different components. More or less simple solutions have only been implemented for the cases
 * of port constraints set to {@link org.eclipse.elk.core.options.PortConstraints#FREE FREE} or
 * {@link org.eclipse.elk.core.options.PortConstraints#FIXED_SIDE FIXED_SIDE}. If the graph contains
 * external ports with port constraints other than these, connected components processing is disabled
 * even if requested by the user.</p>
 * 
 * <h4>Splitting into components</h4>
 * <dl>
 *   <dt>Precondition:</dt>
 *     <dd>an unlayered graph.</dd>
 *   <dt>Postcondition:</dt>
 *     <dd>a list of graphs that represent the connected components of the input graph.</dd>
 * </dl>
 * 
 * <h4>Packing components</h4>
 * <dl>
 *   <dt>Precondition:</dt>
 *     <dd>a list of unlayered graphs with complete layout.</dd>
 *   <dt>Postcondition:</dt>
 *     <dd>a single unlayered graph.</dd>
 * </dl>
 *
 * @author msp
 * @author cds
 */
public final class ComponentsProcessor {
    
    /** Cached instance of a {@link ComponentGroupGraphPlacer}. */
    private final ComponentGroupGraphPlacer componentGroupGraphPlacer = new ComponentGroupGraphPlacer();
    /** Cached instance of a {@link ComponentGroupGraphPlacer}. */
    private final ComponentGroupModelOrderGraphPlacer componentGroupModelOrderGraphPlacer =
            new ComponentGroupModelOrderGraphPlacer();
    /** Cached instance of a {@link ComponentGroupGraphPlacer}. */
    private final ModelOrderRowGraphPlacer modelOrderRowGraphPlacer =
            new ModelOrderRowGraphPlacer();
    /** Cached instance of a {@link SimpleRowGraphPlacer}. */
    private final SimpleRowGraphPlacer simpleRowGraphPlacer = new SimpleRowGraphPlacer();
    /** Graph placer to be used to combine the different components back into a single graph. */
    private AbstractGraphPlacer graphPlacer;
    

    /**
     * Split the given graph into its connected components.
     * 
     * @param graph an input graph with layerless nodes
     * @return a list of components that can be processed one by one
     */
    public List<LGraph> split(final LGraph graph) {
        List<LGraph> result;
        
        // Default to the simple graph placer
        graphPlacer = simpleRowGraphPlacer;
        
        // Whether separate components processing is requested
        Boolean separateProperty = graph.getProperty(LayeredOptions.SEPARATE_CONNECTED_COMPONENTS);
        boolean separate = separateProperty == null || separateProperty.booleanValue();
        
        // Whether the graph contains external ports
        boolean extPorts = graph.getProperty(InternalProperties.GRAPH_PROPERTIES).contains(
                GraphProperties.EXTERNAL_PORTS);
        
        // The graph's external port constraints
        PortConstraints extPortConstraints = graph.getProperty(LayeredOptions.PORT_CONSTRAINTS);
        boolean compatiblePortConstraints = !extPortConstraints.isOrderFixed();
        
        // The graph may only be separated 
        //  1. Separation was requested.
        //  2. If the graph contains external ports, port constraints are set to either
        //     FREE or FIXED_SIDES.
        if (separate && (compatiblePortConstraints || !extPorts)) {
            // Set id of all nodes to 0
            for (LNode node : graph.getLayerlessNodes()) {
                node.id = 0;
            }
            
            // Perform DFS starting on each node, collecting connected components
            result = Lists.newArrayList();
            for (LNode node : graph.getLayerlessNodes()) {
                Pair<List<LNode>, Set<PortSide>> componentData = dfs(node, null);
                
                if (componentData != null) {
                    LGraph newGraph = new LGraph();
                    
                    newGraph.copyProperties(graph);
                    newGraph.setProperty(InternalProperties.EXT_PORT_CONNECTIONS,
                            componentData.getSecond());
                    newGraph.getPadding().copy(graph.getPadding());
                    
                    // If a minimum size was set on the original graph, setting it on the seperated graphs as well
                    // might enlarge them although their combined area might not actually have fallen below the minimum
                    // size; thus, remove the minimum size
                    newGraph.setProperty(LayeredOptions.NODE_SIZE_MINIMUM, null);
                    
                    for (LNode n : componentData.getFirst()) {
                        newGraph.getLayerlessNodes().add(n);
                        n.setGraph(newGraph);
                    }
                    
                    result.add(newGraph);
                }
            }
            
            if (extPorts) {
                // With external port connections, we want to use the more complex components
                // placement algorithm
                if (graph.getProperty(LayeredOptions.CONSIDER_MODEL_ORDER_COMPONENTS)
                        == ComponentOrderingStrategy.GROUP_MODEL_ORDER) {
                    graphPlacer = componentGroupModelOrderGraphPlacer;
                } else if (graph.getProperty(LayeredOptions.CONSIDER_MODEL_ORDER_COMPONENTS)
                        == ComponentOrderingStrategy.MODEL_ORDER) {
                    graphPlacer = modelOrderRowGraphPlacer;
                } else {
                    graphPlacer = componentGroupGraphPlacer;
                }
            }
        } else {
            result = Arrays.asList(graph);
        }
        // FIXME order by the dominant group instead. Needs property dominantComponentOrderingGroup
        // If model order should be preserved the connected components should be ordered by their elements.
        // The component with the node with the smallest order should be first.
        if (graph.getProperty(LayeredOptions.CONSIDER_MODEL_ORDER_COMPONENTS) != ComponentOrderingStrategy.NONE) {
            Collections.sort(result, (g1, g2) -> {
                int g1Order = LGraphUtil.getMinimalModelOrder(g1);
                int g2Order = LGraphUtil.getMinimalModelOrder(g2);
                return Integer.compare(g1Order, g2Order);
            });
        }
        
        return result;
    }
    
    /**
     * Perform a DFS starting on the given node, collect all nodes that are found in the corresponding
     * connected component and return the set of external port sides the component connects to.
     * 
     * @param node a node.
     * @param data pair of nodes in the component and external port sides used to produce the result
     *             during recursive calls. Should be {@code null} when this method is called.
     * @return a pairing of the connected component and the set of port sides of external ports it
     *         connects to, or {@code null} if the node was already visited
     */
    private Pair<List<LNode>, Set<PortSide>> dfs(final LNode node,
            final Pair<List<LNode>, Set<PortSide>> data) {
        
        if (node.id == 0) {
            // Mark the node as visited
            node.id = 1;
            
            // Check if we already have a list of nodes for the connected component
            Pair<List<LNode>, Set<PortSide>> mutableData = data;
            if (mutableData == null) {
                List<LNode> component = Lists.newArrayList();
                Set<PortSide> extPortSides = EnumSet.noneOf(PortSide.class);
                
                mutableData = new Pair<List<LNode>, Set<PortSide>>(component, extPortSides);
            }
            
            // Add this node to the component
            mutableData.getFirst().add(node);
            
            // Check if this node is an external port dummy and, if so, add its side
            if (node.getType() == NodeType.EXTERNAL_PORT) {
                mutableData.getSecond().add(node.getProperty(InternalProperties.EXT_PORT_SIDE));
            }
            
            // DFS
            for (LPort port1 : node.getPorts()) {
                for (LPort port2 : port1.getConnectedPorts()) {
                    dfs(port2.getNode(), mutableData);
                }
            }
            
            return mutableData;
        }
        
        // The node was already visited
        return null;
    }
    
    /**
     * Combine the given components into a single graph by moving them around such that they are
     * placed next and beneath to each other instead of overlapping.
     * 
     * @param components a list of components
     * @param target the target graph into which the others shall be combined
     */
    public void combine(final List<LGraph> components, final LGraph target) {
        graphPlacer.combine(components, target);
    }
    
}
