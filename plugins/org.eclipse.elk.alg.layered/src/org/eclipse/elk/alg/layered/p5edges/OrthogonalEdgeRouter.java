/*******************************************************************************
 * Copyright (c) 2010, 2020 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.alg.layered.p5edges;

import java.util.List;
import java.util.ListIterator;
import java.util.Set;

import org.eclipse.elk.alg.layered.LayeredPhases;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LGraphUtil;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.intermediate.IntermediateProcessorStrategy;
import org.eclipse.elk.alg.layered.options.GraphProperties;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.layered.p5edges.orthogonal.OrthogonalRoutingGenerator;
import org.eclipse.elk.alg.layered.p5edges.orthogonal.direction.RoutingDirection;
import org.eclipse.elk.core.alg.ILayoutPhase;
import org.eclipse.elk.core.alg.LayoutProcessorConfiguration;
import org.eclipse.elk.core.util.IElkProgressMonitor;


/**
 * Edge routing implementation that creates orthogonal bend points. Inspired by
 * <ul>
 *   <li>Georg Sander, Layout of directed hypergraphs with orthogonal hyperedges. In
 *     <i>Proceedings of the 11th International Symposium on Graph Drawing (GD '03)</i>,
 *     LNCS vol. 2912, pp. 381-386, Springer, 2004.</li>
 *   <li>Giuseppe di Battista, Peter Eades, Roberto Tamassia, Ioannis G. Tollis,
 *     <i>Graph Drawing: Algorithms for the Visualization of Graphs</i>,
 *     Prentice Hall, New Jersey, 1999 (Section 9.4, for cycle breaking in the
 *     hyperedge segment graph)
 * </ul>
 * 
 * <dl>
 *   <dt>Precondition:</dt><dd>the graph has a proper layering with
 *     assigned node and port positions; the size of each layer is
 *     correctly set; edges connected to ports on strange sides were
 *     processed</dd>
 *   <dt>Postcondition:</dt><dd>each node is assigned a horizontal coordinate;
 *     the bend points of each edge are set; the width of the whole graph is set</dd>
 * </dl>
 */
public final class OrthogonalEdgeRouter implements ILayoutPhase<LayeredPhases, LGraph> {
    
    /* The basic processing strategy for this phase is empty. Depending on the graph features,
     * dependencies on intermediate processors are added dynamically as follows:
     * 
     * Before phase 1:
     *   - None.
     * 
     * Before phase 2:
     *   - For center edge labels:
     *      - LABEL_DUMMY_INSERTER
     * 
     * Before phase 3:
     *   - For non-free ports:
     *     - NORTH_SOUTH_PORT_PREPROCESSOR
     *     - INVERTED_PORT_PROCESSOR
     *   
     *   - For self-loops:
     *     - SELF_LOOP_PROCESSOR
     *   
     *   - For hierarchical ports:
     *     - HIERARCHICAL_PORT_CONSTRAINT_PROCESSOR
     *   
     *   - For center edge labels:
     *     - LABEL_DUMMY_SWITCHER
     * 
     * Before phase 4:
     *   - For hyperedges:
     *     - HYPEREDGE_DUMMY_MERGER
     *   
     *   - For hierarchical ports:
     *     - HIERARCHICAL_PORT_DUMMY_SIZE_PROCESSOR
     *     
     *   - For edge labels:
     *     - LABEL_SIDE_SELECTOR
     *     
     *   - For end edge labels:
     *     - END_LABEL_PREPROCESSOR
     * 
     * Before phase 5:
     * 
     * After phase 5:
     *   - For non-free ports:
     *     - NORTH_SOUTH_PORT_POSTPROCESSOR
     *   
     *   - For hierarchical ports:
     *     - HIERARCHICAL_PORT_ORTHOGONAL_EDGE_ROUTER
     *     
     *   - For center edge labels:
     *     - LABEL_DUMMY_REMOVER
     *     
     *   - For end edge labels:
     *     - END_LABEL_POSTPROCESSOR
     */
    
    /** additional processor dependencies for graphs with hyperedges. */
    private static final LayoutProcessorConfiguration<LayeredPhases, LGraph> HYPEREDGE_PROCESSING_ADDITIONS =
        LayoutProcessorConfiguration.<LayeredPhases, LGraph>create()
            .addBefore(LayeredPhases.P4_NODE_PLACEMENT, IntermediateProcessorStrategy.HYPEREDGE_DUMMY_MERGER);
    
    /** additional processor dependencies for graphs with possible inverted ports. */
    private static final LayoutProcessorConfiguration<LayeredPhases, LGraph> INVERTED_PORT_PROCESSING_ADDITIONS =
        LayoutProcessorConfiguration.<LayeredPhases, LGraph>create()
            .addBefore(LayeredPhases.P3_NODE_ORDERING, IntermediateProcessorStrategy.INVERTED_PORT_PROCESSOR);
    
    /** additional processor dependencies for graphs with northern / southern non-free ports. */
    private static final LayoutProcessorConfiguration<LayeredPhases, LGraph> NORTH_SOUTH_PORT_PROCESSING_ADDITIONS =
        LayoutProcessorConfiguration.<LayeredPhases, LGraph>create()
            .addBefore(LayeredPhases.P3_NODE_ORDERING, IntermediateProcessorStrategy.NORTH_SOUTH_PORT_PREPROCESSOR)
            .addAfter(LayeredPhases.P5_EDGE_ROUTING, IntermediateProcessorStrategy.NORTH_SOUTH_PORT_POSTPROCESSOR);
    
    /** additional processor dependencies for graphs with hierarchical ports. */
    private static final LayoutProcessorConfiguration<LayeredPhases, LGraph> HIERARCHICAL_PORT_PROCESSING_ADDITIONS =
        LayoutProcessorConfiguration.<LayeredPhases, LGraph>create()
            .addBefore(LayeredPhases.P3_NODE_ORDERING,
                    IntermediateProcessorStrategy.HIERARCHICAL_PORT_CONSTRAINT_PROCESSOR)
            .addBefore(LayeredPhases.P4_NODE_PLACEMENT,
                    IntermediateProcessorStrategy.HIERARCHICAL_PORT_DUMMY_SIZE_PROCESSOR)
            .addAfter(LayeredPhases.P5_EDGE_ROUTING,
                    IntermediateProcessorStrategy.HIERARCHICAL_PORT_ORTHOGONAL_EDGE_ROUTER);
    
    /** additional processor dependencies for graphs with self-loops. */
    private static final LayoutProcessorConfiguration<LayeredPhases, LGraph> SELF_LOOP_PROCESSING_ADDITIONS =
            LayoutProcessorConfiguration.<LayeredPhases, LGraph>create()
                .addBefore(LayeredPhases.P1_CYCLE_BREAKING, IntermediateProcessorStrategy.SELF_LOOP_PREPROCESSOR)
                .addAfter(LayeredPhases.P5_EDGE_ROUTING, IntermediateProcessorStrategy.SELF_LOOP_POSTPROCESSOR)
                .before(LayeredPhases.P4_NODE_PLACEMENT)
                    .add(IntermediateProcessorStrategy.SELF_LOOP_PORT_RESTORER)
                    .add(IntermediateProcessorStrategy.SELF_LOOP_ROUTER);
    
    /** additional processor dependencies for graphs with hypernodes. */
    private static final LayoutProcessorConfiguration<LayeredPhases, LGraph> HYPERNODE_PROCESSING_ADDITIONS =
        LayoutProcessorConfiguration.<LayeredPhases, LGraph>create()
            .addAfter(LayeredPhases.P5_EDGE_ROUTING, IntermediateProcessorStrategy.HYPERNODE_PROCESSOR);
    
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
                LayoutProcessorConfiguration.<LayeredPhases, LGraph>create();
        
        // Additional dependencies
        if (graphProperties.contains(GraphProperties.HYPEREDGES)) {
            configuration.addAll(HYPEREDGE_PROCESSING_ADDITIONS);
            configuration.addAll(INVERTED_PORT_PROCESSING_ADDITIONS);
        }
        
        if (graphProperties.contains(GraphProperties.NON_FREE_PORTS)
                || graph.getProperty(LayeredOptions.FEEDBACK_EDGES)) {
            
            configuration.addAll(INVERTED_PORT_PROCESSING_ADDITIONS);

            if (graphProperties.contains(GraphProperties.NORTH_SOUTH_PORTS)) {
                configuration.addAll(NORTH_SOUTH_PORT_PROCESSING_ADDITIONS);
            }
        }

        if (graphProperties.contains(GraphProperties.EXTERNAL_PORTS)) {
            configuration.addAll(HIERARCHICAL_PORT_PROCESSING_ADDITIONS);
        }

        if (graphProperties.contains(GraphProperties.SELF_LOOPS)) {
            configuration.addAll(SELF_LOOP_PROCESSING_ADDITIONS);
        }
        
        if (graphProperties.contains(GraphProperties.HYPERNODES)) {
            configuration.addAll(HYPERNODE_PROCESSING_ADDITIONS);
        }
        
        if (graphProperties.contains(GraphProperties.CENTER_LABELS)) {
            configuration.addAll(CENTER_EDGE_LABEL_PROCESSING_ADDITIONS);
        }
        
        if (graphProperties.contains(GraphProperties.END_LABELS)) {
            configuration.addAll(END_EDGE_LABEL_PROCESSING_ADDITIONS);
        }
        
        return configuration;
    }
    
    @Override
    public void process(final LGraph layeredGraph, final IElkProgressMonitor monitor) {
        monitor.begin("Orthogonal edge routing", 1);
        
        // Retrieve some generic values
        double nodeNodeSpacing =
                layeredGraph.getProperty(LayeredOptions.SPACING_NODE_NODE_BETWEEN_LAYERS).doubleValue();
        double edgeEdgeSpacing =
                layeredGraph.getProperty(LayeredOptions.SPACING_EDGE_EDGE_BETWEEN_LAYERS).doubleValue();
        double edgeNodeSpacing =
                layeredGraph.getProperty(LayeredOptions.SPACING_EDGE_NODE_BETWEEN_LAYERS).doubleValue();
        
        // Prepare for iteration!
        OrthogonalRoutingGenerator routingGenerator = new OrthogonalRoutingGenerator(
                RoutingDirection.WEST_TO_EAST, edgeEdgeSpacing, "phase5");
        float xpos = 0.0f;
        ListIterator<Layer> layerIter = layeredGraph.getLayers().listIterator();
        Layer leftLayer = null;
        Layer rightLayer = null;
        List<LNode> leftLayerNodes = null;
        List<LNode> rightLayerNodes = null;
        int leftLayerIndex = -1;
        int rightLayerIndex = -1;
        
        // Iterate!
        do {
            int slotsCount;
            
            // Fetch the next layer, if any
            rightLayer = layerIter.hasNext() ? layerIter.next() : null;
            rightLayerNodes = rightLayer == null ? null : rightLayer.getNodes();
            rightLayerIndex = layerIter.previousIndex();
            
            // Place the left layer's nodes, if any
            if (leftLayer != null) {
                LGraphUtil.placeNodesHorizontally(leftLayer, xpos);
                xpos += leftLayer.getSize().x;
            }
            
            // Route edges between the two layers
            double startPos = leftLayer == null ? xpos : xpos + edgeNodeSpacing;
            slotsCount = routingGenerator.routeEdges(monitor, layeredGraph, leftLayerNodes, leftLayerIndex,
                    rightLayerNodes, startPos);
            
            boolean isLeftLayerExternal = leftLayer == null || java.util.stream.StreamSupport.stream(leftLayerNodes.spliterator(), false).allMatch(PolylineEdgeRouter.PRED_EXTERNAL_WEST_OR_EAST_PORT);
            boolean isRightLayerExternal = rightLayer == null || java.util.stream.StreamSupport.stream(rightLayerNodes.spliterator(), false).allMatch(PolylineEdgeRouter.PRED_EXTERNAL_WEST_OR_EAST_PORT);
            
            if (slotsCount > 0) {
                // Compute routing area's width
                double routingWidth = (slotsCount - 1) * edgeEdgeSpacing;
                
                if (leftLayer != null) {
                    routingWidth += edgeNodeSpacing;
                }
                
                if (rightLayer != null) {
                    routingWidth += edgeNodeSpacing;
                }
                
                // If we are between two layers, make sure their minimal spacing is preserved
                if (routingWidth < nodeNodeSpacing && !isLeftLayerExternal && !isRightLayerExternal) {
                    routingWidth = nodeNodeSpacing;
                }
                xpos += routingWidth;
            } else if (!isLeftLayerExternal && !isRightLayerExternal) {
                // If all edges are straight, use the usual spacing 
                xpos += nodeNodeSpacing;
            }
            
            leftLayer = rightLayer;
            leftLayerNodes = rightLayerNodes;
            leftLayerIndex = rightLayerIndex;
        } while (rightLayer != null);
        
        layeredGraph.getSize().x = xpos;
        
        monitor.done();
    }
    
}
