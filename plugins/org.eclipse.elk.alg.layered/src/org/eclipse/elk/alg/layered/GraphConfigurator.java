/*******************************************************************************
 * Copyright (c) 2015, 2020 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.alg.layered;

import java.util.Random;
import java.util.Set;

import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LGraphUtil;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.intermediate.IntermediateProcessorStrategy;
import org.eclipse.elk.alg.layered.options.CrossingMinimizationStrategy;
import org.eclipse.elk.alg.layered.options.GraphCompactionStrategy;
import org.eclipse.elk.alg.layered.options.GraphProperties;
import org.eclipse.elk.alg.layered.options.GreedySwitchType;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.layered.options.NodePromotionStrategy;
import org.eclipse.elk.alg.layered.options.OrderingStrategy;
import org.eclipse.elk.alg.layered.options.Spacings;
import org.eclipse.elk.alg.layered.p5edges.EdgeRouterFactory;
import org.eclipse.elk.core.alg.AlgorithmAssembler;
import org.eclipse.elk.core.alg.LayoutProcessorConfiguration;
import org.eclipse.elk.core.labels.LabelManagementOptions;
import org.eclipse.elk.core.options.Direction;
import org.eclipse.elk.core.options.EdgeRouting;
import org.eclipse.elk.core.options.HierarchyHandling;
import org.eclipse.elk.core.options.PortConstraints;

/**
 * The configurator configures a graph in preparation for layout. This includes making sure that a
 * given graph has sensible properties set as well as assembling the phases and processors required
 * to layout the graph. That list is attached to the graph in its
 * {@link InternalProperties#PROCESSORS} property.
 * 
 * <p>Each phase and processor is cached, so a given instance of this class can safely (and usually
 * should) be reused.</p>
 * 
 * @author cds
 */
final class GraphConfigurator {

    ////////////////////////////////////////////////////////////////////////////////
    // Constants
    
    /** intermediate processing configuration for basic graphs. */
    private static final LayoutProcessorConfiguration<LayeredPhases, LGraph> BASELINE_PROCESSING_CONFIGURATION =
        LayoutProcessorConfiguration.<LayeredPhases, LGraph>create()
            .addBefore(LayeredPhases.P4_NODE_PLACEMENT, IntermediateProcessorStrategy.INNERMOST_NODE_MARGIN_CALCULATOR)
            .addBefore(LayeredPhases.P4_NODE_PLACEMENT, IntermediateProcessorStrategy.LABEL_AND_NODE_SIZE_PROCESSOR)
            .addBefore(LayeredPhases.P5_EDGE_ROUTING,
                    IntermediateProcessorStrategy.LAYER_SIZE_AND_GRAPH_HEIGHT_CALCULATOR)
            .addAfter(LayeredPhases.P5_EDGE_ROUTING, IntermediateProcessorStrategy.END_LABEL_SORTER);
    
    /** intermediate processors for label management. */
    private static final LayoutProcessorConfiguration<LayeredPhases, LGraph> LABEL_MANAGEMENT_ADDITIONS =
        LayoutProcessorConfiguration.<LayeredPhases, LGraph>create()
            .addBefore(LayeredPhases.P4_NODE_PLACEMENT, IntermediateProcessorStrategy.CENTER_LABEL_MANAGEMENT_PROCESSOR)
            .addBefore(LayeredPhases.P4_NODE_PLACEMENT,
                       IntermediateProcessorStrategy.END_NODE_PORT_LABEL_MANAGEMENT_PROCESSOR);
    
    /** intermediate processors for hierarchical layout, i.e. {@link HierarchyHandling#INCLUDE_CHILDREN}. */
    private static final LayoutProcessorConfiguration<LayeredPhases, LGraph> HIERARCHICAL_ADDITIONS = 
        LayoutProcessorConfiguration.<LayeredPhases, LGraph>create()
            .addAfter(LayeredPhases.P5_EDGE_ROUTING, IntermediateProcessorStrategy.HIERARCHICAL_NODE_RESIZER);
    
    ////////////////////////////////////////////////////////////////////////////////
    // Processor Caching
    
    /** The algorithm assembler we use to assemble our algorithm configurations. */
    private final AlgorithmAssembler<LayeredPhases, LGraph> algorithmAssembler =
            AlgorithmAssembler.<LayeredPhases, LGraph>create(LayeredPhases.class);
    
    
    ////////////////////////////////////////////////////////////////////////////////
    // Graph Preprocessing (Property Configuration)

    /** the minimal spacing between edges, so edges won't overlap. */
    private static final double MIN_EDGE_SPACING = 2.0;

    /**
     * Set special layout options for the layered graph.
     * 
     * @param lgraph a new layered graph
     */
    private void configureGraphProperties(final LGraph lgraph) {
        // check the bounds of some layout options
        // TODO Find a new concept for checking validity of bounds
//        lgraph.checkProperties(InternalProperties.SPACING, InternalProperties.BORDER_SPACING,
//                Properties.THOROUGHNESS, InternalProperties.ASPECT_RATIO);
        
        double edgeSpacing = lgraph.getProperty(LayeredOptions.SPACING_EDGE_EDGE);
        if (edgeSpacing < MIN_EDGE_SPACING) {
            // Make sure the resulting edge spacing is at least 2 in order to avoid overlapping edges.
            lgraph.setProperty(LayeredOptions.SPACING_EDGE_EDGE, MIN_EDGE_SPACING);
        }
        
        Direction direction = lgraph.getProperty(LayeredOptions.DIRECTION);
        if (direction == Direction.UNDEFINED) {
            lgraph.setProperty(LayeredOptions.DIRECTION, LGraphUtil.getDirection(lgraph));
        }
        
        // set the random number generator based on the random seed option
        Integer randomSeed = lgraph.getProperty(LayeredOptions.RANDOM_SEED);
        if (randomSeed == 0) {
            lgraph.setProperty(InternalProperties.RANDOM, new Random());
        } else {
            lgraph.setProperty(InternalProperties.RANDOM, new Random(randomSeed));
        }
        
        Boolean favorStraightness = lgraph.getProperty(LayeredOptions.NODE_PLACEMENT_FAVOR_STRAIGHT_EDGES);
        if (favorStraightness == null) {
            lgraph.setProperty(LayeredOptions.NODE_PLACEMENT_FAVOR_STRAIGHT_EDGES,
                    lgraph.getProperty(LayeredOptions.EDGE_ROUTING) == EdgeRouting.ORTHOGONAL);
        }
        
        // copy the port constraints to keep a list of original port constraints
        copyPortContraints(lgraph);
        
        // pre-calculate spacing information
        Spacings spacings = new Spacings(lgraph);
        lgraph.setProperty(InternalProperties.SPACINGS, spacings);
    }
    
    private void copyPortContraints(final LGraph lgraph) {
        lgraph.getLayerlessNodes().stream()
            .forEach(lnode -> copyPortConstraints(lnode));
        lgraph.getLayers().stream()
            .flatMap(layer -> layer.getNodes().stream())
            .forEach(lnode -> copyPortConstraints(lnode));
    }
    
    private void copyPortConstraints(final LNode node) {
        PortConstraints originalPortconstraints = node.getProperty(LayeredOptions.PORT_CONSTRAINTS);
        node.setProperty(InternalProperties.ORIGINAL_PORT_CONSTRAINTS, originalPortconstraints);

        LGraph nestedGraph = node.getNestedGraph();
        if (nestedGraph != null) {
            copyPortContraints(nestedGraph);
        }
    }
    
    ////////////////////////////////////////////////////////////////////////////////
    // Updating the Configuration
    
    /**
     * Rebuilds the configuration to include all processors required to layout the given graph. The list
     * of processors is attached to the graph in the {@link InternalProperties#PROCESSORS} property.
     * 
     * @param lgraph the graph to layout.
     */
    public void prepareGraphForLayout(final LGraph lgraph) {
        // Make sure the graph properties are sensible
        configureGraphProperties(lgraph);
        
        // Setup the algorithm assembler
        algorithmAssembler.reset();
        
        algorithmAssembler.setPhase(LayeredPhases.P1_CYCLE_BREAKING,
                lgraph.getProperty(LayeredOptions.CYCLE_BREAKING_STRATEGY));
        algorithmAssembler.setPhase(LayeredPhases.P2_LAYERING,
                lgraph.getProperty(LayeredOptions.LAYERING_STRATEGY));
        algorithmAssembler.setPhase(LayeredPhases.P3_NODE_ORDERING,
                lgraph.getProperty(LayeredOptions.CROSSING_MINIMIZATION_STRATEGY));
        algorithmAssembler.setPhase(LayeredPhases.P4_NODE_PLACEMENT,
                lgraph.getProperty(LayeredOptions.NODE_PLACEMENT_STRATEGY));
        algorithmAssembler.setPhase(LayeredPhases.P5_EDGE_ROUTING,
                EdgeRouterFactory.factoryFor(lgraph.getProperty(LayeredOptions.EDGE_ROUTING)));
        
        algorithmAssembler.addProcessorConfiguration(getPhaseIndependentLayoutProcessorConfiguration(lgraph));
        
        lgraph.setProperty(InternalProperties.PROCESSORS, algorithmAssembler.build(lgraph));
    }
    
    /**
     * Returns an intermediate processing configuration with processors not tied to specific phases.
     * 
     * @param lgraph the layered graph to be processed. The configuration may vary depending on certain
     *               properties of the graph.
     * @return intermediate processing configuration. May be {@code null}.
     */
    private LayoutProcessorConfiguration<LayeredPhases, LGraph> getPhaseIndependentLayoutProcessorConfiguration(
            final LGraph lgraph) {

        Set<GraphProperties> graphProperties = lgraph.getProperty(InternalProperties.GRAPH_PROPERTIES);

        // Basic configuration
        LayoutProcessorConfiguration<LayeredPhases, LGraph> configuration =
                LayoutProcessorConfiguration.createFrom(BASELINE_PROCESSING_CONFIGURATION);
        
        // Hierarchical layout. 
        //  Note that the recursive graph layout engine made sure that at this point
        //  'INCLUDE_CHILDREN' has been propagated to all parent nodes with 'INHERIT'.
        //  Thus, every lgraph of the hierarchical lgraph structure is configured with the following additions. 
        HierarchyHandling hierarchyHandling = lgraph.getProperty(LayeredOptions.HIERARCHY_HANDLING);
        if (hierarchyHandling == HierarchyHandling.INCLUDE_CHILDREN) {
            configuration.addAll(HIERARCHICAL_ADDITIONS);
        }

        // Port side processor, put to first slot only if requested and routing is orthogonal
        if (lgraph.getProperty(LayeredOptions.FEEDBACK_EDGES)) {
            configuration.addBefore(LayeredPhases.P1_CYCLE_BREAKING, IntermediateProcessorStrategy.PORT_SIDE_PROCESSOR);
        } else {
            configuration.addBefore(LayeredPhases.P3_NODE_ORDERING, IntermediateProcessorStrategy.PORT_SIDE_PROCESSOR);
        }
        
        // If the graph has a label manager, so add label management additions
        if (lgraph.getProperty(LabelManagementOptions.LABEL_MANAGER) != null) {
            configuration.addAll(LABEL_MANAGEMENT_ADDITIONS);
        }
        
        // If the graph should be laid out interactively or you really want it,
        // add the layers and positions to the nodes.
        if (lgraph.getProperty(LayeredOptions.INTERACTIVE_LAYOUT)
                || lgraph.getProperty(LayeredOptions.GENERATE_POSITION_AND_LAYER_IDS)) {
            configuration.addAfter(LayeredPhases.P5_EDGE_ROUTING, 
                    IntermediateProcessorStrategy.CONSTRAINTS_POSTPROCESSOR);
        }

        // graph transformations for unusual layout directions
        switch (lgraph.getProperty(LayeredOptions.DIRECTION)) {
        case LEFT:
        case DOWN:
        case UP:
            configuration
                .addBefore(LayeredPhases.P1_CYCLE_BREAKING, IntermediateProcessorStrategy.DIRECTION_PREPROCESSOR)
                .addAfter(LayeredPhases.P5_EDGE_ROUTING, IntermediateProcessorStrategy.DIRECTION_POSTPROCESSOR);
            break;
        default:
            // This is either RIGHT or UNDEFINED, which is just mapped to RIGHT. Either way, we
            // don't need any processors here
            break;
        }

        // Additional dependencies
        if (graphProperties.contains(GraphProperties.COMMENTS)) {
            configuration
                .addBefore(LayeredPhases.P1_CYCLE_BREAKING, IntermediateProcessorStrategy.COMMENT_PREPROCESSOR)
                .addBefore(LayeredPhases.P4_NODE_PLACEMENT,
                        IntermediateProcessorStrategy.COMMENT_NODE_MARGIN_CALCULATOR)
                .addAfter(LayeredPhases.P5_EDGE_ROUTING, IntermediateProcessorStrategy.COMMENT_POSTPROCESSOR);
        }
        
        // Node-Promotion application for reduction of dummy nodes after layering
        if (lgraph.getProperty(LayeredOptions.LAYERING_NODE_PROMOTION_STRATEGY) != NodePromotionStrategy.NONE) {
            configuration.addBefore(LayeredPhases.P3_NODE_ORDERING, IntermediateProcessorStrategy.NODE_PROMOTION);
        }

        // Preserve certain partitions during layering
        if (graphProperties.contains(GraphProperties.PARTITIONS)) {
            configuration.addBefore(LayeredPhases.P1_CYCLE_BREAKING,
                    IntermediateProcessorStrategy.PARTITION_PREPROCESSOR);
            configuration.addBefore(LayeredPhases.P2_LAYERING,
                    IntermediateProcessorStrategy.PARTITION_MIDPROCESSOR);
            configuration.addBefore(LayeredPhases.P3_NODE_ORDERING,
                    IntermediateProcessorStrategy.PARTITION_POSTPROCESSOR);
        }
        
        // Additional horizontal compaction depends on orthogonal edge routing
        if (lgraph.getProperty(LayeredOptions.COMPACTION_POST_COMPACTION_STRATEGY) != GraphCompactionStrategy.NONE
              && lgraph.getProperty(LayeredOptions.EDGE_ROUTING) != EdgeRouting.POLYLINE) {
            configuration.addAfter(LayeredPhases.P5_EDGE_ROUTING, IntermediateProcessorStrategy.HORIZONTAL_COMPACTOR);
        }
        
        // Move trees of high degree nodes to separate layers
        if (lgraph.getProperty(LayeredOptions.HIGH_DEGREE_NODES_TREATMENT)) {
            configuration.addBefore(LayeredPhases.P3_NODE_ORDERING,
                    IntermediateProcessorStrategy.HIGH_DEGREE_NODE_LAYER_PROCESSOR);
        }
        
        // Introduce in-layer constraints to preserve the order of regular nodes
        if (lgraph.getProperty(LayeredOptions.CROSSING_MINIMIZATION_SEMI_INTERACTIVE)) {
            configuration.addBefore(LayeredPhases.P3_NODE_ORDERING,
                    IntermediateProcessorStrategy.SEMI_INTERACTIVE_CROSSMIN_PROCESSOR);
        }

        // Configure greedy switch
        //  Note that in the case of hierarchical layout, the configuration may further be adjusted by
        //  ElkLayered#reviewAndCorrectHierarchicalProcessors(...)
        if (activateGreedySwitchFor(lgraph)) {
            final GreedySwitchType greedySwitchType;
            if (isHierarchicalLayout(lgraph)) {
                greedySwitchType =
                        lgraph.getProperty(LayeredOptions.CROSSING_MINIMIZATION_GREEDY_SWITCH_HIERARCHICAL_TYPE);
            } else {
                greedySwitchType = lgraph.getProperty(LayeredOptions.CROSSING_MINIMIZATION_GREEDY_SWITCH_TYPE);
            }
            IntermediateProcessorStrategy internalGreedyType = (greedySwitchType == GreedySwitchType.ONE_SIDED)
                    ? IntermediateProcessorStrategy.ONE_SIDED_GREEDY_SWITCH
                    : IntermediateProcessorStrategy.TWO_SIDED_GREEDY_SWITCH;
            configuration.addBefore(LayeredPhases.P4_NODE_PLACEMENT, internalGreedyType);
        }
        
        switch (lgraph.getProperty(LayeredOptions.LAYER_UNZIPPING_STRATEGY)) {
        case ALTERNATING:
            configuration.addBefore(LayeredPhases.P4_NODE_PLACEMENT, IntermediateProcessorStrategy.ALTERNATING_LAYER_UNZIPPER);
            break;
        default:
            break;
        }

        // Wrapping of graphs
        switch (lgraph.getProperty(LayeredOptions.WRAPPING_STRATEGY)) {
        case SINGLE_EDGE:
            configuration.addBefore(LayeredPhases.P4_NODE_PLACEMENT,
                    IntermediateProcessorStrategy.SINGLE_EDGE_GRAPH_WRAPPER);
            break;
        case MULTI_EDGE:
            configuration
                    .addBefore(LayeredPhases.P3_NODE_ORDERING,
                            IntermediateProcessorStrategy.BREAKING_POINT_INSERTER)
                    .addBefore(LayeredPhases.P4_NODE_PLACEMENT, IntermediateProcessorStrategy.BREAKING_POINT_PROCESSOR)
                    .addAfter(LayeredPhases.P5_EDGE_ROUTING, IntermediateProcessorStrategy.BREAKING_POINT_REMOVER);
            break;
        default: // OFF
        }
        
        if (lgraph.getProperty(LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY) != OrderingStrategy.NONE) {
            configuration
                .addBefore(LayeredPhases.P3_NODE_ORDERING, IntermediateProcessorStrategy.SORT_BY_INPUT_ORDER_OF_MODEL);
        }
        
        return configuration;
    }
    
    /**
     * Greedy switch may be activated if the following holds.
     * 
     * <h3>Hierarchical layout</h3>
     * <ol>
     *  <li>the {@link LayeredOptions#CROSSING_MINIMIZATION_GREEDY_SWITCH_HIERARCHICAL_TYPE} is set to something 
     *      different than OFF</li>
     * </ol>
     * 
     * <h3>Non-hierarchical layout</h3>
     * <ol>
     *  <li>no interactive crossing minimization is performed</li>
     *  <li>the {@link LayeredOptions#CROSSING_MINIMIZATION_GREEDY_SWITCH_TYPE} option is set to something different 
     *      than OFF</li>
     *  <li>the activationThreshold is larger than or equal to the graph's number of nodes (or '0')</li>
     * </ol>
     * @return {@code true} if above condition holds, {@code false} otherwise.
     */
    public static boolean activateGreedySwitchFor(final LGraph lgraph) {
        
        // First, check the hierarchical case
        if (isHierarchicalLayout(lgraph)) {
            // Note that we only activate it for the root in the hierarchical case and rely on
            // ElkLayered#reviewAndCorrectHierarchicalProcessors(...) to properly set it for the whole hierarchy.
            // Also, interactive crossing minimization is not allowed/applicable for hierarchical layout
            // (and thus doesn't have to be checked here).
            return lgraph.getParentNode() == null && lgraph.getProperty(
                    LayeredOptions.CROSSING_MINIMIZATION_GREEDY_SWITCH_HIERARCHICAL_TYPE) != GreedySwitchType.OFF;
        }

        // Second, if not hierarchical, check the slightly more complex non-hierarchical case
        GreedySwitchType greedySwitchType = lgraph.getProperty(LayeredOptions.CROSSING_MINIMIZATION_GREEDY_SWITCH_TYPE);
        boolean interactiveCrossMin =
                lgraph.getProperty(LayeredOptions.CROSSING_MINIMIZATION_SEMI_INTERACTIVE) 
                || lgraph.getProperty(
                        LayeredOptions.CROSSING_MINIMIZATION_STRATEGY) == CrossingMinimizationStrategy.INTERACTIVE;
        int activationThreshold =
                lgraph.getProperty(LayeredOptions.CROSSING_MINIMIZATION_GREEDY_SWITCH_ACTIVATION_THRESHOLD);
        int graphSize = lgraph.getLayerlessNodes().size();
        
        return !interactiveCrossMin 
                && greedySwitchType != GreedySwitchType.OFF 
                && (activationThreshold == 0 || activationThreshold > graphSize);
    }

    private static boolean isHierarchicalLayout(final LGraph lgraph) {
        return lgraph.getProperty(LayeredOptions.HIERARCHY_HANDLING) == HierarchyHandling.INCLUDE_CHILDREN;
    }

}
