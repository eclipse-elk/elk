/*******************************************************************************
 * Copyright (c) 2010, 2020 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.alg.layered.intermediate;

import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.intermediate.compaction.HorizontalGraphCompactor;
import org.eclipse.elk.alg.layered.intermediate.unzipping.AlternatingLayerUnzipper;
import org.eclipse.elk.alg.layered.intermediate.wrapping.BreakingPointInserter;
import org.eclipse.elk.alg.layered.intermediate.wrapping.BreakingPointProcessor;
import org.eclipse.elk.alg.layered.intermediate.wrapping.BreakingPointRemover;
import org.eclipse.elk.alg.layered.intermediate.wrapping.SingleEdgeGraphWrapper;
import org.eclipse.elk.alg.layered.p3order.LayerSweepCrossingMinimizer;
import org.eclipse.elk.alg.layered.p3order.LayerSweepCrossingMinimizer.CrossMinType;
import org.eclipse.elk.core.alg.ILayoutProcessor;
import org.eclipse.elk.core.alg.ILayoutProcessorFactory;

/**
 * Definition of available intermediate layout processors for the layered layouter. This enumeration also serves as a
 * factory for intermediate layout processors.
 */
public enum IntermediateProcessorStrategy implements ILayoutProcessorFactory<LGraph> {

    /*
     * In this enumeration, intermediate layout processors are listed by the earliest slot in which
     * they can sensibly be used. The order in which they are listed is determined by the
     * dependencies on other processors.
     */

    // Before Phase 1

    /** Transforms the graph to internally perform a left-to-right drawing. */
    DIRECTION_PREPROCESSOR,
    /** Removes some comment boxes to place them separately in a post-processor. */
    COMMENT_PREPROCESSOR,
    /** Makes sure nodes with layer constraints have only incoming or only outgoing edges. */
    EDGE_AND_LAYER_CONSTRAINT_EDGE_REVERSER,
    /** If one of the phases is set to interactive mode, this processor positions external ports. */
    INTERACTIVE_EXTERNAL_PORT_POSITIONER,
    /** Reverse edges that run from higher-index to lower-index partitions. */
    PARTITION_PREPROCESSOR,
    
    // Before Phase 2

    /** Adds dummy nodes in edges where center labels are present. */
    LABEL_DUMMY_INSERTER,
    /** Takes care of self loop preprocessing. */
    SELF_LOOP_PREPROCESSOR,
    /** Hides {@code FIRST_SEPARATE} and {@code LAST_SEPARATE} nodes. */
    LAYER_CONSTRAINT_PREPROCESSOR,
    /** Adds edges that force layering algorithms to respect layout partitions. */
    PARTITION_MIDPROCESSOR,

    // Before Phase 3
    
    /** Moves trees of high degree nodes to separate layers. */
    HIGH_DEGREE_NODE_LAYER_PROCESSOR,
    /** Node-promotion for prettier graphs, especially algorithms like longest-path are prettified. */
    NODE_PROMOTION,
    /** Restores {@code FIRST_SEPARATE} and {@code LAST_SEPARATE} nodes and enforces layer constraints. */
    LAYER_CONSTRAINT_POSTPROCESSOR,
    /** Remove partition constraint edges. */
    PARTITION_POSTPROCESSOR,
    /** Handles northern and southern hierarchical ports. */
    HIERARCHICAL_PORT_CONSTRAINT_PROCESSOR,
    /** Adds successor constraints between regular nodes before crossing minimization. */
    SEMI_INTERACTIVE_CROSSMIN_PROCESSOR,
    /** Inserts breaking points which are used to 'wrap' the graph after crossing minimization. */
    BREAKING_POINT_INSERTER,
    /** Takes a layered graph and turns it into a properly layered graph. */
    LONG_EDGE_SPLITTER,
    /** Makes sure nodes have at least fixed port sides. */
    PORT_SIDE_PROCESSOR,
    /** Takes a layered graph and inserts dummy nodes for edges connected to inverted ports. */
    INVERTED_PORT_PROCESSOR,
    /** Orders the port lists of nodes with fixed port order. */
    PORT_LIST_SORTER,
    /** Sorts nodes and edges according to their order in the input model. */
    SORT_BY_INPUT_ORDER_OF_MODEL,
    /** Inserts dummy nodes to take care of northern and southern ports. */
    NORTH_SOUTH_PORT_PREPROCESSOR,

    // Before Phase 4
    
    /** Performs 'wrapping' of the graph, potentially executing improvement heuristics. */
    BREAKING_POINT_PROCESSOR,
    /** Hierarchical one-sided greedy switch crossing reduction. */
    ONE_SIDED_GREEDY_SWITCH,
    /** Hierarchical two-sided greedy switch crossing reduction. */
    TWO_SIDED_GREEDY_SWITCH,
    /** Position self loops after phase 3. */
    SELF_LOOP_PORT_RESTORER,
    /** Unzips layers for compaction. */
    ALTERNATING_LAYER_UNZIPPER,
    /** Wraps graphs such that they better fit a given drawing area, allowing only a single edge per cut. */
    SINGLE_EDGE_GRAPH_WRAPPER,
    /** Makes sure that in-layer constraints are handled. */
    IN_LAYER_CONSTRAINT_PROCESSOR,
    /** Manages edge end labels, node labels, and port labels. */
    END_NODE_PORT_LABEL_MANAGEMENT_PROCESSOR,
    /** Sets the positions of ports and labels, and sets the node sizes. */
    LABEL_AND_NODE_SIZE_PROCESSOR,
    /** Calculates the margins of nodes according to the sizes of ports and port labels. */
    INNERMOST_NODE_MARGIN_CALCULATOR,
    /** Calculate the bendpoints for the self-loop edges. */
    SELF_LOOP_ROUTER,
    /** Extends node margin by the space required for comment boxes. */
    COMMENT_NODE_MARGIN_CALCULATOR,
    /** Place end labels of edges and extend node margins accordingly. */
    END_LABEL_PREPROCESSOR,
    /** Tries to switch the label dummy nodes which the middle most dummy node of a long edge. */
    LABEL_DUMMY_SWITCHER,
    /** Tries to shorten edge center labels where necessary. */
    CENTER_LABEL_MANAGEMENT_PROCESSOR,
    /** Decides, on which side of an edge the edge labels should be placed. */
    LABEL_SIDE_SELECTOR,
    /** Merges long edge dummy nodes belonging to the same hyperedge. */
    HYPEREDGE_DUMMY_MERGER,
    /** Adjusts the width of hierarchical port dummy nodes. */
    HIERARCHICAL_PORT_DUMMY_SIZE_PROCESSOR,

    // Before Phase 5

    /** Calculate the size of layers and the graph's height and offset. */
    LAYER_SIZE_AND_GRAPH_HEIGHT_CALCULATOR,
    /** Fix coordinates of hierarchical port dummy nodes. */
    HIERARCHICAL_PORT_POSITION_PROCESSOR,

    // After Phase 5
    
    /** Adds the layer and positions that were computed to the nodes. */
    CONSTRAINTS_POSTPROCESSOR,
    /** Reinserts and places comment boxes that have been removed before. */
    COMMENT_POSTPROCESSOR,
    /** Moves hypernodes horizontally for better placement. */
    HYPERNODE_PROCESSOR,
    /** Routes edges incident to hierarchical ports orthogonally. */
    HIERARCHICAL_PORT_ORTHOGONAL_EDGE_ROUTER,
    /** Takes a properly layered graph and removes the dummy nodes due to proper layering. */
    LONG_EDGE_JOINER,
    /** Add node position to the already routed self-loops. */
    SELF_LOOP_POSTPROCESSOR,
    /** Removes the breaking points that were inserted for 'wrapping', derives edge routes correspondingly. */
    BREAKING_POINT_REMOVER,
    /** Removes dummy nodes inserted by the north south side preprocessor and routes edges. */
    NORTH_SOUTH_PORT_POSTPROCESSOR,
    /** Moves nodes and vertical edge segments in horizontal direction to close some gaps that are a
      * result of the layering. */
    HORIZONTAL_COMPACTOR,
    /** Removes dummy nodes which were introduced for center labels. */
    LABEL_DUMMY_REMOVER,
    /** Processor that calculates final control points for splines as intended by the 
     * {@link org.eclipse.elk.alg.layered.p5edges.splines.SplineEdgeRouter SplineEdgeRouter}. */
    FINAL_SPLINE_BENDPOINTS_CALCULATOR,
    /** Sorts multiple end labels to reflect the order of source or target nodes. */
    END_LABEL_SORTER,
    /** Takes the reversed edges of a graph and restores their original direction. */
    REVERSED_EDGE_RESTORER,
    /** Assign correct coordinates to end labels. */
    END_LABEL_POSTPROCESSOR,
    /** In hierarchical graphs, maps a child graph to its parent node. */
    HIERARCHICAL_NODE_RESIZER,
    /** Transforms the internal left-to-right drawing back to the original direction. */
    DIRECTION_POSTPROCESSOR;

    /**
     * Creates an instance of the layout processor described by this instance.
     * 
     * @return the layout processor.
     */
    // SUPPRESS CHECKSTYLE NEXT MethodLength
    public ILayoutProcessor<LGraph> create() {
        switch (this) {

        case BREAKING_POINT_INSERTER:
            return new BreakingPointInserter();
            
        case BREAKING_POINT_PROCESSOR:
            return new BreakingPointProcessor();
            
        case BREAKING_POINT_REMOVER:
            return new BreakingPointRemover();
            
        case CENTER_LABEL_MANAGEMENT_PROCESSOR:
            return new LabelManagementProcessor(true);
            
        case COMMENT_NODE_MARGIN_CALCULATOR:
            return new CommentNodeMarginCalculator();
            
        case COMMENT_POSTPROCESSOR:
            return new CommentPostprocessor();
            
        case COMMENT_PREPROCESSOR:
            return new CommentPreprocessor();
            
        case CONSTRAINTS_POSTPROCESSOR:
            return new ConstraintsPostprocessor();

        case DIRECTION_POSTPROCESSOR:
            return new GraphTransformer(GraphTransformer.Mode.TO_INTERNAL_LTR);
            
        case DIRECTION_PREPROCESSOR:
            return new GraphTransformer(GraphTransformer.Mode.TO_INPUT_DIRECTION);

        case EDGE_AND_LAYER_CONSTRAINT_EDGE_REVERSER:
            return new EdgeAndLayerConstraintEdgeReverser();
            
        case END_LABEL_POSTPROCESSOR:
            return new EndLabelPostprocessor();
            
        case END_LABEL_PREPROCESSOR:
            return new EndLabelPreprocessor();
        
        case END_NODE_PORT_LABEL_MANAGEMENT_PROCESSOR:
            return new LabelManagementProcessor(false);
            
        case FINAL_SPLINE_BENDPOINTS_CALCULATOR:
            return new FinalSplineBendpointsCalculator();

        case HIERARCHICAL_NODE_RESIZER:
            return new HierarchicalNodeResizingProcessor();

        case HIERARCHICAL_PORT_CONSTRAINT_PROCESSOR:
            return new HierarchicalPortConstraintProcessor();

        case HIERARCHICAL_PORT_DUMMY_SIZE_PROCESSOR:
            return new HierarchicalPortDummySizeProcessor();

        case HIERARCHICAL_PORT_ORTHOGONAL_EDGE_ROUTER:
            return new HierarchicalPortOrthogonalEdgeRouter();

        case HIERARCHICAL_PORT_POSITION_PROCESSOR:
            return new HierarchicalPortPositionProcessor();

        case HIGH_DEGREE_NODE_LAYER_PROCESSOR:
            return new HighDegreeNodeLayeringProcessor();

        case HORIZONTAL_COMPACTOR:
            return new HorizontalGraphCompactor();

        case HYPEREDGE_DUMMY_MERGER:
            return new HyperedgeDummyMerger();

        case HYPERNODE_PROCESSOR:
            return new HypernodesProcessor();

        case IN_LAYER_CONSTRAINT_PROCESSOR:
            return new InLayerConstraintProcessor();

        case INNERMOST_NODE_MARGIN_CALCULATOR:
            return new InnermostNodeMarginCalculator();
            
        case INTERACTIVE_EXTERNAL_PORT_POSITIONER:
            return new InteractiveExternalPortPositioner();
            
        case INVERTED_PORT_PROCESSOR:
            return new InvertedPortProcessor();

        case LABEL_AND_NODE_SIZE_PROCESSOR:
            return new LabelAndNodeSizeProcessor();

        case LABEL_DUMMY_INSERTER:
            return new LabelDummyInserter();

        case LABEL_DUMMY_REMOVER:
            return new LabelDummyRemover();

        case LABEL_DUMMY_SWITCHER:
            return new LabelDummySwitcher();
            
        case LABEL_SIDE_SELECTOR:
            return new LabelSideSelector();
        
        case END_LABEL_SORTER:
            return new EndLabelSorter();

        case LAYER_CONSTRAINT_POSTPROCESSOR:
            return new LayerConstraintPostprocessor();

        case LAYER_CONSTRAINT_PREPROCESSOR:
            return new LayerConstraintPreprocessor();
            
        case LAYER_SIZE_AND_GRAPH_HEIGHT_CALCULATOR:
            return new LayerSizeAndGraphHeightCalculator();

        case LONG_EDGE_JOINER:
            return new LongEdgeJoiner();

        case LONG_EDGE_SPLITTER:
            return new LongEdgeSplitter();
            
        case NODE_PROMOTION:
            return new NodePromotion();

        case NORTH_SOUTH_PORT_POSTPROCESSOR:
            return new NorthSouthPortPostprocessor();

        case NORTH_SOUTH_PORT_PREPROCESSOR:
            return new NorthSouthPortPreprocessor();

        case ONE_SIDED_GREEDY_SWITCH:
            return new LayerSweepCrossingMinimizer(CrossMinType.ONE_SIDED_GREEDY_SWITCH);
            
        case PARTITION_MIDPROCESSOR:
            return new PartitionMidprocessor();

        case PARTITION_POSTPROCESSOR:
            return new PartitionPostprocessor();

        case PARTITION_PREPROCESSOR:
            return new PartitionPreprocessor();

        case PORT_LIST_SORTER:
            return new PortListSorter();

        case PORT_SIDE_PROCESSOR:
            return new PortSideProcessor();

        case REVERSED_EDGE_RESTORER:
            return new ReversedEdgeRestorer();
            
        case SELF_LOOP_PREPROCESSOR:
            return new SelfLoopPreProcessor();
            
        case SELF_LOOP_PORT_RESTORER:
            return new SelfLoopPortRestorer();
        
        case ALTERNATING_LAYER_UNZIPPER:
            return new AlternatingLayerUnzipper();

        case SELF_LOOP_POSTPROCESSOR:
            return new SelfLoopPostProcessor();
            
        case SELF_LOOP_ROUTER:
            return new SelfLoopRouter();
            
        case SEMI_INTERACTIVE_CROSSMIN_PROCESSOR:
            return new SemiInteractiveCrossMinProcessor();

        case SINGLE_EDGE_GRAPH_WRAPPER:
            return new SingleEdgeGraphWrapper();
            
        case SORT_BY_INPUT_ORDER_OF_MODEL:
            return new SortByInputModelProcessor();
            
        case TWO_SIDED_GREEDY_SWITCH:
            return new LayerSweepCrossingMinimizer(CrossMinType.TWO_SIDED_GREEDY_SWITCH);

        default:
            throw new IllegalArgumentException(
                    "No implementation is available for the layout processor " + toString());
        }
    }
}
