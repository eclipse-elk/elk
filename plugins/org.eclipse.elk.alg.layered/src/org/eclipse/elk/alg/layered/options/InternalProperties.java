/*******************************************************************************
 * Copyright (c) 2010, 2015 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.alg.layered.options;

import java.util.ArrayList;
import java.util.EnumSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import org.eclipse.elk.alg.common.nodespacing.cellsystem.LabelCell;
import org.eclipse.elk.alg.layered.compound.CrossHierarchyEdge;
import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LLabel;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LNode.NodeType;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.intermediate.FinalSplineBendpointsCalculator;
import org.eclipse.elk.alg.layered.intermediate.LayerConstraintPreprocessor;
import org.eclipse.elk.alg.layered.intermediate.loops.SelfLoopHolder;
import org.eclipse.elk.alg.layered.intermediate.wrapping.BreakingPointInserter;
import org.eclipse.elk.alg.layered.p5edges.splines.SplineEdgeRouter;
import org.eclipse.elk.alg.layered.p5edges.splines.SplineSegment;
import org.eclipse.elk.core.alg.ILayoutProcessor;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.math.KVectorChain;
import org.eclipse.elk.core.options.LabelSide;
import org.eclipse.elk.core.options.PortConstraints;
import org.eclipse.elk.core.options.PortSide;
import org.eclipse.elk.graph.properties.IProperty;
import org.eclipse.elk.graph.properties.Property;

import com.google.common.collect.Multimap;

/**
 * Container for property definitions for internal use of the algorithm. These properties should
 * not be accessed from outside. Some properties here are redefinitions of layout options in
 * {@link LayeredOptions} to change their defaults.
 *
 * @author msp
 * @author cds
 * @author uru
 */
public final class InternalProperties {

    /**
     * The original object from which a graph element was created.
     */
    public static final IProperty<Object> ORIGIN = new Property<Object>("origin");
    
    /**
     * The {@link LGraph} the coordinates of a cross-hierarchy edge are relative to.
     */
    public static final IProperty<LGraph> COORDINATE_SYSTEM_ORIGIN = new Property<>("coordinateOrigin");

    /**
     * The list of layout processors executed for an input graph.
     */
    public static final IProperty<List<ILayoutProcessor<LGraph>>> PROCESSORS = new Property<>("processors");

    /**
     * Whether the original node an LNode was created from was a compound node or not. This might
     * influence certain layout decisions, such as where to place inside port labels so that they
     * don't overlap edges.
     */
    public static final IProperty<Boolean> COMPOUND_NODE = new Property<Boolean>("compoundNode",
            false);

    /**
     * Whether the original port an LPort was created from was a compound port with connections to
     * or from descendants of its node. This might influence certain layout decisions, such as where
     * to place its inside port label.
     */
    public static final IProperty<Boolean> INSIDE_CONNECTIONS = new Property<Boolean>(
            "insideConnections", false);

    /**
     * The original bend points of an edge.
     */
    public static final IProperty<KVectorChain> ORIGINAL_BENDPOINTS = new Property<KVectorChain>(
            "originalBendpoints");

    /**
     * In interactive layout settings, this property can be set to indicate where a dummy node that
     * represents an edge in a given layer was probably placed in that layer. This information can
     * be calculated during the crossing minimization phase and later be used by an interactive node
     * placement algorithm.
     */
    public static final IProperty<Double> ORIGINAL_DUMMY_NODE_POSITION = new Property<Double>(
            "originalDummyNodePosition");

    /**
     * The edge a label originally belonged to. This property was introduced to remember which
     * cross-hierarchy edge a label originally belonged to.
     */
    public static final IProperty<LEdge> ORIGINAL_LABEL_EDGE = new Property<LEdge>(
            "originalLabelEdge");

    /**
     * Edge labels represented by an edge label dummy node.
     */
    public static final IProperty<List<LLabel>> REPRESENTED_LABELS = new Property<>("representedLabels");
    
    /**
     * If any edge incident to a node has end labels, those are stored in one label cell per port by the end label
     * processors.
     */
    public static final IProperty<Map<LPort, LabelCell>> END_LABELS = new Property<>("endLabels");
    
    /**
     * The edge an end label belongs to.
     */
    public static final IProperty<LEdge> END_LABEL_EDGE = new Property<>("endLabel.origin");

    /**
     * The side (of an edge) a label is placed on.
     */
    public static final IProperty<LabelSide> LABEL_SIDE = new Property<>("labelSide", LabelSide.UNKNOWN);
    
    /**
     * The maximum thickness of any edge connected to a port.
     */
    public static final IProperty<Double> MAX_EDGE_THICKNESS = new Property<>("maxEdgeThickness", 0.0);

    /**
     * Flag for reversed edges.
     */
    public static final IProperty<Boolean> REVERSED = new Property<Boolean>("reversed", false);

    /**
     * Random number generator for the algorithm.
     */
    public static final IProperty<Random> RANDOM = new Property<Random>("random");

    /**
     * The source port of a long edge before it was broken into multiple segments.
     */
    public static final IProperty<LPort> LONG_EDGE_SOURCE = new Property<>("longEdgeSource", null);
    
    /**
     * The target port of a long edge before it was broken into multiple segments.
     */
    public static final IProperty<LPort> LONG_EDGE_TARGET = new Property<>("longEdgeTarget", null);
    
    /**
     * Set on long edge dummies to indicate whether someplace in the long edge is a label dummy. This can for
     * example influence decisions on how to merge long edge dummies.
     */
    public static final IProperty<Boolean> LONG_EDGE_HAS_LABEL_DUMMIES = new Property<>("longEdgeHasLabelDummies",
            false);
    
    /**
     * If a long edge has a label dummy somewhere, this property indicates for every long edge dummy whether or
     * not it comes before that label dummy ({@code true}) or not ({@code false}). Coming before a label dummy
     * in this case means being in a layer left of the label dummy.
     */
    public static final IProperty<Boolean> LONG_EDGE_BEFORE_LABEL_DUMMY = new Property<>("longEdgeBeforeLabelDummy",
            false);

    /**
     * Edge constraints for nodes.
     */
    public static final IProperty<EdgeConstraint> EDGE_CONSTRAINT = new Property<EdgeConstraint>(
            "edgeConstraint", EdgeConstraint.NONE);

    /**
     * The layout unit a node belongs to. This property only makes sense for nodes. A layout unit is
     * a set of nodes between which no nodes belonging to other layout units may be placed. Nodes
     * not belonging to any layout unit may be placed arbitrarily between nodes of a layout unit.
     * Layer layout units are identified through one of their nodes.
     * The idea behind the latter is that for instance long edge dummies can be placed between 
     * an original node and one of its north port dummy nodes. It results in an edge crossing. 
     */
    public static final IProperty<LNode> IN_LAYER_LAYOUT_UNIT = new Property<LNode>(
            "inLayerLayoutUnit");

    /**
     * The in-layer constraint placed on a node. This indicates whether this node should be handled
     * like any other node, or if it must be placed at the top or bottom of a layer. This is
     * important for external port dummy nodes. Crossing minimizers are not required to respect this
     * constraint. If they don't, however, they must include a dependency on
     * {@link org.eclipse.elk.alg.layered.intermediate.InLayerConstraintProcessor}.
     */
    public static final IProperty<InLayerConstraint> IN_LAYER_CONSTRAINT
           = new Property<InLayerConstraint>("inLayerConstraint", InLayerConstraint.NONE);

    /**
     * Indicates that a node {@code x} may only appear inside a layer before the node {@code y} the
     * property is set to. That is, having {@code x} appear after {@code y} would violate this
     * constraint. This property only makes sense for nodes.
     */
    public static final IProperty<List<LNode>> IN_LAYER_SUCCESSOR_CONSTRAINTS =
            new Property<List<LNode>>("inLayerSuccessorConstraint", new ArrayList<LNode>());
    
    /**
     * Set to {@code true} if there are in-layer successor constraints between regular nodes. If that
     * is the case, constraint resolvers may have to know.
     */
    public static final IProperty<Boolean> IN_LAYER_SUCCESSOR_CONSTRAINTS_BETWEEN_NON_DUMMIES =
            new Property<Boolean>("inLayerSuccessorConstraintBetweenNonDummies", false);

    /**
     * A property set on ports indicating a dummy node created for that port. This is not set for
     * all ports that have dummy nodes created for them.
     */
    public static final IProperty<LNode> PORT_DUMMY = new Property<LNode>("portDummy");

    /**
     * Crossing hint used for in-layer cross counting with northern and southern port dummies. This
     * is effectively the number of different ports a northern or southern port dummy represents.
     */
    public static final IProperty<Integer> CROSSING_HINT = new Property<Integer>("crossingHint", 0);

    /**
     * Flags indicating the properties of a graph.
     */
    public static final IProperty<Set<GraphProperties>> GRAPH_PROPERTIES =
            new Property<Set<GraphProperties>>("graphProperties",
                    EnumSet.noneOf(GraphProperties.class));

    /**
     * The side of an external port a dummy node was created for.
     */
    public static final IProperty<PortSide> EXT_PORT_SIDE = new Property<PortSide>(
            "externalPortSide", PortSide.UNDEFINED);

    /**
     * Original size of the external port a dummy node was created for.
     */
    public static final IProperty<KVector> EXT_PORT_SIZE = new Property<KVector>(
            "externalPortSize", new KVector());
    
    /**
     * External port dummies that represent northern or southern external ports are replaced by new
     * dummy nodes during layout. The original ones are kept in this property set on the graph.
     */
    public static final IProperty<List<LNode>> EXT_PORT_REPLACED_DUMMIES = new Property<>(
            "externalPortReplacedDummies");
    
    /**
     * External port dummies that represent northern or southern external ports are replaced by new
     * dummy nodes during layout. In these cases, this property is set to the original dummy node.
     */
    public static final IProperty<LNode> EXT_PORT_REPLACED_DUMMY = new Property<LNode>(
            "externalPortReplacedDummy");

    /**
     * The port sides of external ports a connected component connects to. This property is set on
     * the layered graph that represents a connected component and defaults to no connections. If a
     * connected component connects to an external port on the EAST side and to another external
     * port on the NORTH side, this enumeration will list both sides.
     */
    public static final IProperty<Set<PortSide>> EXT_PORT_CONNECTIONS =
            new Property<Set<PortSide>>("externalPortConnections", EnumSet.noneOf(PortSide.class));

    /**
     * The original position or position-to-node-size ratio of a port. This property has two use
     * cases:
     * <ol>
     * <li>For external port dummies. In this use case, the property gives the original position of
     * the external port (if port constraints are set to {@code FIXED_POS}) or the original
     * position-to-node-size ratio of the external port ((if port constraints are set to
     * {@code FIXED_RATIO}).</li>
     * <li>For ports of regular nodes with port constraints set to {@code FIXED_RATIO}. Since
     * regular nodes may be resized, the original ratio must be remembered for the new port position
     * to be determined.</li>
     * </ol>
     * <p>
     * This is a one-dimensional value since the side of the port determines the other dimension.
     * (For eastern and western ports, the x coordinate is determined automatically; for northern
     * and southern ports, the y coordinate is determined automatically)
     * </p>
     */
    public static final IProperty<Double> PORT_RATIO_OR_POSITION = new Property<Double>(
            "portRatioOrPosition", 0.0);

    /**
     * A list of nodes whose barycenters should go into the barycenter calculation of the node this
     * property is set on. Nodes in this list are expected to be in the same layer as the node the
     * property is set on. This is primarily used when edges are rerouted from a node to dummy
     * nodes.
     * <p>
     * This property is currently not declared as one of the layout options offered by ELK Layered
     * and should be considered highly experimental.
     */
    public static final IProperty<List<LNode>> BARYCENTER_ASSOCIATES = new Property<List<LNode>>(
            "barycenterAssociates");

    /**
     * List of comment boxes that are placed on top of a node.
     */
    public static final IProperty<List<LNode>> TOP_COMMENTS = new Property<List<LNode>>(
            "TopSideComments");

    /**
     * List of comment boxes that are placed in the bottom of of a node.
     */
    public static final IProperty<List<LNode>> BOTTOM_COMMENTS = new Property<List<LNode>>(
            "BottomSideComments");

    /**
     * The port of a node that originally connected a comment box with that node.
     */
    public static final IProperty<LPort> COMMENT_CONN_PORT = new Property<LPort>(
            "CommentConnectionPort");

    /**
     * Whether a port is used to collect all incoming edges of a node.
     */
    public static final IProperty<Boolean> INPUT_COLLECT = new Property<Boolean>("inputCollect",
            false);
    /**
     * Whether a port is used to collect all outgoing edges of a node.
     */
    public static final IProperty<Boolean> OUTPUT_COLLECT = new Property<Boolean>("outputCollect",
            false);

    /**
     * Property of a LayeredGraph. Whether the graph has been processed by the cycle breaker and the
     * cycle breaker has detected cycles and reverted edges.
     */
    public static final IProperty<Boolean> CYCLIC = new Property<Boolean>("cyclic", false);

    /**
     * Map of original hierarchy crossing edges to a set of dummy edges by which the original edge
     * has been replaced.
     */
    public static final IProperty<Multimap<LEdge, CrossHierarchyEdge>> CROSS_HIERARCHY_MAP =
            new Property<Multimap<LEdge, CrossHierarchyEdge>>("crossHierarchyMap");

    /**
     * Offset to be added to the target anchor point of an edge when the layout is applied back to
     * the origin.
     */
    public static final IProperty<KVector> TARGET_OFFSET = new Property<KVector>("targetOffset");

    /**
     * Combined size of all edge labels of a spline self loop.
     */
    public static final IProperty<KVector> SPLINE_LABEL_SIZE =
            new Property<KVector>("splineLabelSize", new KVector());

    /**
     * Internal container for all possible spacing variations that we support.
     */
    public static final IProperty<Spacings> SPACINGS =
            new Property<Spacings>("spacings");

    /**
     * Specifies if the corresponding LGraph element was added by the
     * {@link org.eclipse.elk.alg.layered.intermediate.PartitionPreprocessor PartitionPreprocessor}.
     */
    public static final IProperty<Boolean> PARTITION_DUMMY = new Property<Boolean>(
            "partitionConstraint", false);

    /**
     * Holds information about a breaking point, e.g. the start dummy node, the end dummy node 
     * and the original edge that was split.
     * 
     * @see org.eclipse.elk.alg.layered.intermediate.wrapping.BreakingPointInserter BreakingPointInserter
     */
    public static final IProperty<BreakingPointInserter.BPInfo> BREAKING_POINT_INFO =
            new Property<BreakingPointInserter.BPInfo>("breakingPoint.info");

    /**
     * The edge to which to assign the final bendpoints when edges are routed as splines. By default this is the first
     * edge of the {@link InternalProperties#SPLINE_EDGE_CHAIN}, in which case this property should be {@code null}.
     * However, in certain scenarios, e.g. if graph wrapping is used, the edge must be specified explicitly. This must
     * be done using this property.
     */
    public static final IProperty<LEdge> SPLINE_SURVIVING_EDGE = new Property<>("splines.survivingEdge");
    
    /**
     * Collection of all spline routes created by the {@link SplineEdgeRouter}. Set on an {@link LEdge}, to be read by
     * the {@link FinalSplineBendpointsCalculator}.
     */
    public static final IProperty<List<SplineSegment>> SPLINE_ROUTE_START = new Property<>("splines.route.start");

    /**
     * Collection of all edge chains as computed by the {@link SplineEdgeRouter}. Set on an {@link LEdge}, to be read by
     * the {@link FinalSplineBendpointsCalculator}.
     */
    public static final IProperty<List<LEdge>> SPLINE_EDGE_CHAIN = new Property<>("splines.edgeChain");
    
    /**
     * Holds the information about the originally user set port constraints.
     */
    public static final IProperty<PortConstraints> ORIGINAL_PORT_CONSTRAINTS =  
            new Property<PortConstraints>("originalPortConstraints");
    
    /**
     * Holds all of the information necessary to route self loops around a node.
     */
    public static final IProperty<SelfLoopHolder> SELF_LOOP_HOLDER =
            new Property<>("selfLoopHolder");

    /**
     * Holds the y-coordinate of a deleted {@link NodeType#NORTH_SOUTH_PORT} dummy node. To be read by the
     * {@link FinalSplineBendpointsCalculator}.
     */
    public static final IProperty<Double> SPLINE_NS_PORT_Y_COORD = new Property<>("splines.nsPortY");
    
    /**
     * Set for nodes and edges to preserve the order in the model file.
     */
    public static final IProperty<Integer> MODEL_ORDER = new Property<>("modelOrder");

    /**
     * Set on ports to save their (long edge) target node.
     * Only used during port sorting if the node/edge order shall be preserved.
     */
    public static final IProperty<LNode> LONG_EDGE_TARGET_NODE = new Property<>("longEdgeTargetNode");
    
    /**
     * Set during crossing minimization on a graph to signal that crossing minimization should not shuffle the free
     * nodes.
     * This also fixes the initial sweep direction to forward.
     */
    public static final IProperty<Boolean> FIRST_TRY_WITH_INITIAL_ORDER =
            new Property<>("firstTryWithInitialOrder", false);
    
    /**
     * Set during crossing minimization on a graph to signal that crossing minimization should not shuffle the free
     * nodes.
     * This also fixes the initial sweep direction to backwards (since the previous one was forward).
     */
    public static final IProperty<Boolean> SECOND_TRY_WITH_INITIAL_ORDER =
            new Property<>("firstTryWithInitialOrder", false);

    /**
     * List of nodes hidden by the {@link LayerConstraintPreprocessor}. Set on graphs.
     */
    public static final IProperty<List<LNode>> HIDDEN_NODES = new Property<>("layerConstraints.hiddenNodes");
    
    /**
     * The port a hidden edge was originally connected to before {@link LayerConstraintPreprocessor} temporarily
     * disconnected it. Set on edges.
     */
    public static final IProperty<LPort> ORIGINAL_OPPOSITE_PORT = new Property<>("layerConstraints.opposidePort");
    
    /**
     * Holds the model order of the target nodes of a port.
     * This is used to cache this once it was calculated.
     */
    public static final IProperty<Map<LNode, Integer>> TARGET_NODE_MODEL_ORDER = new Property<>("targetNode.modelOrder");
    
    /**
     * The weight of a node as used by the MedianHeuristic class.
     */
    public static final IProperty<Double> WEIGHT = new Property<>("medianHeuristic.weight");
    
    /**
     * Hidden default constructor.
     */
    private InternalProperties() {
    }
}
