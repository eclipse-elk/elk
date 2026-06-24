/*******************************************************************************
 * Copyright (c) 2016, 2017 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.alg.layered.p4nodes;

import static org.eclipse.elk.alg.layered.options.NodeFlexibility.getNodeFlexibility;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.ListIterator;
import java.util.Map;
import java.util.Queue;
import java.util.Set;
import java.util.Stack;
import java.util.function.Predicate;
import java.util.stream.StreamSupport;

import org.eclipse.elk.alg.common.networksimplex.NEdge;
import org.eclipse.elk.alg.common.networksimplex.NGraph;
import org.eclipse.elk.alg.common.networksimplex.NNode;
import org.eclipse.elk.alg.common.networksimplex.NetworkSimplex;
import org.eclipse.elk.alg.layered.LayeredPhases;
import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LGraphElement;
import org.eclipse.elk.alg.layered.graph.LLabel;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LNode.NodeType;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.intermediate.IntermediateProcessorStrategy;
import org.eclipse.elk.alg.layered.intermediate.LabelAndNodeSizeProcessor;
import org.eclipse.elk.alg.layered.intermediate.NorthSouthPortPreprocessor;
import org.eclipse.elk.alg.layered.options.GraphProperties;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.layered.options.NodeFlexibility;
import org.eclipse.elk.alg.layered.options.Spacings;
import org.eclipse.elk.core.alg.ILayoutPhase;
import org.eclipse.elk.core.alg.LayoutProcessorConfiguration;
import org.eclipse.elk.core.math.ElkMargin;
import org.eclipse.elk.core.options.NodeLabelPlacement;
import org.eclipse.elk.core.options.PortConstraints;
import org.eclipse.elk.core.options.PortSide;
import org.eclipse.elk.core.util.IElkProgressMonitor;

import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.google.common.math.DoubleMath;
/**
 * Implements the node placement strategy as described by Gansner et al. in the following paper. It
 * is based on the idea to convert the problem into an auxiliary graph which is layered using the
 * network simplex algorithm.
 * <ul>
 *   <li>Emden R. Gansner, Eleftherios Koutsofios, Stephen C. North, Kiem-Phong Vo, A technique for
 *       drawing directed graphs. <i>Software Engineering</i> 19(3), pp. 214-230, 1993.</li>
 * </ul>
 * 
 * <h2>Restrictions</h2> 
 * In this approach node positions are represented by layers, and layers are integral. 
 * This means in particular that we cannot support rational positions, neither for nodes nor for ports.
 * In case such values are present in the graph we resolve them as follows:
 * <ul>
 *   <li>Nodes: The position, height, margin, and spacing of a node are summed and ceiled.</li> 
 *   <li>Ports: A port's position is altered to the closest integer position.</li>
 * </ul>
 * 
 * <h2>Favor Straight Edges </h2>
 * The way the auxiliary graph is built allows for symmetric optimal solutions, which may result in edge "stair-cases".
 * Consider the following example. n1 and n4 cannot be moved closer to each other because of spacing restrictions. 
 * If the edge (n4,n3) is straight, several (vertical) positions of n2 result in minimal edge length.  
 * <pre>
 *   __
 *  |  |
 *  |n1|--.   __
 *  |__|  ^--|n2|--.   ____
 *                 ^--|    | 
 *   __               | n3 |
 *  |n4|--------------|____|
 * </pre>
 * 
 * With the {@link LayeredOptions#NODE_PLACEMENT_FAVOR_EDGE_STRAIGHTNESS} option set,
 * we <em>try</em> to improve this by additionally reducing the number of bends, like this:
 * <pre>
 *   __
 *  |  |      __
 *  |n1|-----|n2|--.
 *  |__|           |   ____
 *                 ^--|    | 
 *   __               | n3 |
 *  |n4|--------------|____|  
 * </pre>
 *
 * <h2>North/South Ports and Inverted Ports</h2> 
 * <p>Inverted ports are connected to a dummy node within the same layer by an in-layer edge. In-layer edges are mostly
 * ignored. However, to keep the edges short an edge is added to the network simplex graph. 
 * Since the relative order of original node and dummy node is already known, a simple edge (as opposed to the 
 * realization of an absolute value in the objective) is sufficient.</p> 
 * <p>North/south ports are realized using dummy nodes within the same layer as the port's parent node. The dummy 
 * is <b>not</b> connected to the parent (See {@link NorthSouthPortPreprocessor}). This is addressed here 
 * by adding extra edges to the network simplex graph with a rather small weight. When straight edges are preferred, 
 * vertical segments of north/south edges may be elongated in exchange for straighter horizontal edges.</p>
 * 
 * <h2>Flexible Ports</h2>
 * <p>The example above would allow all edges to be straight if n3's height is increased. Usually we consider 
 * the size of a node to be fixed and raise the port constraint to FIXED_POS before we start the node placement phase.
 * Setting the {@link LayeredOptions#NODE_PLACEMENT_NETWORK_SIMPLEX_NODE_FLEXIBILITY} layout option on a node
 * changes this assumption and allows to either reposition the ports on the node's border or even to increase 
 * the node's height. The option only takes effect if initial port constraints are at most FIXED_ORDER and,
 * in case of {@link NodeFlexibility#PORT_POSITION}, if the node is tall enough for the ports to be moved at all.
 * The behavior is deliberately <b>not</b> affected by {@link LayeredOptions#NODE_SIZE_CONSTRAINTS}.</p>
 * 
 * <h3>Implementation Details</h3>
 * <p>The flexibility of port positions and node sizes is realized by adding further nodes to the auxiliary graph
 * that represent the top and bottom border of a node. These nodes are connected by an edge of minimal length
 * guaranteeing a minimal size of the node. The network simplex method doesn't allow to specify a maximum length 
 * for an edge. If only the port positions are allowed to move but the node's size must not be changed, 
 * a very large weight is associated with an edge connecting the two aforementioned nodes. Thus it is likely,
 * but not guaranteed that the node's size doen't increase.</p>  
 * <p>Additionally, the various weights associated with the edges of the auxiliary graph must be carefully 
 * calibrated. For instance, in the flexible case, having a non-straight edge must be more costly than increasing 
 * the size of a node or moving a port.
 * </p> 
 *
 * <ul>
 * <li>"Regular edge": weight 4, 8, 32 (as opposed to 1, 2, 8 of the original paper).
 * <li>"Node size edge": flexible weight 1, static weight 10000 
 * <li>"Port spacing edge": weight 0
 * <li>"In-layer edge": weight 4
 * <li>"North-south edge": weight 0.1
 * <ul>
 * 
 * <h3>Current Limitations</h3>
 * <ul>
 *  <li>Port positions are not 'balanced' if they are allowed to move.</li> 
 *  <li>"Stair-cases" of long edges are possible, even with the 
 *      {@link LayeredOptions#NODE_PLACEMENT_FAVOR_STRAIGHT_EDGES} option set.</li>
 * <ul> 
 */
public class NetworkSimplexPlacer implements ILayoutPhase<LayeredPhases, LGraph> {
    
    /** Additional processor dependencies for graphs with hierarchical ports. */
    private static final LayoutProcessorConfiguration<LayeredPhases, LGraph> HIERARCHY_PROCESSING_ADDITIONS =
        LayoutProcessorConfiguration.<LayeredPhases, LGraph>create()
            .addBefore(LayeredPhases.P5_EDGE_ROUTING,
                    IntermediateProcessorStrategy.HIERARCHICAL_PORT_POSITION_PROCESSOR);

    /** The {@link LGraph} to be handled. */
    private LGraph lGraph;
    /** The spacings of the graph. */
    private Spacings spacings;
    /** The internally used network simplex graph. */
    private NGraph nGraph;
    
    /** Mapping of the internal representations of nodes, indexed by {@link LNode#id}. */
    private NodeRep[] nodeReps;
    /** Mapping of the internal representations of edges, indexed by {@link LEdge#id}. */
    private EdgeRep[] edgeReps;
    /** Mapping of network simplex nodes to the graph elements they represent. */
    private Map<LGraphElement, NNode> portMap = Maps.newHashMap();
    
    /** Node count of {@link #lGraph}. The field is uninitialized until {@link #buildInitialAuxiliaryGraph()} 
     *  has been executed. */
    private int nodeCount;
    private int edgeCount;

    // used for edge straightening
    /** Internal book keeping array for the case that {@link LayeredOptions#NODE_PLACEMENT_FAVOR_EDGE_STRAIGHTNESS} is
     *  used. Indexed by {@link LNode#id}. */
    private int[] nodeState;
    private List<Path> twoPaths;
    private boolean[] crossing;

    // used for a special version of node flexibility
    private Set<NEdge> flexibleWhereSpacePermitsEdges = Sets.newHashSet();
    
    // - - - - - - edge weights used in the auxiliary network simplex graph - - - - - -  
    /** Basis for the weight of edges. */
    private static final double EDGE_WEIGHT_BASE = 4;
    /** Smaller weight than default, since horizontal edges are more important. */
    private static final double SMALL_EDGE_WEIGHT = 0.1d;
    /** If this factor is smaller than one straight long edges are deemed more important
     *  than straight (node) paths (longer than two) identified during 
     *  {@link #preferStraightEdges(IElkProgressMonitor)}. */
    private static final double LONG_EDGE_VS_PATH_FACTOR = 2;
    
    /** Large weight to be applied if nodes must not change in size. */
    private static final double NODE_SIZE_WEIGHT_STATIC = 10000;
    /** Large weight to be applied if nodes must not change in size. */
    private static final double NODE_SIZE_WEIGHT_FLEXIBLE = 1;
    
    /** Epsilon for double equality testing. */
    private static final double EPSILON = 0.00001d;
    
    @Override
    public LayoutProcessorConfiguration<LayeredPhases, LGraph> getLayoutProcessorConfiguration(final LGraph graph) {
        if (graph.getProperty(InternalProperties.GRAPH_PROPERTIES)
                .contains(GraphProperties.EXTERNAL_PORTS)) {
            return HIERARCHY_PROCESSING_ADDITIONS;
        } else {
            return null;
        }
    }
    
    @Override
    public void process(final LGraph layeredGraph, final IElkProgressMonitor progressMonitor) {
        progressMonitor.begin("Network simplex node placement", 1);
        
        this.lGraph = layeredGraph;
        this.spacings = layeredGraph.getProperty(InternalProperties.SPACINGS);
        
        // -------------------------------
        // #1 build the auxiliary graph
        // -------------------------------
        prepare();
        
        buildInitialAuxiliaryGraph();
        
        insertNorthSouthAuxiliaryEdges();
        insertInLayerEdgeAuxiliaryEdges();
        
        if (lGraph.getProperty(LayeredOptions.NODE_PLACEMENT_FAVOR_STRAIGHT_EDGES)) {
            IElkProgressMonitor pm = progressMonitor.subTask(1);
            pm.begin("Straight Edges Pre-Processing", 1);
            preferStraightEdges();
            pm.done();
        }
        
        // make sure the ngraph is connected. Cases where this doesn't have to be the case include
        //  hierarchical nodes with unconnected ports that are (in the case of hierarchical layout) 
        //  converted into unconnected dummy nodes
        nGraph.makeConnected();
        
        // --------------------------------
        // #2 execute the network simplex
        // --------------------------------
        // compared to {@link NetworkSimplexLayerer} a significantly larger iteration limit 
        // is selected here because the node placement uses an auxiliary graph with
        // larger node and edge count
        int iterLimit = layeredGraph.getProperty(LayeredOptions.THOROUGHNESS) * nGraph.nodes.size();
        
        NetworkSimplex.forGraph(nGraph)
            .withIterationLimit(iterLimit)
            .withBalancing(false)
            .execute(progressMonitor.subTask(1));
        
        // every individual node can be 'flexible where space permits'.
        // thus we cannot check for the property here but must rely on the fact that the 
        // affected nodes have been collected before. In fact, the edges that represent 
        // the size of the node were collected
        if (!flexibleWhereSpacePermitsEdges.isEmpty()) {

            IElkProgressMonitor pm = progressMonitor.subTask(1);
            pm.begin("Flexible Where Space Processing", 1);
            
            insertFlexibleWhereSpaceAuxiliaryEdges();
            // now the nodes may resize -> alter the weights
            for (NEdge edge : flexibleWhereSpacePermitsEdges) {
                edge.weight = NODE_SIZE_WEIGHT_FLEXIBLE;
            }

            // run network simplex a second time
            NetworkSimplex.forGraph(nGraph)
                .withIterationLimit(iterLimit)
                .withBalancing(false)
                .execute(pm.subTask(1));
            
            pm.done();
        }
      
        // post process 'two paths' that have been identified during the previous call of #preferStraightEdges()
        if (layeredGraph.getProperty(LayeredOptions.NODE_PLACEMENT_FAVOR_STRAIGHT_EDGES)) {
            IElkProgressMonitor pm = progressMonitor.subTask(1);
            pm.begin("Straight Edges Post-Processing", 1);
            postProcessTwoPaths();
            pm.done();
        }
        
        // --------------------------------
        // #3 apply positions
        // --------------------------------
        applyPositions();
        
        cleanup();
        
        progressMonitor.done();
    }
    
    // ------------------------------------------------------------------------------------------------
    //                                         Preparation
    // ------------------------------------------------------------------------------------------------
    private void prepare() {
        
        this.nGraph = new NGraph();
        
        // "integerify" port anchor and port positions
        //   note that margin.top and margin.bottom are not required to be integral
        //   since they do not influence the offset calculation for the edges
        // ... while we're at it, we assign ids to the nodes and edges 
        int nodeIdx = 0;
        int edgeIdx = 0;
        for (Layer l : lGraph) {
            for (LNode lNode : l) {
                lNode.id = nodeIdx++;
                for (LEdge e : lNode.getOutgoingEdges()) {
                    e.id = edgeIdx++;
                }

                // if a node is flexible, an edge attaches to the port itself within 
                //  the auxiliary graph, thus the anchor must be integer
                // otherwise the port position can be altered such that it accounts for the anchor's position as well
                boolean anchorMustBeInteger = isFlexibleNode(lNode);
                for (LPort p : lNode.getPorts()) {
                    if (anchorMustBeInteger) {
                        // anchor
                        double y = p.getAnchor().y;
                        if (y != Math.floor(y)) {
                            double offset = y - Math.round(y);
                            p.getAnchor().y -= offset;
                        }
                    }
                    
                    // port + anchor
                    double y = p.getPosition().y + p.getAnchor().y;
                    if (y != Math.floor(y)) {
                        double offset = y - Math.round(y);
                        p.getPosition().y -= offset;
                    }
                }
            }
        }

        this.nodeCount = nodeIdx;
        this.edgeCount = edgeIdx;
        this.nodeReps = new NodeRep[nodeIdx];
        this.edgeReps = new EdgeRep[edgeIdx];
        this.flexibleWhereSpacePermitsEdges.clear();
    }
    
    private void cleanup() {
        this.lGraph = null;
        this.nGraph = null;
        
        this.nodeReps = null;
        this.edgeReps = null;
        this.portMap.clear();

        this.nodeState = null;
        this.crossing = null;
        this.twoPaths = null;
        this.flexibleWhereSpacePermitsEdges.clear();
    }
    
    // ------------------------------------------------------------------------------------------------
    //                                      Auxiliary Graph
    // ------------------------------------------------------------------------------------------------  

    /**
     * Builds the auxiliary graph on which the network simplex algorithm is run.
     */
    private void buildInitialAuxiliaryGraph() {
        for (Layer l : lGraph) {
            transformLayer(l);
        }
        transformEdges();
    }
    
    /**
     * Transforms all nodes of the passed {@link Layer} according to their level of {@link PortConstraints} either into
     * their corresponding {@link NodeRep}. Additionally adds separating edges to assert the node order within each
     * layer and proper spacing.
     */
    private void transformLayer(final Layer layer) {
        NodeRep lastRep = null;
        for (LNode lNode : layer) {

            NodeRep nodeRep;
            if (isFlexibleNode(lNode)) {
                nodeRep = transformFixedOrderNode(lNode);
            } else {
                nodeRep = transformFixedPosNode(lNode);
            }
            nodeReps[lNode.id] = nodeRep;
            
            // if there is a previous node in the layer, we have to create a separation edge
            if (lastRep != null) {
                double spacing =
                          lastRep.origin.getMargin().bottom
                        + spacings.getVerticalSpacing(lastRep.origin, lNode)
                        + lNode.getMargin().top;
                
                if (!lastRep.isFlexible) {
                    // for non-flexible nodes their height must be included 
                    // in the minimal length of the separation edge
                    spacing += lastRep.origin.getSize().y;
                }

                NEdge.of()
                    .delta((int) Math.ceil(spacing))
                    .weight(0)
                    .source(lastRep.tail)
                    .target(nodeRep.head)
                    .create();
            }
            
            // remember current elements for next iteration
            lastRep = nodeRep;
        }
    }
    
    /**
     * @return a {@link NodeRep} instance that basically contains two references to {@link LNode}.
     */
    private NodeRep transformFixedPosNode(final LNode lNode) {
        
        NNode singleNode = NNode.of()
            .origin(lNode)
            .type("non-flexible")
            .create(nGraph);
        
        // register the ports with the node
        lNode.getPorts().stream()
            .filter(p -> PortSide.SIDES_EAST_WEST.contains(p.getSide()))
            .forEach(p -> portMap.put(p, singleNode));
        
        return new NodeRep(lNode, false, singleNode, singleNode);
    }
    
    /**
     * @return the corners that were created for {@link LNode}.
     */
    private NodeRep transformFixedOrderNode(final LNode lNode) {

        // -----------------------------------
        //          corner creation
        // -----------------------------------
        NNode topLeft = NNode.of()
                .origin(lNode)
                .type("flexible-head")
                .create(nGraph); 
        NNode bottomLeft = NNode.of()
                .origin(lNode)
                .type("flexible-tail")
                .create(nGraph);
        NodeRep corners = new NodeRep(lNode, true, topLeft, bottomLeft);
        
        // -----------------------------------
        //       weight & minimum length
        // -----------------------------------
        // the LabelAndNodeSizeProcessor ensures that any node size constraint is satisfied,
        //  i.e. the node's size is at least the specified minimum size and, if desired,
        //  the node is resized according to owned ports and labels
        // as such we can safely use the current height as minimum height
        double minHeight = lNode.getSize().y;

        NodeFlexibility nf = getNodeFlexibility(lNode);
        double sizeWeight = NODE_SIZE_WEIGHT_STATIC;
        if (nf.isFlexibleSize()) {
            // we are allowed to enlarge to node
            //  nevertheless, a little weight is good, otherwise the node can become arbitrarily tall
            //  especially since the edges that connect the node's do not carry nay weight
            sizeWeight = NODE_SIZE_WEIGHT_FLEXIBLE; 
        }
        
        NEdge nodeSizeEdge = NEdge.of()
            .weight(sizeWeight)
            .delta((int) Math.ceil(minHeight))
            .source(topLeft)
            .target(bottomLeft)
            .create();

        if (nf == NodeFlexibility.NODE_SIZE_WHERE_SPACE_PERMITS) {
            flexibleWhereSpacePermitsEdges.add(nodeSizeEdge);
        }
        
        // -----------------------------------
        //          port transformation
        // -----------------------------------
        // convert the ports to NNodes, note that the list of westward ports 
        // must be reversed since their original order is from bottom to top
        transformPorts(Lists.reverse(lNode.getPortSideView(PortSide.WEST)), corners);
        transformPorts(lNode.getPortSideView(PortSide.EAST), corners);

        // return the corners for further processing
        return corners;
    }
    
    private void transformPorts(final Iterable<LPort> ports, final NodeRep corners) {
        
        if (Iterables.isEmpty(ports)) {
            // nothing to do ...
            // the top and bottom border of the node are already safely spaced apart
            return;
        }
         
        final double portSpacing = Spacings.getIndividualOrDefault(corners.origin, LayeredOptions.SPACING_PORT_PORT);
        ElkMargin portSurrounding =
                Spacings.getIndividualOrDefault(corners.origin, LayeredOptions.SPACING_PORTS_SURROUNDING);
        if (portSurrounding == null) {
            // No additional port spacing set
            portSurrounding = new ElkMargin();
        }
        
        NNode lastNNode = corners.head;
        LPort lastPort = null;
        for (LPort port : ports) {

            // spacing between the current pair of ports (or to the top border of the node)
            double spacing = 0;
            if (lastPort == null) {
                spacing = portSurrounding.top;
            } else {
                spacing = portSpacing;
                spacing += lastPort.getSize().y;
                // TODO only if PORT_LABELS is set?
                // + lastPort.getMargin().bottom
                // + port.getMargin().top;
            }
            
            // create NNode for the port
            NNode nNode = NNode.of()
                    .origin(port)
                    .type("port")
                    .create(nGraph);
            portMap.put(port, nNode);
            
            // connect with previous NNode
            NEdge.of()
                .weight(0)
                .delta((int) Math.ceil(spacing))
                .source(lastNNode)
                .target(nNode)
                .create();
            
            lastPort = port;
            lastNNode = nNode;
        }
        
        // and connect to the bottom border 
        NEdge.of()
            .weight(0)
            .delta((int) Math.ceil(portSurrounding.bottom + lastPort.getSize().y))
            .source(lastNNode)
            .target(corners.tail)
            .create();
    }
    
    private void transformEdges() {
        lGraph.getLayers().stream()
            .flatMap(layer -> layer.getNodes().stream())
            .flatMap(node -> StreamSupport.stream(node.getOutgoingEdges().spliterator(), false))
            .filter(e -> isHandledEdge(e))
            .forEach(e -> transformEdge(e));
    }
    
    private void transformEdge(final LEdge lEdge) {
        // a dummy node
        NNode dummy = NNode.of()
                .type("edge")
                .create(nGraph);

        // calculate port offsets
        NodeRep srcRep = nodeReps[lEdge.getSource().getNode().id];
        NodeRep tgtRep = nodeReps[lEdge.getTarget().getNode().id];
        LPort srcPort = lEdge.getSource();
        LPort tgtPort = lEdge.getTarget();

        double srcOffset = srcPort.getAnchor().y;
        double tgtOffset = tgtPort.getAnchor().y;
        // for non-flexible nodes, ports are relative to node positions
        if (!srcRep.isFlexible) {
            srcOffset += srcPort.getPosition().y;
        }
        if (!tgtRep.isFlexible) {
            tgtOffset += tgtPort.getPosition().y;
        } 
        assert DoubleMath.fuzzyEquals(srcOffset - tgtOffset, Math.round(srcOffset - tgtOffset),
                EPSILON) : "Port positions must be integral";
        
        int tgtDelta = (int) Math.max(0, srcOffset - tgtOffset);
        int srcDelta = (int) Math.max(0, tgtOffset - srcOffset);

        double weight = getEdgeWeight(lEdge);
        
        // an edge to the source
        NEdge left = NEdge.of(lEdge)
            .weight(weight)
            .delta(srcDelta)
            .source(dummy)
            .target(portMap.get(lEdge.getSource()))
            .create();

        // an edge to the target
        NEdge right = NEdge.of(lEdge)
            .weight(weight)
            .delta(tgtDelta)
            .source(dummy)
            .target(portMap.get(lEdge.getTarget()))
            .create();
        
        // remember 
        EdgeRep edgeRep = new EdgeRep(lEdge, dummy, left, right);
        edgeReps[lEdge.id] = edgeRep;
    }
    
    /** Insert {@link NEdge}s to keep edges connected to inverted ports short. */
    private void insertInLayerEdgeAuxiliaryEdges() {
        lGraph.getLayers().stream()
            .flatMap(l -> l.getNodes().stream())
            .filter(n -> n.getType() == NodeType.NORMAL)
            .flatMap(n -> StreamSupport.stream(n.getConnectedEdges().spliterator(), false))
            .filter(e -> e.isInLayerEdge()) 
            .forEach(inLayerEdge -> {
                
                boolean srcIsDummy = inLayerEdge.getSource().getNode().getType() != NodeType.NORMAL; 

                LPort thePort = srcIsDummy ? inLayerEdge.getTarget() : inLayerEdge.getSource();
                LNode dummyNode = inLayerEdge.getOther(thePort).getNode();
                
                NNode portRep = portMap.get(thePort);
                NNode dummyRep = nodeReps[dummyNode.id].head; // head/tail doesn't matter since it's a dummy node

                final NNode src, tgt;
                if (thePort.getNode().getIndex() < dummyNode.getIndex()) {
                    // port --> dummy
                    src = portRep;
                    tgt = dummyRep;
                } else {
                    // dummy --> port
                    src = dummyRep;
                    tgt = portRep;
                }
                
                NEdge.of()
                    .delta(0)
                    .weight(EDGE_WEIGHT_BASE)
                    .source(src)
                    .target(tgt) 
                    .create();
            });
    }
    
    /** Insert {@link NEdge}s to keep north and south port edges short. */
    private void insertNorthSouthAuxiliaryEdges() {
        lGraph.getLayers().stream()
            .flatMap(l -> l.getNodes().stream())
            .forEach(n -> {
                for (LPort sp : n.getPortSideView(PortSide.SOUTH)) {
                    LNode other = sp.getProperty(InternalProperties.PORT_DUMMY);
                    // if no edge was attached to the port, no dummy was created ...
                    if (other != null) {
                        NEdge.of()
                            .delta(0) // doesn't matter, separation is already taken care off
                            .weight(SMALL_EDGE_WEIGHT)
                            .source(nodeReps[n.id].tail)
                            .target(nodeReps[other.id].head)
                            .create();
                    }
                }
                
                for (LPort sp : n.getPortSideView(PortSide.NORTH)) {
                    LNode other = sp.getProperty(InternalProperties.PORT_DUMMY);
                    // if no edge was attached to the port, no dummy was created ...
                    if (other != null) {
                        NEdge.of()
                            .delta(0) // doesn't matter
                            .weight(SMALL_EDGE_WEIGHT)
                            .source(nodeReps[other.id].tail)
                            .target(nodeReps[n.id].head)
                            .create();
                    }
                }
            });
    }
    
    /** Inserts auxiliary edges for the case that {@link NodeFlexibility#NODE_SIZE_WHERE_SPACE_PERMITS} node exist. */
    private void insertFlexibleWhereSpaceAuxiliaryEdges() {

        int minLayer = nGraph.nodes.stream().map(n -> n.layer).min(Integer::compare).get();
        int maxLayer = nGraph.nodes.stream().map(n -> n.layer).max(Integer::compare).get();
        int usedLayers = maxLayer - minLayer;
        
        NNode globalSource = NNode.of().create(nGraph);
        NNode globalSink = NNode.of().create(nGraph);

        // make sure the distance between source and sink is preserved
        NEdge.of()
            .weight(NODE_SIZE_WEIGHT_STATIC * 2)
            .delta(usedLayers)
            .source(globalSource)
            .target(globalSink)
            .create();
        
        // fix the position of most non-flexible nodes and make sure the flexible nodes 
        // can only increase in size
        Arrays.stream(nodeReps)
            .filter(nr -> nr.origin.getType() == NodeType.NORMAL)
            .filter(nr -> nr.origin.getPorts().size() > 1) // allow leaves to move
            .forEach(nr -> {
                NEdge.of()
                    .weight(0)
                    .delta(nr.tail.layer - minLayer)
                    .source(globalSource)
                    .target(nr.tail)
                    .create();
            
                NEdge.of()
                    .weight(0)
                    .delta(usedLayers - nr.head.layer)
                    .source(nr.head)
                    .target(globalSink)
                    .create();
            });
    }
    
    // ------------------------------------------------------------------------------------------------
    //                                       Apply Layout
    // ------------------------------------------------------------------------------------------------
    
    private void applyPositions() {
        for (Layer l : lGraph) {
            for (LNode lNode : l) {
                // find the node's corners
                NodeRep nodeRep = nodeReps[lNode.id];
                double minY = nodeRep.head.layer;
                double maxY = nodeRep.tail.layer;
                
                // set new position and size
                lNode.getPosition().y = minY;

                double sizeDelta = (maxY - minY) - lNode.getSize().y; 
                
                boolean flexibleNode = isFlexibleNode(lNode);
                NodeFlexibility nf = getNodeFlexibility(lNode);
                
                // modify the size?
                if (flexibleNode && nf.isFlexibleSizeWhereSpacePermits()) {
                    lNode.getSize().y += sizeDelta;
                }

                // reposition ports if allowed
                if (flexibleNode && nf.isFlexiblePorts()) {
                    for (LPort p : lNode.getPorts()) {
                        if (PortSide.SIDES_EAST_WEST.contains(p.getSide())) {
                            NNode nNode = portMap.get(p);
                            p.getPosition().y = nNode.layer - minY;
                        }
                    }
                    // when the node got resized, the positions of labels and south ports have to be adjusted 
                    for (LLabel label : lNode.getLabels()) {
                        adjustLabelPosition(lNode, label, sizeDelta);
                    }
                    if (nf.isFlexibleSizeWhereSpacePermits()) {
                        lNode.getPortSideView(PortSide.SOUTH).forEach(p -> p.getPosition().y += sizeDelta);
                    }
                }
            }
        }
    }
    
    private void adjustLabelPosition(final LNode node, final LLabel label, final double sizeDelta) {
        Set<NodeLabelPlacement> placement = node.getProperty(LayeredOptions.NODE_LABELS_PLACEMENT);
        if (placement.contains(NodeLabelPlacement.V_BOTTOM)) {
            label.getPosition().y += sizeDelta;
        } else if (placement.contains(NodeLabelPlacement.V_CENTER)) {
            label.getPosition().y += sizeDelta / 2.0;
        }
        // V_TOP placement does not require adjustment
    }
    
    // ------------------------------------------------------------------------------------------------
    //                                        Convenience
    // ------------------------------------------------------------------------------------------------
    
    /**
     * @return for the passed edge, the individual edge weight as specified by Gansner et al. The
     *         idea is to use a higher weight for long edge dummies such that long edges are drawn
     *         straight with higher priority. Additionally, we consider the
     *         {@link org.eclipse.elk.alg.layered.options.PRIORITY PRIORITY} layout option.
     */
    private double getEdgeWeight(final LEdge edge) {
        int priority = Math.max(1, edge.getProperty(LayeredOptions.PRIORITY_STRAIGHTNESS));
        double edgeTypeWeight =
                getEdgeWeight(edge.getSource().getNode().getType(), edge.getTarget().getNode().getType());
        return priority * edgeTypeWeight;
    }
    
    private double getEdgeWeight(final NodeType nodeType1, final NodeType nodeType2) {
        if (nodeType1 == NodeType.NORMAL && nodeType2 == NodeType.NORMAL) {
            return 1 * EDGE_WEIGHT_BASE;
        } else if (nodeType1 == NodeType.NORMAL || nodeType2 == NodeType.NORMAL) {
            return 2 * EDGE_WEIGHT_BASE;
        } else {
            return 8 * EDGE_WEIGHT_BASE; // SUPPRESS CHECKSTYLE MagicNumber
        }
    }
    
    /**
     * @return true if the edge is neither a self loop, nor an in-layer edge.
     */
    private static boolean isHandledEdge(final LEdge edge) {
        return !edge.isSelfLoop() && !edge.isInLayerEdge();
    }
    
    // SUPPRESS CHECKSTYLE NEXT 10 Javadoc|VisibilityModifier
    private static class NodeRep {
        public LNode origin;
        /** The 'head' of the node, that is the border with the lower y coordinate. */
        public NNode head;
        /** The 'tail' of a node. That is, the border with the larger y coordinate. */
        public NNode tail;
        /** True if {@link #origin}'s {@link NodeFlexibility} doesn't equal {@link NodeFlexibility#NONE} .*/
        public boolean isFlexible;

        NodeRep(final LNode origin, final boolean isFlexible, final NNode top, final NNode bottom) {
            this.origin = origin;
            this.isFlexible = isFlexible;
            this.head = top;
            this.tail = bottom;
        }
    }
    
    // SUPPRESS CHECKSTYLE NEXT 5 Javadoc|VisibilityModifier
    private static class EdgeRep {
        @SuppressWarnings("unused")
        public LEdge origin;
        public NEdge left; 
        public NEdge right;
        
        EdgeRep(final LEdge origin, final NNode dummy, final NEdge left, final NEdge right) {
            this.origin = origin;
            this.left = left; 
            this.right = right;
        }
        
        /**
         * @return {@code true} if the edge's source and the edge's target have the same y coordinate. {@code false}
         *         otherwise.
         */
        public boolean isStraight() {
            return notStraightBy() == 0;
        }

        /**
         * @return an integer value indication the y coordinate difference between the edge's source and the edge's
         *         target.
         *         <ul>
         *         <li>A positive value indicates that the source has a larger y coordinate than the target.</li>
         *         <li>A negative value indicates that the source has a lower y coordinate than the target.</li>
         *         <li>Zero indicates that the edge is straight.</li>
         *         </ul>
         */
        public int notStraightBy() {
            return (left.target.layer - left.delta) - (right.target.layer - right.delta);
        }
    }
    
    // ------------------------------------------------------------------------------------------------
    //                                       Flexible Ports
    // ------------------------------------------------------------------------------------------------
    
    /**
     * A node is <em>flexible</em> if
     * <ul>
     *  <li>its port constraints are <b>not</b> {@link PortConstraints#FIXED_POS}, and</li>
     *  <li>the node's height is large enough to give the ports 
     *      on both sides (WEST and EAST) enough room to potentially alter their position.</li>
     * </ul>
     * 
     * The rationale for the latter case is the following:
     * if the height of a node is not enough to accommodate all ports 
     * and we are not allowed to change the node's size,
     * we simply use the port positions that have been computed by the {@link LabelAndNodeSizeProcessor}.
     * The processor is able to handle this situation properly.
     * 
     * @param lNode a node
     * @return whether the node is regarded flexible.
     */
    public static boolean isFlexibleNode(final LNode lNode) {

        // dummies are not flexible!
        if (lNode.getType() != NodeType.NORMAL) {
            return false;
        }
        
        // at least two ports are required ...
        if (lNode.getPorts().size() <= 1) {
            return false;
        }
        
        // if ports may not be moved there's no use in enlarging the node
        PortConstraints pc = lNode.getProperty(LayeredOptions.PORT_CONSTRAINTS);
        if (pc.isPosFixed()) {
            return false;
        }
        
        NodeFlexibility nf = getNodeFlexibility(lNode);
        if (nf == NodeFlexibility.NONE) {
            return false;
        }
        
        // if we cannot resize the node, and the given height is not enough to properly accommodate 
        // all ports, reuse the existing port positions
        if (!nf.isFlexibleSizeWhereSpacePermits()) {
            // we are not allowed to increase the node's size
            double portSpacing = Spacings.getIndividualOrDefault(lNode, LayeredOptions.SPACING_PORT_PORT);
            ElkMargin additionalPortSpacing = lNode.getProperty(LayeredOptions.SPACING_PORTS_SURROUNDING);
            if (additionalPortSpacing == null) {
                additionalPortSpacing = new ElkMargin(portSpacing, portSpacing, portSpacing, portSpacing);
            }

            // check west side
            List<LPort> westPorts = lNode.getPortSideView(PortSide.WEST);
            double requiredWestHeight = additionalPortSpacing.top 
                                      + additionalPortSpacing.bottom
                                      + (westPorts.size() - 1) * portSpacing;
            if (requiredWestHeight > lNode.getSize().y) {
                return false;
            }
            
            // check east side
            List<LPort> eastPorts = lNode.getPortSideView(PortSide.EAST);
            double requiredEastHeight = additionalPortSpacing.top 
                                      + additionalPortSpacing.bottom 
                                      + (eastPorts.size() - 1) * portSpacing;
            if (requiredEastHeight > lNode.getSize().y) {
                return false;
            }
        }
        return true;
    }
    
    // ------------------------------------------------------------------------------------------------
    //                                      Edge Straightening
    // ------------------------------------------------------------------------------------------------
    
    /** Convenience class, basically a type def for {@code ArrayList<LEdge>}. */
    @SuppressWarnings("serial")
    private static class Path extends ArrayList<LEdge> {

        /**
         * @return {@code true} if any node in this path is of type {@link NodeType#LONG_EDGE}. {@code false} otherwise.
         */
        public boolean containsLongEdgeDummy() {
            if (this.isEmpty()) {
                return false;
            }
            if (get(0).getSource().getNode().getType() == NodeType.LONG_EDGE) {
                return true;
            }
            return this.stream()
                    .map(e -> e.getTarget().getNode().getType())
                    .anyMatch(t -> t == NodeType.LONG_EDGE);
        }
        
        /**
         * @return {@code true} is any node of this path fulfills the passed predicate {@code p}.
         */
        public boolean containsFlexibleNode(final Predicate<NodeFlexibility> p) {
            if (this.isEmpty()) {
                return false;
            }
            NodeFlexibility nf = getNodeFlexibility(this.get(0).getSource().getNode());
            if (p.test(nf)) {
                return true;
            }
            return this.stream()
                    .map(e -> e.getTarget().getNode())
                    .anyMatch(n -> p.test(getNodeFlexibility(n)));
        }
        
        /**
         * Orders the two edges of this path such that the shared node is the target of the first edge 
         * and the source of the second edge.
         */
        public void orderTwoPath() {
            if (size() != 2) {
                throw new IllegalStateException("Order only allowed for two paths.");
            }
            LEdge first = get(0);
            LEdge second = get(1);
            if (first.getTarget().getNode() != second.getSource().getNode()) {
                clear();
                add(second);
                add(first);
            }
        }
        
        /**
         * @return {@code true} if the target node of the first edge {@link NetworkSimplexPlacer#isFlexibleNode(LNode)
         *         isFlexibleNode}. Call this method only if the path is ordered, see {@link #orderTwoPath()}. 
         */
        public boolean isTwoPathCenterNodeFlexible() {
            return isFlexibleNode(get(0).getTarget().getNode());
        }
    }
    
    /**
     * Potentially alters the weights of paths in a similar fashion to long edges. Also collects 'two paths' 
     * for later post processing.
     */
    private void preferStraightEdges() {

        // the nodes were counted and indexed during #prepare
        nodeState = new int[nodeCount];
        twoPaths = Lists.newArrayList();
        
        // record node states
        lGraph.getLayers().stream()
            .flatMap(l -> l.getNodes().stream())
            .forEach(n -> nodeState[n.id] = getNodeState(n));
        
        markEdgeCrossings();
        List<Path> identifiedPaths = identifyPaths();
        
        // essentially 'long paths' are treated like 'long edges': 
        // the weights of the connecting edges are altered to straighten the path as much as possible
        for (Path path : identifiedPaths) {

            if (path.size() <= 1) {
                continue;
            }
            
            // remember 'two paths' for processing after network simplex has been executed
            if (path.size() == 2) {
                path.orderTwoPath();
                if (!path.isTwoPathCenterNodeFlexible()) {
                    twoPaths.add(path);
                }
                continue;
            }

            // ignore paths that contain long edge dummies, and paths that contain flexible nodes that allow resizing
            if (path.containsLongEdgeDummy() 
                || path.containsFlexibleNode(nf -> nf.isFlexibleSizeWhereSpacePermits())) {
                continue;
            }
            
            Iterator<LEdge> pathIt = path.iterator();
            LEdge last = null;
            while (pathIt.hasNext()) {
                LEdge cur = pathIt.next();
                EdgeRep curRep = edgeReps[cur.id];

                double weight; 
                if (last == null || !pathIt.hasNext()) {
                    // first or last segment
                    weight = getEdgeWeight(NodeType.NORMAL, NodeType.LONG_EDGE);
                } else {
                    weight = getEdgeWeight(NodeType.LONG_EDGE, NodeType.LONG_EDGE);
                }
                
                // at this point one can decide whether long edges are more important than "paths"
                weight *= LONG_EDGE_VS_PATH_FACTOR;
                
                double oldLeftWeight = curRep.left.weight;
                curRep.left.weight = Math.max(oldLeftWeight, oldLeftWeight + (weight - oldLeftWeight));
                double oldRightWeight = curRep.right.weight;
                curRep.right.weight = Math.max(oldRightWeight, oldRightWeight + (weight - oldRightWeight));
                
                last = cur;
            }
        }
    }
    
    /**
     * Iteratively post processes 'two paths' until no further bends can be removed. 
     */
    private void postProcessTwoPaths() {
        Queue<Path> q = Lists.newLinkedList();
        q.addAll(twoPaths);
        
        Stack<Path> s = new Stack<>();
        while (!q.isEmpty()) {
            Path path = q.poll();
            boolean tryAgain = improveTwoPath(path, true);
            if (tryAgain) {
                s.add(path);
            }
        }
        
        while (!s.isEmpty()) {
            Path path = s.pop();
            improveTwoPath(path, false);
        }
    }
    
    /**
     * A two path <code>(u, center, v)</code> has the following form:
     * <pre>
     * u  <--leftEdge.left-- d1 --leftEdge.right--> center <--rightEdge.left-- d2 --rightEdge.right--> v
     * </pre>
     * 
     * The path is assumed to be ordered 'left to right', i.e. the center node is the target node of the first edge. 
     * See {@link Path#orderTwoPath()}. Additionally, the center node must not be flexible in any way.
     * 
     * @return {@code true} if nothing was changed and the path should be checked again later. 
     *         {@code false} if no further processing is required.
     */
    private boolean improveTwoPath(final Path path, final boolean probe) {
     
        EdgeRep leftEdge = edgeReps[path.get(0).id];
        EdgeRep rightEdge = edgeReps[path.get(1).id];

        // is the edge already straight?
        if (leftEdge.isStraight() && rightEdge.isStraight()) {
            return false;
        }
        
        // get center node
        Object centerOrigin = leftEdge.right.target.origin;
        if (!(centerOrigin instanceof LNode)) {
            return false;
        }
        // only two paths without flexible nodes are allowed here
        assert centerOrigin instanceof LNode;
        // can be a node or a port
        final LNode centerNode = (LNode) centerOrigin;
        final NodeRep nNode = nodeReps[centerNode.id];
        
        // identify on which side there is more space
        int nodeIndex = centerNode.getIndex();
        double aboveDist = Double.POSITIVE_INFINITY;
        if (nodeIndex > 0) {
            LNode above = centerNode.getLayer().getNodes().get(nodeIndex - 1);
            NodeRep aboveRep = nodeReps[above.id];
            double spacing = Math.ceil(spacings.getVerticalSpacing(above, centerNode));
            aboveDist = (nNode.head.layer - centerNode.getMargin().top) 
                        - (aboveRep.head.layer + above.getSize().y + above.getMargin().bottom)
                        - spacing;
        }
        double belowDist = Double.POSITIVE_INFINITY;
        if (nodeIndex < centerNode.getLayer().getNodes().size() - 1) {
            LNode below = centerNode.getLayer().getNodes().get(nodeIndex + 1);
            NodeRep belowRep = nodeReps[below.id];
            double spacing = Math.ceil(spacings.getVerticalSpacing(below, centerNode));
            belowDist = (belowRep.head.layer - below.getMargin().top) 
                        - (nNode.head.layer + centerNode.getSize().y + centerNode.getMargin().bottom)
                        - spacing;
        }
        
        // same space on both sides, check again later 
        if (probe && DoubleMath.fuzzyEquals(aboveDist, belowDist, EPSILON)) {
            return true; 
        }
        
        // from the way the auxiliary graph is built, one can differentiate four cases (A-D): 
        // Case A:          Case B:         Case C:         Case D: 
        //          o-v        u-o                             
        //         /              \            u-o   o-v         o-c-o
        //      o-c                c-o            \ /           /     \
        //     /                      \            c           u       v
        //    u                        v            
        // Note that case C and D only occur if node c is blocked either above or below, respectively. 
        // Otherwise the edge would be straight. 
        // 
        // In every case there are two non-straight edges. When routing edges orthogonally later on, 
        // one saves a bend point if only one of the four edges is non-straight. To straighten one of the 
        // two candidate edges, node c has to be moved up or down. It depends on the adjacent nodes in c's layer
        // if this is possible. 
        
        // the following variables represent the length of each of the four edges,
        // i.e. if possible node c has to be moved by this distance (either up or down)
        // O--a--o--b--O--c--o--d--O
        int a = +length(leftEdge.left);
        int b = -length(leftEdge.right);
        int c = -length(rightEdge.left);
        int d = +length(rightEdge.right);

        boolean caseD = (leftEdge.notStraightBy() > 0 && rightEdge.notStraightBy() < 0);
        boolean caseC = (leftEdge.notStraightBy() < 0 && rightEdge.notStraightBy() > 0);
        boolean caseB =
                leftEdge.left.target.layer + leftEdge.right.delta < rightEdge.right.target.layer + rightEdge.left.delta;
        boolean caseA =
                leftEdge.left.target.layer + leftEdge.right.delta > rightEdge.right.target.layer + rightEdge.left.delta;

        int move = 0;
        if (!caseD && !caseC) {
            if (caseA) {
                if (aboveDist + c > 0) {
                    move = c;
                } else if (belowDist - a > 0) {
                    move = a;
                }
            } else if (caseB) {
                if (aboveDist + b > 0) {
                    move = b;
                } else if (belowDist - d > 0) {
                    move = d;
                }
            }
        }

        // move the center node
        nNode.head.layer += move;
        if (nNode.isFlexible) {
            nNode.tail.layer += move;
        }
        
        return false;
    }

    /**
     * @return the length of the passed edge.
     */
    private static int length(final NEdge edge) {
        return Math.abs(edge.source.layer - edge.target.layer) - edge.delta;
    }
    
    // ------------------------------------------------------------------------------------------------
    //                                      Path identification
    // ------------------------------------------------------------------------------------------------
    
    private List<Path> identifyPaths() {
        final List<Path> paths = Lists.newArrayList();
        lGraph.getLayers().stream()
            .flatMap(l -> l.getNodes().stream()) 
            .filter(n -> nodeState[n.id] == JUNCTION)
            .forEach(junction -> {
                for (LEdge e : junction.getConnectedEdges()) {
                    if (!isHandledEdge(e)) {
                        continue;
                    }
                    Path path = follow(e, junction, new Path());
                    if (path.size() > 1) {
                        paths.add(path);
                    }
                }    
            });
        return paths;
    }
    
    private Path follow(final LEdge edge, final LNode current, final Path path) {
        LNode other = edge.getOther(current);
        path.add(edge); 

        // stop criteria
        if (nodeState[other.id] == VISITED 
                || nodeState[other.id] == JUNCTION 
                || crossing[edge.id]) {
            return path;    
        } 
        
        // recurse
        nodeState[other.id] = VISITED;
        for (LEdge incident : other.getConnectedEdges()) {
            if (!isHandledEdge(incident) || (incident == edge)) {
                continue;
            }
            return follow(incident, other, path);
        }
        
        return path;
    }
    
    /** Indicates that the node has not been visited yet. */
    private static final int VISITED = -1;
    /** Indicates that a node is not a {@link #JUNCTION}. */
    private static final int OTHER = 0;
    /** A junction has either an incoming degree larger than two, an outgoing degree larger than two, or both.
     *  Or it is a leaf node, i.e. has exactly one incident edge. */
    private static final int JUNCTION = 2;
    
    private static int getNodeState(final LNode node) {
        int inco = 0;
        int ouco = 0;
        for (LPort p : node.getPorts()) {
            inco += p.getIncomingEdges().stream().filter(e -> !e.isSelfLoop()).count();
            ouco += p.getOutgoingEdges().stream().filter(e -> !e.isSelfLoop()).count();
            if (inco > 1 || ouco > 1) {
                return JUNCTION;
            }
        }
        if (inco + ouco == 1) {
            return JUNCTION;
        }
        return OTHER;
    }
    
    // ------------------------------------------------------------------------------------------------
    //                                      Mark Crossings
    // ------------------------------------------------------------------------------------------------
    
    /**
     * Mark all edges that are involved in edge crossings, which is stored in the {@link #crossing} array.
     */
    private void markEdgeCrossings() {
        crossing = new boolean[edgeCount];
        lGraph.getLayers().stream().reduce((left, right) -> {
            markCrossingEdges(left, right); 
            return right;   
        });
    }
    
    /**
     * Marks the crossings of edges between the two passed layers.
     */
    private void markCrossingEdges(final Layer left, final Layer right) {
        
        final List<LEdge> openEdges = Lists.newArrayList();

        // add all edges in the order they occur in the left layer
        for (LNode node : left) {
            for (LPort port : node.getPortSideView(PortSide.EAST)) {
                for (LEdge edge : port.getOutgoingEdges()) {
                    if (edge.isInLayerEdge() || edge.isSelfLoop() 
                            || edge.getTarget().getNode().getLayer() != right) {
                        continue;
                    }
                    openEdges.add(edge);
                }
            }
        }
        
        // close the edges one after another, recording edge crossings
        for (LNode node : Lists.reverse(right.getNodes())) {
            for (LPort port : node.getPortSideView(PortSide.WEST)) { // don't reverse, bottom up is correct
                for (LEdge edge : port.getIncomingEdges()) {
                    if (edge.isInLayerEdge() || edge.isSelfLoop() 
                            || edge.getSource().getNode().getLayer() != left) {
                        continue;
                    }
                    if (!openEdges.isEmpty()) {
                        ListIterator<LEdge> openEdgesIt = openEdges.listIterator(openEdges.size()); // start at the end
                        LEdge last = openEdgesIt.previous();
                        while (last != edge && openEdgesIt.hasPrevious()) {
                            // mark both edges as being part of an edge crossing
                            crossing[last.id] = true;
                            crossing[edge.id] = true;
                            last = openEdgesIt.previous();
                        }
                        if (openEdgesIt.hasPrevious()) {
                            openEdgesIt.remove();
                        }
                    }
                }
            }
        }
    }
    
}
