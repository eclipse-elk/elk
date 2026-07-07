/*******************************************************************************
 * Copyright (c) 2012, 2015 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.layered.p4nodes.bk;

import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.eclipse.elk.alg.layered.LayeredPhases;
import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LNode.NodeType;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.intermediate.IntermediateProcessorStrategy;
import org.eclipse.elk.alg.layered.options.FixedAlignment;
import org.eclipse.elk.alg.layered.options.GraphProperties;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.layered.p4nodes.bk.BKAlignedLayout.HDirection;
import org.eclipse.elk.alg.layered.p4nodes.bk.BKAlignedLayout.VDirection;
import org.eclipse.elk.core.alg.ILayoutPhase;
import org.eclipse.elk.core.alg.LayoutProcessorConfiguration;
import org.eclipse.elk.core.util.IElkProgressMonitor;
import org.eclipse.elk.core.util.Pair;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;

/**
 * This algorithm is an implementation for solving the node placement problem
 * which is posed in phase 4 of the ELK Layered algorithm. Inspired by:
 * <ul>
 *   <li> Ulrik Brandes and Boris K&ouml;pf, Fast and simple horizontal coordinate assignment.
 *     In <i>Proceedings of the 9th International Symposium on Graph Drawing (GD'01)</i>,
 *     LNCS vol. 2265, pp. 33-36, Springer, 2002. </li>
 * </ul>
 * 
 * <p>The original algorithm was extended to be able to cope with ports, node sizes, and node margins,
 * and was made more stable in general. The algorithm is structured in five steps, which include two new
 * steps which were not included in the original algorithm by Brandes and Koepf. The middle three steps
 * are executed four times, traversing the graph in all combinations of TOP or BOTTOM and LEFT or
 * RIGHT.</p>
 * 
 * <p>In ELK Layered we have the general idea of layouting from left to right and
 * transforming in the desired direction later. We decided to translate the terminology of the original
 * algorithm which thinks of a layout from top to bottom. When placing coordinates, we have to differ
 * from the original algorithm, since node placement in ELK Layered has to assign y-coordinates and not
 * x-coordinates.</p>
 * 
 * <p>The variable naming in this code is mostly chosen for an iteration direction within our
 * left to right convention. Arrows describe iteration direction.
 * 
 * LEFT                  RIGHT
 * <----------           ------->
 * 
 * UP    ^              DOWN |
 *       |                   |
 *       |                   v
 * </p>   
 *      
 * <h4>The algorithm:</h4>
 * 
 * <p>The first step checks the graphs' edges and marks short edges which cross long edges (called
 * type 1 conflict). The algorithm indents to draw long edges as straight as possible, thus trying to
 * solve the marked conflicts in a way which keep the long edge straight.</p>
 * 
 * <p>============ UP, DOWN x LEFT, RIGHT ============</p>
 * 
 * <p>The second step traverses the graph in the given directions and tries to group connected nodes
 * into (horizontal) blocks. These blocks, respectively the contained nodes, will be drawn straight when
 * the algorithm is finished. Here, type 1 conflicts are resolved, so that the dummy nodes of a long
 * edge share the same block if possible, such that the long edge is drawn straightly.</p>
 * 
 * <p>The third step contains the addition of node size and port positions to the original algorithm.
 * Each block is investigated from top to bottom. Nodes are moved inside the blocks, such that the port
 * of the edge going to the next node is on the same level as that next node. Furthermore, the size of
 * the block is calculated, regarding node sizes and new space needed due to node movement.</p>
 * 
 * <p>In the fourth step, actual y-coordinates are determined. The blocks are placed, start block and
 * direction determined by the directions of the current iteration. It is tried to place the blocks as
 * compact as possible by grouping blocks.</p>
 *  
 * <p>======================= END =======================</p>
 * 
 * <p>The action of the last step depends on a layout option. If "fixedAlignment" is not set to 
 * BALANCED, one of the four calculated layouts is selected and applied, choosing the layout which 
 * uses the least space. If it is false, a balanced layout is chosen by calculating a median layout 
 * of all four layouts.</p>
 * 
 * <p>In rare cases, it is possible that one or more layouts is not correct, e.g. having nodes which
 * overlap each other or violating the layer ordering constraint. If the algorithm detects that, the
 * respective layout is discarded and another one is chosen.</p>
 * 
 * <dl>
 *   <dt>Preconditions:</dt>
 *     <dd>The graph has a proper layering with optimized nodes ordering</dd>
 *     <dd>Ports are properly arranged</dd>
 *   <dt>Postconditions:</dt>
 *     <dd>Each node is assigned a vertical coordinate such that no two nodes overlap</dd>
 *     <dd>The size of each layer is set according to the area occupied by its nodes</dd>
 *     <dd>The height of the graph is set to the maximal layer height</dd>
 * </dl>
 */
public final class BKNodePlacer implements ILayoutPhase<LayeredPhases, LGraph> {
    
    /** Additional processor dependencies for graphs with hierarchical ports. */
    private static final LayoutProcessorConfiguration<LayeredPhases, LGraph> HIERARCHY_PROCESSING_ADDITIONS =
            LayoutProcessorConfiguration.<LayeredPhases, LGraph>create()
                .addBefore(LayeredPhases.P5_EDGE_ROUTING,
                        IntermediateProcessorStrategy.HIERARCHICAL_PORT_POSITION_PROCESSOR);

    @Override
    public LayoutProcessorConfiguration<LayeredPhases, LGraph> getLayoutProcessorConfiguration(final LGraph graph) {
        if (graph.getProperty(InternalProperties.GRAPH_PROPERTIES).contains(
                GraphProperties.EXTERNAL_PORTS)) {
            return HIERARCHY_PROCESSING_ADDITIONS;
        } else {
            return null;
        }
    }
    
    private LGraph lGraph;
    /** List of edges involved in type 1 conflicts (see above). */
    private final Set<LEdge> markedEdges = Sets.newHashSet();
    /**  Precalculated information on nodes' neighborhoods etc. */
    private NeighborhoodInformation ni;

    /** Whether to produce a balanced layout or not. */
    private boolean produceBalancedLayout = false;
    
    @Override
    public void process(final LGraph layeredGraph, final IElkProgressMonitor monitor) {
        monitor.begin("Brandes & Koepf node placement", 1);

        this.lGraph = layeredGraph;
        
        // Precalculate some information that we require during the 
        // following processes. 
        ni = NeighborhoodInformation.buildFor(layeredGraph);

        // a balanced layout is desired if
        //  a) no specific alignment is set and straight edges are not desired
        //  b) a balanced alignment is enforced
        FixedAlignment align = layeredGraph.getProperty(LayeredOptions.NODE_PLACEMENT_BK_FIXED_ALIGNMENT);
        boolean favorStraightEdges = layeredGraph.getProperty(LayeredOptions.NODE_PLACEMENT_FAVOR_STRAIGHT_EDGES);
        produceBalancedLayout = 
                (align == FixedAlignment.NONE && !favorStraightEdges) 
                || align == FixedAlignment.BALANCED;

        // Phase which marks type 1 conflicts, no difference between the directions so only
        // one run is required.
        markConflicts(layeredGraph);

        // Initialize four layouts which result from the two possible directions respectively.
        BKAlignedLayout rightdown = null, rightup = null, leftdown = null, leftup = null;
        // SUPPRESS CHECKSTYLE NEXT MagicNumber
        List<BKAlignedLayout> layouts = Lists.newArrayListWithCapacity(4);
        switch (layeredGraph.getProperty(LayeredOptions.NODE_PLACEMENT_BK_FIXED_ALIGNMENT)) {
            case LEFTDOWN:
                leftdown =
                      new BKAlignedLayout(layeredGraph, ni.nodeCount, VDirection.DOWN, HDirection.LEFT);
                layouts.add(leftdown);
                break;
            case LEFTUP:
                leftup = 
                      new BKAlignedLayout(layeredGraph, ni.nodeCount, VDirection.UP, HDirection.LEFT);
                layouts.add(leftup);
                break;
            case RIGHTDOWN:
                rightdown = 
                      new BKAlignedLayout(layeredGraph, ni.nodeCount, VDirection.DOWN, HDirection.RIGHT);
                layouts.add(rightdown);
                break;
            case RIGHTUP:
                rightup = 
                      new BKAlignedLayout(layeredGraph, ni.nodeCount, VDirection.UP, HDirection.RIGHT);
                layouts.add(rightup); 
                break;
            default:
                leftdown = 
                   new BKAlignedLayout(layeredGraph, ni.nodeCount, VDirection.DOWN, HDirection.LEFT);
                leftup = 
                   new BKAlignedLayout(layeredGraph, ni.nodeCount, VDirection.UP, HDirection.LEFT);
                rightdown = 
                   new BKAlignedLayout(layeredGraph, ni.nodeCount, VDirection.DOWN, HDirection.RIGHT);
                rightup = 
                   new BKAlignedLayout(layeredGraph, ni.nodeCount, VDirection.UP, HDirection.RIGHT);
                layouts.add(rightdown);
                layouts.add(rightup);
                layouts.add(leftdown);
                layouts.add(leftup); 
        }
        
        BKAligner aligner = new BKAligner(layeredGraph, ni);
        for (BKAlignedLayout bal : layouts) {
            // Phase which determines the nodes' memberships in blocks. This happens in four different
            // ways, either from processing the nodes from the first layer to the last or vice versa.
            aligner.verticalAlignment(bal, markedEdges);
            
            // Additional phase which is not included in the original Brandes-Koepf Algorithm.
            // It makes sure that the connected ports within a block are aligned to avoid unnecessary
            // bend points. Also, the required size of each block is determined.
            aligner.insideBlockShift(bal);
        }

        ICompactor compacter = new BKCompactor(layeredGraph, ni);
        for (BKAlignedLayout bal : layouts) {
            // This phase determines the y coordinates of the blocks and thus the vertical coordinates
            // of all nodes.
            compacter.horizontalCompaction(bal);
        }

        // Debug output
        if (monitor.isLoggingEnabled()) {
            for (BKAlignedLayout bal : layouts) {
                monitor.log(bal + " size is " + bal.layoutSize());
            }
        }

        // Choose a layout from the four calculated layouts. Layouts that contain errors are skipped.
        // The layout with the smallest size is selected. If more than one smallest layout exists,
        // the first one of the competing layouts is selected.
        BKAlignedLayout chosenLayout = null;

        // If layout options chose to use the balanced layout, it is calculated and added here.
        // If it is broken for any reason, one of the four other layouts is selected by the
        // given criteria.
        if (produceBalancedLayout) {
            BKAlignedLayout balanced = createBalancedLayout(layouts, ni.nodeCount);
            if (checkOrderConstraint(layeredGraph, balanced, monitor)) {
                chosenLayout = balanced;
            }
        } 
        
        // Either if no balanced layout is requested, or, if the balanced layout
        // violates order constraints, pick the one with the smallest height
        if (chosenLayout == null) {
            for (BKAlignedLayout bal : layouts) {
                if (checkOrderConstraint(layeredGraph, bal, monitor)) {
                    if (chosenLayout == null || chosenLayout.layoutSize() > bal.layoutSize()) {
                        chosenLayout = bal;
                    }
                }
            }
        }
             
        // If no layout is correct (which should never happen but is not strictly impossible),
        // the RIGHTDOWN layout is chosen by default.
        if (chosenLayout == null) {
            chosenLayout = layouts.get(0); // there has to be at least one layout in the list
        }
        
        // Apply calculated positions to nodes.
        for (Layer layer : layeredGraph.getLayers()) {
            for (LNode node : layer.getNodes()) {
                node.getPosition().y = chosenLayout.y[node.id] + chosenLayout.innerShift[node.id];
            }
        }

        // Debug output
        if (monitor.isLoggingEnabled()) {
            monitor.log("Chosen node placement: " + chosenLayout);
            monitor.log("Blocks: " + getBlocks(chosenLayout));
            monitor.log("Classes: " + getClasses(chosenLayout, monitor));
            monitor.log("Marked edges: " + markedEdges);
        }
        
        // cleanup
        for (BKAlignedLayout bal : layouts) {
            bal.cleanup();
        }
        ni.cleanup();
        markedEdges.clear();
        
        monitor.done();
    }
    

    /////////////////////////////////////////////////////////////////////////////////////////////////////
    // Conflict Detection
    
    /** The minimum number of layers we need to have conflicts. */
    private static final int MIN_LAYERS_FOR_CONFLICTS = 3;

    /**
     * This phase of the node placer marks all type 1 and type 2 conflicts.
     * 
     * <p>The conflict types base on the distinction of inner segments and non-inner segments of edges.
     * A inner segment is present if an edge is drawn between two dummy nodes and thus is part of
     * a long edge. A non-inner segment is present if one of the connected nodes is not a dummy
     * node.</p>
     * 
     * <p>Type 0 conflicts occur if two non-inner segments cross each other. Type 1 conflicts happen 
     * when a non-inner segment and a inner segment cross. Type 2 conflicts are present if two
     * inner segments cross.</p>
     * 
     * <p>The markers are later used to solve conflicts in favor of long edges. In case of type 2
     * conflicts, the marker favors the earlier node in layout order.</p>
     * 
     * @param layeredGraph The layered graph to be layouted
     */
    private void markConflicts(final LGraph layeredGraph) {
        int numberOfLayers = layeredGraph.getLayers().size();
        
        // Check if there are enough layers to detect conflicts
        if (numberOfLayers < MIN_LAYERS_FOR_CONFLICTS) {
            return;
        }
        
        // We'll need the number of nodes in the different layers quite often in this method, so save
        // them up front
        int[] layerSize = new int[numberOfLayers];
        int layerIndex = 0;
        for (Layer layer : layeredGraph.getLayers()) {
            layerSize[layerIndex++] = layer.getNodes().size();
        }
        
        // The following call succeeds since there are at least 3 layers in the graph
        Iterator<Layer> layerIterator = layeredGraph.getLayers().listIterator(2);
        for (int i = 1; i < numberOfLayers - 1; i++) {
            // The variable naming here follows the notation of the corresponding paper
            // Normally, underscores are not allowed in local variable names, but since there
            // is no way of properly writing indices beside underscores, Checkstyle will be
            // disabled here and in future methods containing indexed variables
            // CHECKSTYLEOFF Local Variable Names
            Layer currentLayer = layerIterator.next();
            Iterator<LNode> nodeIterator = currentLayer.getNodes().iterator();
            
            int k_0 = 0;
            int l = 0;
            
            for (int l_1 = 0; l_1 < layerSize[i + 1]; l_1++) {
                // In the paper, l and i are indices for the layer and the position in the layer
                LNode v_l_i = nodeIterator.next(); // currentLayer.getNodes().get(l_1);
                
                if (l_1 == ((layerSize[i + 1]) - 1) || incidentToInnerSegment(v_l_i, i + 1, i)) {
                    int k_1 = layerSize[i] - 1;
                    if (incidentToInnerSegment(v_l_i, i + 1, i)) {
                        k_1 = ni.nodeIndex[ni.leftNeighbors.get(v_l_i.id).get(0).getFirst().id];
                    }
                    
                    while (l <= l_1) {
                        LNode v_l = currentLayer.getNodes().get(l);
                        
                        if (!incidentToInnerSegment(v_l, i + 1, i)) {
                            for (Pair<LNode, LEdge> upperNeighbor : ni.leftNeighbors.get(v_l.id)) {
                                int k = ni.nodeIndex[upperNeighbor.getFirst().id];
                                
                                if (k < k_0 || k > k_1) {
                                    // Marked edge can't return null here, because the upper neighbor
                                    // relationship between v_l and upperNeighbor enforces the existence
                                    // of at least one edge between the two nodes
                                    markedEdges.add(upperNeighbor.getSecond());
                                }
                            }
                        }
                        
                        l++;
                    }
                    
                    k_0 = k_1;
                }
            }
            // CHECKSTYLEON Local Variable Names
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////
    // Layout Balancing

    /**
     * Calculates a balanced layout by determining the median of the four layouts.
     * 
     * <p>First, the layout with the smallest height, meaning the difference between the highest and the
     * lowest y-coordinate placement, is used as a starting point. Then, the median position of each of
     * the four layouts is used to determine the final position.</p>
     * 
     * <p>During this process, a node's inner shift value is regarded.</p>
     * 
     * @param layouts The four calculated layouts
     * @param nodeCount The number of nodes in the graph
     * @return A balanced layout, the median of the four layouts
     */
    private BKAlignedLayout createBalancedLayout(final List<BKAlignedLayout> layouts,
            final int nodeCount) {
        
        final int noOfLayouts = layouts.size();
        BKAlignedLayout balanced = new BKAlignedLayout(lGraph, nodeCount, null, null);
        double[] width = new double[noOfLayouts];
        double[] min = new double[noOfLayouts];
        double[] max = new double[noOfLayouts];
        int minWidthLayout = 0;

        // Find the smallest layout
        for (int i = 0; i < noOfLayouts; i++) {
            min[i] = Integer.MAX_VALUE;
            max[i] = Integer.MIN_VALUE;
        }
        
        for (int i = 0; i < noOfLayouts; i++) {
            BKAlignedLayout bal = layouts.get(i);
            width[i] = bal.layoutSize();
            if (width[minWidthLayout] > width[i]) {
                minWidthLayout = i;
            }
            
            for (Layer l : lGraph) {
                for (LNode n : l) {
                    double nodePosY = bal.y[n.id] + bal.innerShift[n.id];
                    min[i] = Math.min(min[i], nodePosY);
                    max[i] = Math.max(max[i], nodePosY + n.getSize().y);
                }
            }
        }

        // Find the shift between the smallest and the four layouts
        double[] shift = new double[noOfLayouts];
        for (int i = 0; i < noOfLayouts; i++) {
            if (layouts.get(i).vdir == VDirection.DOWN) {
                shift[i] = min[minWidthLayout] - min[i];
            } else {
                shift[i] = max[minWidthLayout] - max[i];
            }
        }

        // Calculated y-coordinates for a balanced placement
        double[] calculatedYs = new double[noOfLayouts];
        for (Layer layer : lGraph.getLayers()) {
            for (LNode node : layer.getNodes()) {
                for (int i = 0; i < noOfLayouts; i++) {
                    // it's important to include the innerShift here!
                    calculatedYs[i] =
                            layouts.get(i).y[node.id] + layouts.get(i).innerShift[node.id] + shift[i];
                }
               
                Arrays.sort(calculatedYs);
                balanced.y[node.id] = (calculatedYs[1] + calculatedYs[2]) / 2.0;
                // since we include the inner shift in the calculation of a balanced y 
                // coordinate we don't need it any more
                // note that after this step no further processing of the graph that 
                // would include the inner shift is possible
                balanced.innerShift[node.id] = 0d;
            }
        }

        return balanced;
    }
    

    /////////////////////////////////////////////////////////////////////////////////////////////////////
    // Utility Methods

    /**
     * Checks whether the given node is part of a long edge between the two given layers.
     * At this 'layer2' is left, or before, 'layer1'.
     * 
     * @param node Possible long edge node
     * @param layer1 The first layer, the layer of the node
     * @param layer2 The second layer
     * @return True if the node is part of a long edge between the layers, false else
     */
    private boolean incidentToInnerSegment(final LNode node, final int layer1, final int layer2) {
        if (node.getType() == NodeType.LONG_EDGE) {
            for (LEdge edge : node.getIncomingEdges()) {
                NodeType sourceNodeType = edge.getSource().getNode().getType();
                
                if (sourceNodeType == NodeType.LONG_EDGE
                        && ni.layerIndex[edge.getSource().getNode().getLayer().id] == layer2
                        && ni.layerIndex[node.getLayer().id] == layer1) {
                    
                    return true;
                }
            }
        }
        return false;
    }
  
    /**
     * Find an edge between two given nodes.
     * 
     * @param source The source node of the edge
     * @param target The target node of the edge
     * @return The edge between source and target, or null if there is none
     */
    static LEdge getEdge(final LNode source, final LNode target) {
        for (LEdge edge : source.getConnectedEdges()) {
            if (edge.getTarget().getNode().equals(target) || edge.getSource().getNode().equals(target)) {
                return edge;
            }
        }
        
        return null;
    }

    /**
     * Finds all blocks of a given layout.
     * 
     * @param bal The layout of which the blocks shall be found
     * @return The blocks of the given layout
     */
    static Map<LNode, List<LNode>> getBlocks(final BKAlignedLayout bal) {
        Map<LNode, List<LNode>> blocks = Maps.newLinkedHashMap();
        
        for (Layer layer : bal.layeredGraph.getLayers()) {
            for (LNode node : layer.getNodes()) {
                LNode root = bal.root[node.id];
                List<LNode> blockContents = blocks.get(root);
                
                if (blockContents == null) {
                    blockContents = Lists.newArrayList();
                    blocks.put(root, blockContents);
                }
                
                blockContents.add(node);
            }
        }
        
        return blocks;
    }

    /**
     * Finds all classes of a given layout. Only used for debug output.
     * 
     * @param bal The layout whose classes to find
     * @param monitor progress monitor for debug output
     * @return The classes of the given layout
     */
    static Map<LNode, List<LNode>> getClasses(final BKAlignedLayout bal, final IElkProgressMonitor monitor) {
        Map<LNode, List<LNode>> classes = Maps.newLinkedHashMap();
        
        // We need to enumerate all block roots
        Set<LNode> roots = Sets.newLinkedHashSet(Arrays.asList(bal.root));
        for (LNode root : roots) {
            if (root == null) {
                monitor.log("There are no classes in a balanced layout.");
                break;
            }
            LNode sink = bal.sink[root.id];
            List<LNode> classContents = classes.get(sink);
            
            if (classContents == null) {
                classContents = Lists.newArrayList();
                classes.put(sink, classContents);
            }
            
            classContents.add(root);
        }
        
        return classes;
    }
        
    /**
     * Checks whether all nodes are placed in the correct order in their layers and do not overlap
     * each other.
     * 
     * @param layeredGraph the containing layered graph.
     * @param bal the layout which shall be checked.
     * @param monitor progress monitor for logging.
     * @return {@code true} if the order is preserved and no nodes overlap, {@code false} otherwise.
     */
    private boolean checkOrderConstraint(final LGraph layeredGraph, final BKAlignedLayout bal,
            final IElkProgressMonitor monitor) {
        
        // Flag indicating whether the layout is feasible or not
        boolean feasible = true;
        
        // Iterate over the layers
        for (Layer layer : layeredGraph.getLayers()) {
            // Current Y position in the layer
            double pos = Double.NEGATIVE_INFINITY;
            
            // We remember the previous node for debug output
            LNode previous = null;
            
            // Iterate through the layer's nodes
            for (LNode node : layer.getNodes()) {
                // For the layout to be correct, both the node's top border and its bottom border must
                // be beyond the current position in the layer
                double top = bal.y[node.id] + bal.innerShift[node.id] - node.getMargin().top;
                double bottom = bal.y[node.id] + bal.innerShift[node.id] + node.getSize().y
                        + node.getMargin().bottom;
                
                if (top > pos && bottom > pos) {
                    previous = node;
                    
                    // Update the position inside the layer
                    pos = bal.y[node.id] + bal.innerShift[node.id] + node.getSize().y
                            + node.getMargin().bottom;
                } else {
                    // We've found an overlap
                    feasible = false;
                    if (monitor.isLoggingEnabled()) {
                        monitor.log("bk node placement breaks on " + node
                                + " which should have been after " + previous);
                    }
                    break;
                }
            }
            
            // Don't bother continuing if we've already determined that the layout is infeasible
            if (!feasible) {
                break;
            }
        }
        
        if (monitor.isLoggingEnabled()) {
            monitor.log(bal + " is feasible: " + feasible);
        }
        
        return feasible;
    }
}
