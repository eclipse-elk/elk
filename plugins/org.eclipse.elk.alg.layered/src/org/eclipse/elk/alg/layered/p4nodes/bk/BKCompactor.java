/*******************************************************************************
 * Copyright (c) 2015 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.layered.p4nodes.bk;

import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Queue;

import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.options.EdgeStraighteningStrategy;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.Spacings;
import org.eclipse.elk.alg.layered.p4nodes.bk.BKAlignedLayout.HDirection;
import org.eclipse.elk.alg.layered.p4nodes.bk.BKAlignedLayout.VDirection;
import org.eclipse.elk.alg.layered.p4nodes.bk.ThresholdStrategy.NullThresholdStrategy;
import org.eclipse.elk.alg.layered.p4nodes.bk.ThresholdStrategy.SimpleThresholdStrategy;
import org.eclipse.elk.alg.layered.options.LayeredOptions;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

/**
 * For documentation see {@link BKNodePlacer}.
 * 
 * As opposed to the default {@link BKCompactor} this version
 * trades maximal compactness with straight edges. In other words,
 * where possible it favors additional straight edges over compactness. 
 */
public class BKCompactor implements ICompactor {
    
    /** The graph to process. */
    private LGraph layeredGraph;
    /** Specific {@link ThresholdStrategy} to be used for execution. */
    private ThresholdStrategy threshStrategy;
    /** Information about a node's neighbors and index within its layer. */
    private NeighborhoodInformation ni;
    /** Spacings. */
    private Spacings spacings;
    /** Representation of the class graph. */
    private Map<LNode, ClassNode> sinkNodes = Maps.newHashMap();
    
    /**
     * @param layeredGraph the graph to handle.
     * @param ni precalculated information about a node's neighbors.
     */
    public BKCompactor(final LGraph layeredGraph, final NeighborhoodInformation ni) {
        this.layeredGraph = layeredGraph;
        this.ni = ni;
        spacings = layeredGraph.getProperty(InternalProperties.SPACINGS);
        
        // configure the requested threshold strategy
        if (layeredGraph.getProperty(LayeredOptions.NODE_PLACEMENT_BK_EDGE_STRAIGHTENING) 
                == EdgeStraighteningStrategy.IMPROVE_STRAIGHTNESS) {
            
            threshStrategy = new SimpleThresholdStrategy();
        } else {
            // mimics the original compaction strategy without additional straightening
            threshStrategy = new NullThresholdStrategy();
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////
    // Block Placement
    
    /**
     * In this step, actual coordinates are calculated for blocks and its nodes.
     * 
     * <p>First, all blocks are placed, trying to avoid any crossing of the blocks. Then, the blocks are
     * shifted towards each other if there is any space for compaction.</p>
     * 
     * @param bal One of the four layouts which shall be used in this step
     */
    public void horizontalCompaction(final BKAlignedLayout bal) {
        // Initialize fields with basic values, partially depending on the direction
        for (Layer layer : layeredGraph.getLayers()) {
            for (LNode node : layer.getNodes()) {
                bal.sink[node.id] = node;
                bal.shift[node.id] = bal.vdir == VDirection.UP
                        ? Double.NEGATIVE_INFINITY : Double.POSITIVE_INFINITY;
            }
        }
        // clear any previous sinks
        sinkNodes.clear();

        // If the horizontal direction is LEFT, the layers are traversed from right to left, thus
        // a reverse iterator is needed (note that this does not change the original list of layers)
        List<Layer> layers = layeredGraph.getLayers();
        if (bal.hdir == HDirection.LEFT) {
            layers = Lists.reverse(layers);
        }

        // init threshold strategy
        threshStrategy.init(bal, ni);
        // mark all blocks as unplaced
        Arrays.fill(bal.y, null);
        
        for (Layer layer : layers) {
            // As with layers, we need a reversed iterator for blocks for different directions
            List<LNode> nodes = layer.getNodes();
            if (bal.vdir == VDirection.UP) {
                nodes = Lists.reverse(nodes);
            }
            
            // Do an initial placement for all blocks
            for (LNode v : nodes) {
                if (bal.root[v.id].equals(v)) {
                    placeBlock(v, bal);
                }
            }
        }
        
        // Try to compact classes by shifting them towards each other if there is space between them.
        // Other than the original algorithm we use a "class graph" here in conjunction with a longest
        // path layering based on previously calculated separations between any pair of adjacent classes.
        // This allows to have different node sizes and disconnected graphs.
        placeClasses(bal);
        
        // apply final coordinates
        for (Layer layer : layers) {
            for (LNode v : layer.getNodes()) {
                bal.y[v.id] = bal.y[bal.root[v.id].id];
                
                // If this is the root node of the block, check if the whole block can be shifted to
                // further compact the drawing (the block's non-root nodes will be processed later by
                // this loop and will thus use the updated y position calculated here)
                if (v.equals(bal.root[v.id])) {
                    double sinkShift = bal.shift[bal.sink[v.id].id];
                    
                    if ((bal.vdir == VDirection.UP && sinkShift > Double.NEGATIVE_INFINITY)
                     || (bal.vdir == VDirection.DOWN  && sinkShift < Double.POSITIVE_INFINITY)) {
                        
                        bal.y[v.id] = bal.y[v.id] + sinkShift;
                    }
                }
            }
        }
        
        // all blocks were placed, shift latecomers
        threshStrategy.postProcess();
        
    }

    /**
     * Blocks are placed based on their root node. This is done by going through all layers the block
     * occupies and moving the whole block upwards / downwards if there are blocks that it overlaps
     * with.
     * 
     * @param root The root node of the block (usually called {@code v})
     * @param bal One of the four layouts which shall be used in this step
     */
    // SUPPRESS CHECKSTYLE NEXT 1 MethodLength
    private void placeBlock(final LNode root, final BKAlignedLayout bal) { 
        // Skip if the block was already placed
        if (bal.y[root.id] != null) {
            return;
        }
        
        // Initial placement
        // As opposed to the original algorithm we cannot rely on the fact that 
        //  0.0 as initial block position is always feasible. This is due to 
        //  the inside shift allowing for negative block positions in conjunction with
        //  a RIGHT (bottom-to-top) traversal direction. Computing the minimum with 
        //  an initial position of 0.0 thus leads to wrong results.
        // The wrong behavior is documented in KIPRA-1426
        boolean isInitialAssignment = true;
        bal.y[root.id] = 0.0;
        
        // Iterate through block and determine, where the block can be placed (until we arrive at the
        // block's root node again)
        LNode currentNode = root;
        double thresh =
                bal.vdir == VDirection.DOWN ? Double.NEGATIVE_INFINITY : Double.POSITIVE_INFINITY;
        do {
            int currentIndexInLayer = ni.nodeIndex[currentNode.id];
            int currentLayerSize = currentNode.getLayer().getNodes().size();

            // If the node is the top or bottom node of its layer, it can be placed safely since it is
            // the first to be placed in its layer. If it's not, we'll have to check its neighbours
            if ((bal.vdir == VDirection.DOWN && currentIndexInLayer > 0)
                    || (bal.vdir == VDirection.UP && currentIndexInLayer < (currentLayerSize - 1))) {

                // Get the node which is above / below the current node as well as the root of its block
                LNode neighbor = null;
                LNode neighborRoot = null;
                if (bal.vdir == VDirection.UP) {
                    neighbor = currentNode.getLayer().getNodes().get(currentIndexInLayer + 1);
                } else {
                    neighbor = currentNode.getLayer().getNodes().get(currentIndexInLayer - 1);
                }
                neighborRoot = bal.root[neighbor.id];
                
                // Ensure the neighbor was already placed
                placeBlock(neighborRoot, bal);
                
                // calculate threshold value for additional straight edges
                // this call has to be _after_ place block, otherwise the 
                // order of the elements in the postprocessing queue is wrong 
                thresh = threshStrategy.calculateThreshold(thresh, root, currentNode);
                
                // Note that the two nodes and their blocks form a unit called class in the original
                // algorithm. These are combinations of blocks which play a role in the final compaction
                if (bal.sink[root.id].equals(root)) {
                    bal.sink[root.id] = bal.sink[neighborRoot.id];
                }

                // Check if the blocks of the two nodes are members of the same class
                if (bal.sink[root.id].equals(bal.sink[neighborRoot.id])) {
                    // They are part of the same class
                    
                    // The minimal spacing between the two nodes depends on their node type
                    double spacing = spacings.getVerticalSpacing(currentNode, neighbor);
                    
                    // Determine the block's final position
                    if (bal.vdir == VDirection.UP) {
                        
                        double currentBlockPosition = bal.y[root.id];
                        double newPosition = bal.y[neighborRoot.id]
                                + bal.innerShift[neighbor.id]
                                - neighbor.getMargin().top
                                - spacing
                                - currentNode.getMargin().bottom
                                - currentNode.getSize().y
                                - bal.innerShift[currentNode.id];

                        if (isInitialAssignment) {
                            isInitialAssignment = false;
                            bal.y[root.id] = Math.min(newPosition, thresh);
                        } else {
                            bal.y[root.id] = Math.min(currentBlockPosition, 
                                                     Math.min(newPosition, thresh));
                        }
                        
                    } else { // DOWN
                        
                        double currentBlockPosition = bal.y[root.id];
                        double newPosition = bal.y[neighborRoot.id]
                                + bal.innerShift[neighbor.id]
                                + neighbor.getSize().y
                                + neighbor.getMargin().bottom
                                + spacing
                                + currentNode.getMargin().top
                                - bal.innerShift[currentNode.id];
                        
                        if (isInitialAssignment) {
                            isInitialAssignment = false;
                            bal.y[root.id] = Math.max(newPosition, thresh);
                        } else {
                            bal.y[root.id] = Math.max(currentBlockPosition, 
                                                     Math.max(newPosition, thresh));
                        }
                    }
                    
                } else { // CLASSES
                    
                    // They are not part of the same class. Compute how the two classes can be compacted
                    // later. Hence we determine a minimal required space between the two classes 
                    // relative two the two class sinks.
                    double spacing = layeredGraph.getProperty(LayeredOptions.SPACING_NODE_NODE).doubleValue();
                    
                    ClassNode sinkNode = getOrCreateClassNode(bal.sink[root.id], bal);
                    ClassNode neighborSink = getOrCreateClassNode(bal.sink[neighborRoot.id], bal);
                        
                    if (bal.vdir == VDirection.UP) {
                        
                        //  possible setup:
                        //  root         --> currentNode  
                        //  neighborRoot --> neighbor
                        double requiredSpace = 
                                bal.y[root.id]
                                + bal.innerShift[currentNode.id]
                                + currentNode.getSize().y
                                + currentNode.getMargin().bottom
                                + spacing
                                - (bal.y[neighborRoot.id]
                                   + bal.innerShift[neighbor.id]
                                   - neighbor.getMargin().top
                                   );

                        // add an edge to the class graph
                        sinkNode.addEdge(neighborSink, requiredSpace);
                        
                        // original algorithms procedure here:
                        // bal.shift[bal.sink[neighborRoot.id].id] =
                        // Math.max(bal.shift[bal.sink[neighborRoot.id].id], requiredSpace);

                    } else { // DOWN
                        //  possible setup:
                        //  neighborRoot --> neighbor 
                        //  root         --> currentNode
                        double requiredSpace =
                                bal.y[root.id]
                                + bal.innerShift[currentNode.id]
                                - currentNode.getMargin().top
                                - bal.y[neighborRoot.id]
                                - bal.innerShift[neighbor.id]
                                - neighbor.getSize().y
                                - neighbor.getMargin().bottom
                                - spacing;
                        
                        // add an edge to the class graph
                        sinkNode.addEdge(neighborSink, requiredSpace);
                        
                        // original algorithms procedure here:
                        // bal.shift[bal.sink[neighborRoot.id].id] =
                        // Math.min(bal.shift[bal.sink[neighborRoot.id].id], requiredSpace);
                    }
                }
            } else {
                thresh = threshStrategy.calculateThreshold(thresh, root, currentNode);
            }
            
            // Get the next node in the block
            currentNode = bal.align[currentNode.id];
        } while (currentNode != root);
        
        threshStrategy.finishBlock(root);
    }
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    // Class Placement
    
    private void placeClasses(final BKAlignedLayout bal) {
        
        // collect sinks of the class graph
        Queue<ClassNode> sinks = Lists.newLinkedList();
        for (ClassNode n : sinkNodes.values()) {
            if (n.indegree == 0) {
                sinks.add(n);
            }
        }
        
        // propagate shifts in a longest path layering fashion
        while (!sinks.isEmpty()) {
            ClassNode n = sinks.poll();
            
            // position the root of the class node tree
            if (n.classShift == null) {
                n.classShift = 0d;
            } 
            
            for (ClassEdge e : n.outgoing) {
                
                // initial position of a target does not depend on previous positions
                // (we need this as we cannot assume the top-most position to be 0)
                if (e.target.classShift == null) {
                    e.target.classShift = n.classShift + e.separation;
                } else if (bal.vdir == VDirection.DOWN) {
                    e.target.classShift = Math.min(e.target.classShift, n.classShift + e.separation);
                } else {
                    e.target.classShift = Math.max(e.target.classShift, n.classShift + e.separation);
                }
                 
                e.target.indegree--;
                
                if (e.target.indegree == 0) {
                    sinks.add(e.target);
                }
            }
        }
        
        // remember final shifts for all classes such that they 
        // can be applied as absolute coordinates
        for (ClassNode n : sinkNodes.values()) {
            bal.shift[n.node.id] = n.classShift;
        }
    }
    
    private ClassNode getOrCreateClassNode(final LNode sinkNode, final BKAlignedLayout bal) {
        ClassNode node = sinkNodes.get(sinkNode);
        if (node == null) {
            node = new ClassNode();
            node.node = sinkNode;
            sinkNodes.put(node.node, node);
        }
        return node;
    }
    
    /**
     * A node of the class graph.
     */
    private static class ClassNode {
        // SUPPRESS CHECKSTYLE NEXT 5 VisibilityModifier
        Double classShift = null;
        LNode node;
        List<ClassEdge> outgoing = Lists.newArrayList();
        int indegree = 0;
        
        private void addEdge(final ClassNode target, final double separation) {
            ClassEdge se = new ClassEdge();
            se.target = target;
            se.separation = separation;
            target.indegree++;
            outgoing.add(se);
        }
    }
    
    /** 
     * An edge of the class graph, holds the required separation
     * between the connected classes.
     */
    private static class ClassEdge {
        // SUPPRESS CHECKSTYLE NEXT 2 VisibilityModifier
        double separation = 0;
        ClassNode target;
    }
}
