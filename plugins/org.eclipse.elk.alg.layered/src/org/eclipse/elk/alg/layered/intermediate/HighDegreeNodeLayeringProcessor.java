/*******************************************************************************
 * Copyright (c) 2016 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.alg.layered.intermediate;

import java.util.List;
import java.util.ListIterator;

import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.core.alg.ILayoutProcessor;
import org.eclipse.elk.core.util.IElkProgressMonitor;
import org.eclipse.elk.core.util.Pair;

import java.util.function.Function;
import java.util.ArrayList;

/**
 * Tries to improve the legibility of graphs that contain nodes of a very high degree. With such
 * node it can happen that high degree nodes end up in layers with many leaf nodes of other high
 * degree nodes resulting in a very tall drawing. This processor tackles this by identifying high
 * degree nodes and trees connected to the high degree node. Trees that to not exceed a certain
 * height are moved to newly introduced layers next to the layer with a high degree node. If a layer
 * contains multiple high degree nodes the newly introduced layers are shared amongst the high
 * degree nodes which remain in their original layer.
 * 
 * Note that leafs of paths of nodes are special kinds of trees and thus treated as well.
 * 
 * <dl>
 *   <dt>Precondition:</dt>
 *     <dd>a properly layered graph.</dd>
 *   <dt>Postcondition:</dt>
 *     <dd>a properly layered graph with additional layers for trees connected to high degree nodes.</dd>
 *   <dt>Slots:</dt>
 *     <dd>After phase 2.</dd>
 *   <dt>Same-slot dependencies:</dt>
 *     <dd>Should most likey run before any other intermediate processor in this slot 
 *         since it is still some kind of layering algorithm.</dd>
 * </dl>
 */
public class HighDegreeNodeLayeringProcessor implements ILayoutProcessor<LGraph> {

    private static final Function<LNode, Iterable<LEdge>> INCOMING_EDGES = n -> n.getIncomingEdges();
    private static final Function<LNode, Iterable<LEdge>> OUTGOING_EDGES = n -> n.getOutgoingEdges();
    
    private LGraph layeredGraph;
    
    /** Minimum degree for a node to be considered 'high degree'. */
    private int degreeThreshold;
    /** The maximum height of a tree to be considered. Trees with larger height are neglected. */
    private int treeHeightThreshold;
    
    @Override
    public void process(final LGraph graph, final IElkProgressMonitor progressMonitor) {
  
        this.layeredGraph = graph;
        
        // retrieve some properties
        degreeThreshold = graph.getProperty(LayeredOptions.HIGH_DEGREE_NODES_THRESHOLD);
        treeHeightThreshold = graph.getProperty(LayeredOptions.HIGH_DEGREE_NODES_TREE_HEIGHT);
        // translate 0 to 'arbitrary height'
        if (treeHeightThreshold == 0) {
            treeHeightThreshold = Integer.MAX_VALUE;
        }

        // now iterate through all layer
        final ListIterator<Layer> layerIt = graph.getLayers().listIterator();
        while (layerIt.hasNext()) {
            Layer lay = layerIt.next();

            // -----------------------------------------------------------------------
            //  #1 find high degree nodes and find their incoming and outgoing trees
            // -----------------------------------------------------------------------
            List<Pair<LNode, HighDegreeNodeInformation>> highDegreeNodes = new ArrayList<>();
            int incMax = -1;
            int outMax = -1;
            for (LNode n : lay.getNodes()) {
                if (isHighDegreeNode(n)) {
                    HighDegreeNodeInformation hdni = calculateInformation(n);
                    incMax = Math.max(incMax, hdni.incTreesMaxHeight);
                    outMax = Math.max(outMax, hdni.outTreesMaxHeight);
                    highDegreeNodes.add(Pair.of(n, hdni));
                }
            }

            // -----------------------------------------------------------------------
            //  #2 insert layers before the current layer and move the trees
            // -----------------------------------------------------------------------
            List<Layer> preLayers = new ArrayList<>();
            for (int i = 0; i < incMax; ++i) {
                preLayers.add(0, prependLayer(layerIt));
            }
            for (Pair<LNode, HighDegreeNodeInformation> highDegreeNode : highDegreeNodes) {
                List<LNode> incRoots = highDegreeNode.getSecond().incTreeRoots;

                if (incRoots == null) {
                    continue;
                }

                for (LNode incRoot : incRoots) {
                    moveTree(incRoot, INCOMING_EDGES, preLayers);
                }
            }

            // -----------------------------------------------------------------------
            //  #2 insert layers after the current layer and move the trees
            // -----------------------------------------------------------------------
            List<Layer> afterLayers = new ArrayList<>();
            for (int i = 0; i < outMax; ++i) {
                afterLayers.add(appendLayer(layerIt));
            }
            for (Pair<LNode, HighDegreeNodeInformation> highDegreeNode : highDegreeNodes) {
                List<LNode> outRoots = highDegreeNode.getSecond().outTreeRoots;

                if (outRoots == null) {
                    continue;
                }

                for (LNode outRoot : outRoots) {
                    moveTree(outRoot, OUTGOING_EDGES, afterLayers);
                }
            }

        }

        // -----------------------------------------------------------------------
        //  #2 it can happen that layers became empty 
        //     after the previous processing, remove them
        // -----------------------------------------------------------------------
        ListIterator<Layer> layerIt2 = graph.getLayers().listIterator();
        while (layerIt2.hasNext()) {
            Layer l = layerIt2.next();
            if (l.getNodes().isEmpty()) {
                layerIt2.remove();
            }
        }
    }
    
    private boolean isHighDegreeNode(final LNode node) {
        return degree(node) >= degreeThreshold;
    }
    
    /**
     * @return the calculated trees for the passed node. That is a list of incoming trees and a list
     *         of outgoing trees and the maximum height of both of them.
     */
    private HighDegreeNodeInformation calculateInformation(final LNode hdn) {
        
        HighDegreeNodeInformation hdni = new HighDegreeNodeInformation();
        
        // check for incoming trees
        for (LEdge incEdge : hdn.getIncomingEdges()) {

            if (incEdge.isSelfLoop()) {
                continue;
            }
            
            LNode src = incEdge.getSource().getNode();
            if (hasSingleConnection(src, OUTGOING_EDGES)) {

                int treeHeight =
                        isTreeRoot(src, OUTGOING_EDGES, INCOMING_EDGES);
                if (treeHeight == -1) {
                    continue;
                }

                // update the maximum tree height
                hdni.incTreesMaxHeight = Math.max(hdni.incTreesMaxHeight, treeHeight);

                // remember the tree root for later processing
                if (hdni.incTreeRoots == null) {
                    hdni.incTreeRoots = new ArrayList<>();
                }
                hdni.incTreeRoots.add(src);
            }

        }

        // outgoing trees
        for (LEdge outEdge : hdn.getOutgoingEdges()) {

            if (outEdge.isSelfLoop()) {
                continue;
            }
            
            LNode tgt = outEdge.getTarget().getNode();
            if (hasSingleConnection(tgt, INCOMING_EDGES)) {
                
                int treeHeight =
                        isTreeRoot(tgt, INCOMING_EDGES, OUTGOING_EDGES);

                if (treeHeight == -1) {
                    continue;
                }
                hdni.outTreesMaxHeight = Math.max(hdni.outTreesMaxHeight, treeHeight);

                if (hdni.outTreeRoots == null) {
                    hdni.outTreeRoots = new ArrayList<>();
                }
                hdni.outTreeRoots.add(tgt);
            }
        }
        
        return hdni;
    }
    
    /**
     * @return The newly added layer. The iterator itself is left where it was before.
     */
    private Layer prependLayer(final ListIterator<Layer> layerIt) {
        assert layerIt.hasPrevious();
        layerIt.previous();
        
        Layer l = new Layer(layeredGraph);
        layerIt.add(l);
        layerIt.next();
        
        return l;
    }
    
    /**
     * @return The newly added layer. The iterator is left where it was before.
     */
    private Layer appendLayer(final ListIterator<Layer> layerIt) {
        Layer l = new Layer(layeredGraph);
        layerIt.add(l);
        
        return l;
    }
    
    /**
     * Moves the tree defined by {@code root} and {@code edgesFun} to the {@code layers}. The layers
     * should be ordered according to edgesFun. In other words, if the tree expands to the left and
     * the tree's root should be inserted into l_i, the list of layers should be [l_i, l_i-1, l_i-2,
     * ...].
     */
    private void moveTree(final LNode root, 
            final Function<LNode, Iterable<LEdge>> edgesFun,
            final List<Layer> layers) {
        assert layers.size() > 0;
        
        root.setLayer(layers.get(0));
        
        List<Layer> subList = layers.subList(1, layers.size());
        
        for (LEdge e : edgesFun.apply(root)) {
            LNode other = other(e, root);
            moveTree(other, edgesFun, subList);
        }
        
    }
    
    private static int degree(final LNode node) {
        return degree(node, n -> n.getConnectedEdges());
    }
    
    private static int degree(final LNode node, final Function<LNode, Iterable<LEdge>> edgeSelector) {
        return ((int) java.util.stream.StreamSupport.stream(edgeSelector.apply(node).spliterator(), false).count());
    }

    /**
     * @param edgeSelector
     *            the edges to check
     * @return whether the passed edges connect node to a single other node, other than checking for
     *         degree(node) == 1, this allows for multiple edges between the two involved nodes.
     */
    private boolean hasSingleConnection(final LNode node,
            final Function<LNode, Iterable<LEdge>> edgeSelector) {

        LNode connection = null;

        for (LEdge e : edgeSelector.apply(node)) {
            if (connection == null) {
                connection = other(e, node);
            } else {
                if (other(e, node) != connection) {
                    return false;
                }
            }
        }
        
        return true;
    }
    
    /**
     * @return for an edge (u, v) returns u if node == v and v if node == u.
     */
    private LNode other(final LEdge edge, final LNode node) {
        if (edge.getSource().getNode() == node) {
            return edge.getTarget().getNode();
        } else {
            return edge.getSource().getNode();
        }
    }
    
    /**
     * Checks whether {@code root} is the root node of a subtree with height at most {@link #treeHeightThreshold}. 
     * This is the case if: (1) root is connected to at most one node, i.e. it can be detached from a larger graph
     * by removing exactly one edge, (2) root is either a leaf, or all of the connected nodes are subtrees again.
     * 
     * Note that the directions of {@code upEdges} and {@code downEdges} are not related to the
     * actual directions in the layered graph. They refer to the conceptual direction when imagining
     * a graph with its ancestors and descendants.
     * 
     * @param root
     *            the node to be tested
     * @param ancestorEdges
     *            the "upward" pointing edges in a tree, i.e. the ancestor. For a valid tree root,
     *            all of these edges must connect to a unique ancestor.
     * @param descendantEdges
     *            the "downward" pointing edges in a tree, i.e. the descendants.
     * @return either the height of the subtree of which {@code root} is the root, or -1 if it's not
     *         part of a tree.
     */
    private int isTreeRoot(final LNode root,
            final Function<LNode, Iterable<LEdge>> ancestorEdges,
            final Function<LNode, Iterable<LEdge>> descendantEdges) {
        
        // exclude the high degree nodes themselves from the 'tree moving' process,
        //  i.e. leave them where they are
        if (isHighDegreeNode(root)) {
            return -1;
        }
        
        // is it a proper (sub-)root, i.e. does the node have exactly 
        // one parent?
        if (!hasSingleConnection(root, ancestorEdges)) {
            return -1;
        }
        
        // is it a leaf?
        if (!descendantEdges.apply(root).iterator().hasNext()) {
            return 1;
        }
        
        // recursively check subtrees
        int currentHeight = 0;
        for (LEdge e : descendantEdges.apply(root)) {
            LNode other = other(e, root);
            int height = isTreeRoot(other, ancestorEdges, descendantEdges);
            
            // other is not the root of a subtree
            if (height == -1) {
                return -1;
            }
            
            currentHeight = Math.max(currentHeight, height);
            
            // if we hit the threshold ... stop
            if (currentHeight > treeHeightThreshold - 1) {
                return -1;
            }
        }
        
        return currentHeight + 1;
    }
    
    /**
     * Stores the incoming and outgoing tree roots for nodes 
     * with high degrees. 
     */
    private class HighDegreeNodeInformation {
        // SUPPRESS CHECKSTYLE NEXT 4 VisibilityModifier 
        int incTreesMaxHeight = -1;
        List<LNode> incTreeRoots;
        int outTreesMaxHeight = -1;
        List<LNode> outTreeRoots;
    }
    
}
