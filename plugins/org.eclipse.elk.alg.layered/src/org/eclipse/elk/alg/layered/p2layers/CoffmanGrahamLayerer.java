/*******************************************************************************
 * Copyright (c) 2016, 2020 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.alg.layered.p2layers;

import java.util.Arrays;
import java.util.List;
import java.util.ListIterator;
import java.util.PriorityQueue;

import org.eclipse.elk.alg.layered.LayeredPhases;
import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.intermediate.IntermediateProcessorStrategy;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.core.alg.ILayoutPhase;
import org.eclipse.elk.core.alg.LayoutProcessorConfiguration;
import org.eclipse.elk.core.util.IElkProgressMonitor;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;

/**
 * A layering algorithm that places nodes in layers subject to a bound on the maximum number of (original) nodes per 
 * layer. The bound is set via the {@link LayeredOptions#LAYERING_COFFMAN_GRAHAM_LAYER_BOUND} layout option.
 * 
 * The algorithm is described in 
 * <ul><li>
 *   Coffman, E. G., Jr. and Graham, R. L., "Optimal scheduling for two-processor systems", 
 *   <i>Acta Informatica</i> (1972).
 * </li></ul> 
 *  
 * <dl>
 *   <dt>Precondition:</dt><dd>the graph has no cycles.</dd>
 *   <dt>Postcondition:</dt><dd>all nodes have been assigned a layer such that
 *     edges connect only nodes from layers with increasing indices and the specified 
 *     layer bound is met.</dd>
 * </dl>
 */
public class CoffmanGrahamLayerer implements ILayoutPhase<LayeredPhases, LGraph> {

    /** Used during dfs of the transitive reduction algorithm. */
    private boolean[] nodeMark;
    /** A marked edge is one that is not part of the transitive reduction. */
    private boolean[] edgeMark;

    /** Used to identigy sources during computation of the topological ordering. */
    private int[] inDeg;
    /** Used to identify sinks during layer assignment. */
    private int[] outDeg;
    
    /** The calculated topological order which forms the basis for the layering. */
    private int[] topoOrd;
    
    /** For each node, stores the positions of incoming nodes in the topological ordering. 
     *  Due to the way the positions are added to the lists, one can assume that the 
     *  lists are sorted. */
    private ListMultimap<LNode, Integer> inTopo = ArrayListMultimap.create(); 
    
    @Override
    public void process(final LGraph layeredGraph, final IElkProgressMonitor progressMonitor) {
        progressMonitor.begin("Coffman-Graham Layering", 1);
        
        if (layeredGraph.getLayerlessNodes().isEmpty()) {
            progressMonitor.done();
            return;
        }
        
        // the maximum number of allowed nodes per layer
        int w = layeredGraph.getProperty(LayeredOptions.LAYERING_COFFMAN_GRAHAM_LAYER_BOUND);
        
        // initialization, assign indexes and initialize arrays
        int index = 0;
        int edgeIndex = 0;
        for (LNode n : layeredGraph.getLayerlessNodes()) {
            n.id = index++;
            for (LEdge e : n.getOutgoingEdges()) {
                e.id = edgeIndex++;
            }
        }
        nodeMark = new boolean[index];
        edgeMark = new boolean[edgeIndex];
        inDeg = new int[index];
        outDeg = new int[index];
        topoOrd = new int[index];
        inTopo.clear();
        
        // --------------------------
        // #1 Remove transitive edges
        // --------------------------
        transitiveReduction(layeredGraph);
        
        // -------------------------------
        // #2 Compute topological ordering
        // -------------------------------
        PriorityQueue<LNode> sources = new PriorityQueue<LNode>(this::compareNodesInTopo);

        // for each node, determine its current in-degree and remember initial sources
        for (LNode v : layeredGraph.getLayerlessNodes()) {
            for (LEdge e : v.getIncomingEdges()) {
                if (!edgeMark[e.id]) {
                    inDeg[v.id]++;
                }
            }
            if (inDeg[v.id] == 0) {
                sources.add(v);
            }
        }
        
        // compute topological ordering
        int i = 0;
        while (!sources.isEmpty()) {
            LNode v = sources.poll();
            // assign number of topological order
            topoOrd[v.id] = i++;

            // update the rest of the graph
            for (LEdge e : v.getOutgoingEdges()) {
                if (edgeMark[e.id]) {
                    continue; 
                }
                LNode tgt = e.getTarget().getNode();
                inDeg[tgt.id]--;
                inTopo.put(tgt, topoOrd[v.id]);
                if (inDeg[tgt.id] == 0) {
                    sources.add(tgt); // 'tgt' is added according to its priority
                }
            }
        }
        
        // --------------------------
        // #3 Actual layer assignment
        // --------------------------
        // note that this time we start with sinks and work our way to the original graph's sources
        
        // set of candidate nodes, that is nodes for which all successors have been placed
        // highest priority first
        PriorityQueue<LNode> sinks =
                new PriorityQueue<LNode>((n1, n2) -> -Integer.compare(topoOrd[n1.id], topoOrd[n2.id]));
        for (LNode v : layeredGraph.getLayerlessNodes()) {
            for (LEdge e : v.getOutgoingEdges()) {
                if (!edgeMark[e.id]) {
                    outDeg[v.id]++;
                }
            }
            if (outDeg[v.id] == 0) {
                sinks.add(v);
            }
        }

        // assign the layers
        List<Layer> layers = Lists.newArrayList();
        Layer currentLayer = createLayer(layeredGraph, layers);
        while (!sinks.isEmpty()) {
            // select a node for which all outgoing nodes have been placed,
            // and with maximum value in the topological sort
            LNode u = sinks.poll();

            // start a new layer if the current one is already full,
            // or an in-layer edge would be introduced
            if (isLayerFull(currentLayer, w) || !canAdd(u, currentLayer)) {
                currentLayer = createLayer(layeredGraph, layers);
            }

            // place the node in the layer
            u.setLayer(currentLayer);

            // update out-degrees and collect the new sinks
            for (LEdge e : u.getIncomingEdges()) {
                if (edgeMark[e.id]) {
                    continue;
                }
                LNode src = e.getSource().getNode();
                outDeg[src.id]--;
                if (outDeg[src.id] == 0) {
                    sinks.add(src);
                }
            }
        }
        
        // the layers were created in inverse order
        for (int j = layers.size() - 1; j >= 0; --j) {
            layeredGraph.getLayers().add(layers.get(j));
        }
        
        // clear layerless nodes
        layeredGraph.getLayerlessNodes().clear();
        
        progressMonitor.done();
    }
    
    private boolean isLayerFull(final Layer layer, final int w) {
        return layer.getNodes().size() >= w;
    }
    
    /**
     * @return true if node {@code n} can be added to layer {@code l} without introducing in-layer edges.
     */
    private boolean canAdd(final LNode n, final Layer l) {
        for (LEdge e : n.getOutgoingEdges()) {
            LNode v = e.getTarget().getNode();
            if (v.getLayer() == l) {
                return false;
            }
        }
        return true;
    }

    /**
     * @return a newly created layer for the passed graph. The layer is <b>appended</b> to the passed list of layers.
     *         Note that the layer is not yet added to the graph's list of layers.
     */
    private Layer createLayer(final LGraph graph, final List<Layer> layers) {
        Layer aLayer = new Layer(graph);
        layers.add(aLayer);
        return aLayer;
    }
    
    private int compareNodesInTopo(final LNode u, final LNode v) {
        List<Integer> inListU = inTopo.get(u);
        List<Integer> inListV = inTopo.get(v);
        ListIterator<Integer> itU = inListU.listIterator(inListU.size());
        ListIterator<Integer> itV = inListV.listIterator(inListV.size());
        
        // find the node with the lower associated maximum value in 'inTopo'
        // break ties by ignoring all equal maxima
        while (itU.hasPrevious() && itV.hasPrevious()) {
            Integer iu = itU.previous();
            Integer iv = itV.previous();
            if (iu != iv) {
                return Integer.compare(iu, iv);
            }
        }
        
        if (!itU.hasNext() && !itV.hasNext()) {
            // If the two nodes are the same, a secondary criterion is needed.
            // Else use the node id as a last resort.
            if (u.id < v.id) {
                return -1; // u < v
            } else if (u.id > v.id) {
                return 1; // u > v
            }
            // Else they are equal.
            return 0; // u == v
        } else if (!itU.hasNext()) {
            return -1; // u < v
        } else {
            return 1; // u > v
        }
        
    }
    
    /**
     * Remove transitive edges.
     */
    private void transitiveReduction(final LGraph graph) {
        for (LNode start : graph.getLayerlessNodes()) {
            Arrays.fill(nodeMark, false);
            for (LEdge out : start.getOutgoingEdges()) {
                dfs(start, out.getTarget().getNode());
            }
        }
    }
    
    private void dfs(final LNode start, final LNode v) {
        if (nodeMark[v.id]) {
            return;
        }
        
        for (LEdge out : v.getOutgoingEdges()) {
            LNode w = out.getTarget().getNode();
            for (LEdge transitive : w.getIncomingEdges()) {
                if (transitive.getSource().getNode() == start) {
                    edgeMark[transitive.id] = true;
                }
            }
            dfs(start, w);
        }
        
        nodeMark[v.id] = true;
    }
    
    /* ----------------------------------------------------------------------------------------------------------
     *                               Intermediate Processor Configuration 
     * ---------------------------------------------------------------------------------------------------------- */
    
    /** basic intermediate processing configuration. */
    private static final LayoutProcessorConfiguration<LayeredPhases, LGraph> BASELINE_PROCESSING_CONFIGURATION =
        LayoutProcessorConfiguration.<LayeredPhases, LGraph>create()
            .addBefore(LayeredPhases.P1_CYCLE_BREAKING,
                    IntermediateProcessorStrategy.EDGE_AND_LAYER_CONSTRAINT_EDGE_REVERSER)
            .addBefore(LayeredPhases.P2_LAYERING, IntermediateProcessorStrategy.LAYER_CONSTRAINT_PREPROCESSOR)
            .addBefore(LayeredPhases.P3_NODE_ORDERING, IntermediateProcessorStrategy.LAYER_CONSTRAINT_POSTPROCESSOR);
    
    @Override
    public LayoutProcessorConfiguration<LayeredPhases, LGraph> getLayoutProcessorConfiguration(final LGraph graph) {
        return BASELINE_PROCESSING_CONFIGURATION;
    }
    
}
