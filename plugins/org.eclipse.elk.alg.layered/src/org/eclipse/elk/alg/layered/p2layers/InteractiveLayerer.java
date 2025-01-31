/*******************************************************************************
 * Copyright (c) 2011, 2020 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.alg.layered.p2layers;

import java.util.LinkedHashSet;
import java.util.List;
import java.util.ListIterator;
import org.eclipse.elk.alg.layered.LayeredPhases;
import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.intermediate.IntermediateProcessorStrategy;
import org.eclipse.elk.alg.layered.options.PortType;
import org.eclipse.elk.core.alg.ILayoutPhase;
import org.eclipse.elk.core.alg.LayoutProcessorConfiguration;
import org.eclipse.elk.core.util.IElkProgressMonitor;

import com.google.common.collect.Lists;

/**
 * A node layerer that allows user interaction by respecting previous node positions.
 * These positions could be contrary to edge directions, so the resulting layering must
 * be checked for consistency.
 * 
 * <dl>
 *   <dt>Precondition:</dt><dd>the graph has no cycles</dd>
 *   <dt>Postcondition:</dt><dd>all nodes have been assigned a layer such that
 *     edges connect only nodes from layers with increasing indices</dd>
 * </dl>
 */
public final class InteractiveLayerer implements ILayoutPhase<LayeredPhases, LGraph> {

    @Override
    public LayoutProcessorConfiguration<LayeredPhases, LGraph> getLayoutProcessorConfiguration(final LGraph graph) {
        return LayoutProcessorConfiguration.<LayeredPhases, LGraph>create()
                .addBefore(LayeredPhases.P1_CYCLE_BREAKING,
                        IntermediateProcessorStrategy.INTERACTIVE_EXTERNAL_PORT_POSITIONER)
                .addBefore(LayeredPhases.P2_LAYERING,
                        IntermediateProcessorStrategy.LAYER_CONSTRAINT_PREPROCESSOR)
                .addBefore(LayeredPhases.P3_NODE_ORDERING,
                        IntermediateProcessorStrategy.LAYER_CONSTRAINT_POSTPROCESSOR);
    }
    
    /** Utility class for marking horizontal regions that are already covered by some nodes. */
    private static class LayerSpan {
        private double start;
        private double end;
        private List<LNode> nodes = Lists.newArrayList();
    }

    @Override
    public void process(final LGraph layeredGraph, final IElkProgressMonitor monitor) {
        monitor.begin("Interactive node layering", 1);

        // create layers with a start and an end position, merging when they overlap with others
        List<LayerSpan> currentSpans = Lists.newArrayList();
        for (LNode node : layeredGraph.getLayerlessNodes()) {
            double minx = node.getPosition().x;
            double maxx = minx + node.getSize().x;
            // for the following code we have to guarantee that every node has a width,
            // which, for instance, might not be the case for external dummy nodes
            maxx = Math.max(minx + 1, maxx);
            // look for a position in the sorted list where the node can be inserted
            ListIterator<LayerSpan> spanIter = currentSpans.listIterator();
            LayerSpan foundSpan = null;
            while (spanIter.hasNext()) {
                LayerSpan span = spanIter.next();
                if (span.start >= maxx) {
                    // the next layer span is further right, so insert the node here
                    spanIter.previous();
                    break;
                } else if (span.end > minx) {
                    // the layer span has an intersection with the node
                    if (foundSpan == null) {
                        // add the node to the current layer span
                        span.nodes.add(node);
                        span.start = Math.min(span.start, minx);
                        span.end = Math.max(span.end, maxx);
                        foundSpan = span;
                    } else {
                        // merge the previously found layer span with the current one
                        foundSpan.nodes.addAll(span.nodes);
                        foundSpan.end = Math.max(foundSpan.end, span.end);
                        spanIter.remove();
                    }
                }
            }
            if (foundSpan == null) {
                // no intersecting span was found, so create a new one
                foundSpan = new LayerSpan();
                foundSpan.start = minx;
                foundSpan.end = maxx;
                spanIter.add(foundSpan);
                foundSpan.nodes.add(node);
            }
        }
        
        // create real layers from the layer spans
        List<Layer> layers = layeredGraph.getLayers();
        int nextIndex = 0;
        for (LayerSpan span : currentSpans) {
            Layer layer = new Layer(layeredGraph);
            layer.id = nextIndex++;
            layers.add(layer);
            for (LNode node : span.nodes) {
                node.setLayer(layer);
                node.id = 0;
            }
        }
        
        // correct the layering respecting the graph topology, so edges point from left to right
        for (LNode node : layeredGraph.getLayerlessNodes()) {
            if (node.id == 0) {
                LinkedHashSet<LNode> shiftedNodes = checkNode(node, layeredGraph);
                // Since shiftedNodes might require other nodes to be shifted, do it again until all shiftedNodes
                // do no longer require another shift.
                while (!shiftedNodes.isEmpty()) {
                    LNode nodeToCheck = shiftedNodes.iterator().next();
                    shiftedNodes.remove(nodeToCheck);
                    shiftedNodes.addAll(checkNode(nodeToCheck, layeredGraph));
                } 
            }
        }
        
        // remove empty layers, which can happen when the layering has to be corrected
        ListIterator<Layer> layerIterator = layers.listIterator();
        while (layerIterator.hasNext()) {
            if (layerIterator.next().getNodes().isEmpty()) {
                layerIterator.remove();
            }
        }
        
        // clear the list of nodes that have no layer, since now they all have one
        layeredGraph.getLayerlessNodes().clear();
        monitor.done();
    }
    
    /**
     * Check the layering of the given node by comparing the layer index of all successors.
     * 
     * @param node1 a node
     * @param graph the layered graph
     */
    private LinkedHashSet<LNode> checkNode(final LNode node1, final LGraph graph) {
        node1.id = 1;
        Layer layer1 = node1.getLayer();
        LinkedHashSet<LNode> shiftNodes = new LinkedHashSet<>();
        for (LPort port : node1.getPorts(PortType.OUTPUT)) {
            for (LEdge edge : port.getOutgoingEdges()) {
                LNode node2 = edge.getTarget().getNode();
                if (node1 != node2) {
                    Layer layer2 = node2.getLayer();
                    if (layer2.id <= layer1.id) {
                        // a violation was detected - move the target node to the next layer
                        int newIndex = layer1.id + 1;
                        if (newIndex == graph.getLayers().size()) {
                            Layer newLayer = new Layer(graph);
                            newLayer.id = newIndex;
                            graph.getLayers().add(newLayer);
                            node2.setLayer(newLayer);
                        } else {
                            Layer newLayer = graph.getLayers().get(newIndex);
                            node2.setLayer(newLayer);
                        }
                        shiftNodes.add(node2);
                    }
                }
            }
        }
        return shiftNodes;
    }
}
