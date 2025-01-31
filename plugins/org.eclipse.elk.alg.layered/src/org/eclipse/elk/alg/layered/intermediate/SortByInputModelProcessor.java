/*******************************************************************************
 * Copyright (c) 2020 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 * 
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.alg.layered.intermediate;

import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LNode.NodeType;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.intermediate.preserveorder.ModelOrderNodeComparator;
import org.eclipse.elk.alg.layered.intermediate.preserveorder.ModelOrderPortComparator;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.core.alg.ILayoutProcessor;
import org.eclipse.elk.core.options.PortConstraints;
import org.eclipse.elk.core.util.IElkProgressMonitor;

/**
 * Sorts all child nodes and their edges of the given graph by their {@link InternalProperties#MODEL_ORDER},
 * which represents the order in input graph.
 * Outgoing ports are sorted by the order of their edges, incoming ports are sorted by the order of their nodes in the
 * previous layer.
 * <dl><dl>
 *   <dt>Precondition:</dt>
 *     <dd>a proper layered graph</dd>
 *     <dd>nodes have fixed port sides.</dd>
 *   <dt>Postcondition:</dt>
 *     <dd>Nodes and ports (edges) are sorted to respect the order in the input graph</dd>
 *   <dt>Slots:</dt>
 *     <dd>Before phase 3.</dd>
 *   <dt>Same-slot dependencies:</dt>
 *     <dd>Since this processor assumes a proper layered graph it has to be executed after
 *     the {@link LongEdgeSplitter}</dd>
 *     <dd>Since the {@link NorthSouthPortPreprocessor} introduces same layer edges it has to be executed after
 *     this processor</dd>
 * </dl>
 */
public class SortByInputModelProcessor implements ILayoutProcessor<LGraph> {

    @Override
    public void process(final LGraph graph, final IElkProgressMonitor progressMonitor) {
        progressMonitor.begin("Sort By Input Model "
                + graph.getProperty(LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY), 1);
        int layerIndex = 0;
        for (Layer layer : graph) {
            // The layer.id is necessary to check whether nodes really connect to the previous layer.
            // Feedback long edge dummies do might have a long-edge dummy in the same layer.
            layer.id = layerIndex;
            final int previousLayerIndex = layerIndex == 0 ? 0 : layerIndex - 1;
            Layer previousLayer = graph.getLayers().get(previousLayerIndex);
            // Sort nodes before port sorting to have sorted nodes for in-layer feedback edge dummies.
            ModelOrderNodeComparator comparator = new ModelOrderNodeComparator(previousLayer,
                    graph.getProperty(LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY),
                    graph.getProperty(LayeredOptions.CONSIDER_MODEL_ORDER_LONG_EDGE_STRATEGY), true);
            SortByInputModelProcessor.insertionSort(layer.getNodes(), comparator);
            for (LNode node : layer.getNodes()) {
                if (node.getProperty(LayeredOptions.PORT_CONSTRAINTS) != PortConstraints.FIXED_ORDER
                        && node.getProperty(LayeredOptions.PORT_CONSTRAINTS) != PortConstraints.FIXED_POS) {
                    // Special case:
                    // If two edges (of the same node) have the same target node they should be next to each other.
                    // Therefore all ports that connect to the same node should have the same
                    // (their minimal) model order.
                    // Get minimal model order for target node
                    Collections.sort(node.getPorts(),
                            new ModelOrderPortComparator(previousLayer,
                                    graph.getProperty(LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY),
                                    longEdgeTargetNodePreprocessing(node),
                                    graph.getProperty(LayeredOptions.CONSIDER_MODEL_ORDER_PORT_MODEL_ORDER)));
                    progressMonitor.log("Node " + node + " ports: " + node.getPorts());
                }
            }
            // Sort nodes after port sorting to also sort dummy feedback nodes from the current layer correctly.
            comparator = new ModelOrderNodeComparator(previousLayer,
                    graph.getProperty(LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY),
                    graph.getProperty(LayeredOptions.CONSIDER_MODEL_ORDER_LONG_EDGE_STRATEGY), false);
            SortByInputModelProcessor.insertionSort(layer.getNodes(), comparator);
                    
            progressMonitor.log("Layer " + layerIndex + ": " + layer);
            layerIndex++;
        }
        progressMonitor.done();
    }
    
    /**
     * Calculate long edge target of port and saves it and adds a map of all long edge targets to the node.
     * @param node the node
     * @return A map of all long edge targets of a node
     */
    public static Map<LNode, Integer> longEdgeTargetNodePreprocessing(final LNode node) {
        Map<LNode, Integer> targetNodeModelOrder = new HashMap<>();
        if (node.hasProperty(InternalProperties.TARGET_NODE_MODEL_ORDER)) {
            return node.getProperty(InternalProperties.TARGET_NODE_MODEL_ORDER);
        }
        node.getPorts().stream().filter(p -> !p.getOutgoingEdges().isEmpty()).forEach(p -> {
            LNode targetNode = getTargetNode(p);
            p.setProperty(InternalProperties.LONG_EDGE_TARGET_NODE, targetNode);
            if (targetNode != null) {
                int previousOrder = Integer.MAX_VALUE;
                if (targetNodeModelOrder.containsKey(targetNode)) {
                    previousOrder = targetNodeModelOrder.get(targetNode);
                }
                LEdge edge = p.getOutgoingEdges().get(0);
                if (!edge.getProperty(InternalProperties.REVERSED)) {
                    targetNodeModelOrder.put(targetNode, Math
                            .min(edge.getProperty(InternalProperties.MODEL_ORDER), previousOrder));
                }
            }
        });
        node.setProperty(InternalProperties.TARGET_NODE_MODEL_ORDER, targetNodeModelOrder);
        return targetNodeModelOrder;
        
    }

    /**
     * Returns the target node of a port considering long edges.
     * @param port The port
     * @return the target node of the long edge connecting to the port or null if none exist.
     */
    public static LNode getTargetNode(final LPort port) {
        LNode node = null;
        LEdge edge = port.getOutgoingEdges().get(0);
        do {
            node = edge.getTarget().getNode();
            // If the dummy node has a target return it.
            if (node.hasProperty(InternalProperties.LONG_EDGE_TARGET)) {
                return node.getProperty(InternalProperties.LONG_EDGE_TARGET).getNode();
            }
            // It not the current node might be the target node or one has to iterate manually through to it.
            if (node.getType() != NodeType.NORMAL && node.getOutgoingEdges().iterator().hasNext()) {
                edge = node.getOutgoingEdges().iterator().next();
            } else if (node.getType() != NodeType.NORMAL) {
                return null;
            }
        } while (node != null && node.getType() != NodeType.NORMAL);
        return node;
    }
    
    public static void insertionSort(final List<LNode> layer,
            final ModelOrderNodeComparator comparator) {
        LNode temp;
        for (int i = 1; i < layer.size(); i++) {
            temp = layer.get(i);
            int j = i;
            while (j > 0 && comparator.compare(layer.get(j - 1), temp) > 0) {
                layer.set(j, layer.get(j - 1));
                j--;
            }
            layer.set(j, temp);
        }
        comparator.clearTransitiveOrdering();
    }
    
    public static void insertionSortPort(final List<LPort> layer,
            final ModelOrderPortComparator comparator) {
        LPort temp;
        for (int i = 1; i < layer.size(); i++) {
            temp = layer.get(i);
            int j = i;
            while (j > 0 && comparator.compare(layer.get(j - 1), temp) > 0) {
                layer.set(j, layer.get(j - 1));
                j--;
            }
            layer.set(j, temp);
        }
        comparator.clearTransitiveOrdering();
    }
}

