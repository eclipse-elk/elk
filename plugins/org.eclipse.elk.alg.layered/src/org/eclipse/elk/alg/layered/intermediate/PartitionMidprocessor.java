/*******************************************************************************
 * Copyright (c) 2016, 2020 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.alg.layered.intermediate;

import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.stream.Collectors;

import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LNode.NodeType;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.core.alg.ILayoutProcessor;
import org.eclipse.elk.core.options.PortSide;
import org.eclipse.elk.core.util.IElkProgressMonitor;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;

/**
 * Add constraint edges between partitions that force layering algorithms to adhere to the partitions.
 * 
 * <p>
 * In partitioned graphs, nodes may have an assigned partition. For every such node it must hold that every other node
 * with a smaller partition is placed left to it. This is accomplished by adding partition constraint edges. For each
 * node, we add an edge to every node in the next greater partition. During layering, the constraint edges assure the
 * aforementioned condition. The constraint edges are removed later by the {@link PartitionPostprocessor}.
 * </p>
 *
 * <dl>
 *   <dt>Precondition:</dt>
 *     <dd>an unlayered, acyclic graph.</dd>
 *   <dt>Postcondition:</dt>
 *     <dd>unlayered graph with partition constraint edges.</dd>
 *   <dt>Slots:</dt>
 *     <dd>Before phase 2.</dd>
 *   <dt>Same-slot dependencies:</dt>
 *     <dd>None.</dd>
 * </dl>
 * 
 * @see PartitionPreprocessor
 * @see PartitionPostprocessor
 */
public class PartitionMidprocessor implements ILayoutProcessor<LGraph> {

    @Override
    public void process(final LGraph lGraph, final IElkProgressMonitor monitor) {
        monitor.begin("Partition midprocessing", 1);
        
        // Collect nodes which have a partition set
        Multimap<Integer, LNode> partitionToNodesMap = HashMultimap.create();
        
        // Go through all layerless nodes and collect all partition IDs in use
        lGraph.getLayerlessNodes().stream()
            .filter(node -> node.hasProperty(LayeredOptions.PARTITIONING_PARTITION))
            .forEach(node -> partitionToNodesMap.put(node.getProperty(LayeredOptions.PARTITIONING_PARTITION), node));
        
        if (partitionToNodesMap.isEmpty()) {
            // This shouldn't happen, but if it does, there's nothing to do
            return;
        }
        
        // Collect a sorted list of partition IDs
        List<Integer> sortedPartitionIDs = partitionToNodesMap.keySet().stream()
            .sorted()
            .collect(Collectors.toList());
        
        // For each pair of consecutive partition IDs (a, b), connect all nodes from partition a to partition b; we use
        // Integer here to avoid unnecessary unboxing.
        Iterator<Integer> idIterator = sortedPartitionIDs.iterator();
        
        Integer firstId = idIterator.next();
        while (idIterator.hasNext()) {
            Integer secondId = idIterator.next();
            connectNodes(partitionToNodesMap.get(firstId), partitionToNodesMap.get(secondId));
            firstId = secondId;
        }
        
        monitor.done();
    }

    /**
     * Connects all nodes from the first collection to all nodes from the second collection.
     */
    private void connectNodes(final Collection<LNode> firstPartition, final Collection<LNode> secondPartition) {
        for (LNode node : firstPartition) {
            LPort sourcePort = new LPort();
            sourcePort.setNode(node);
            sourcePort.setSide(PortSide.EAST);
            sourcePort.setProperty(InternalProperties.PARTITION_DUMMY, true);
            boolean labelConnected = false;
            
            for (LEdge edge : node.getOutgoingEdges()) {
                if (edge.getTarget().getNode().getType() == NodeType.LABEL) {
                    labelConnected = true;
                }
            }
            
            for (LNode otherNode : secondPartition) {
                LPort targetPort = new LPort();
                targetPort.setNode(otherNode);
                targetPort.setSide(PortSide.WEST);
                targetPort.setProperty(InternalProperties.PARTITION_DUMMY, true);
                
                LEdge edge = new LEdge();
                edge.setProperty(InternalProperties.PARTITION_DUMMY, true);
                edge.setSource(sourcePort);
                edge.setTarget(targetPort);
                
                if (labelConnected) {
                    // Dummy node to compensate for potential edge label and its influence on the layering.
                    LNode dummyNode = new LNode(node.getGraph());
                    node.getGraph().getLayerlessNodes().add(dummyNode);
                    dummyNode.setType(NodeType.LONG_EDGE);
                    dummyNode.setProperty(InternalProperties.PARTITION_DUMMY, true);
                    LPort dummyIn = new LPort();
                    dummyIn.setNode(dummyNode);
                    dummyIn.setSide(PortSide.WEST);
                    LPort dummyOut = new LPort();
                    dummyOut.setNode(dummyNode);
                    dummyOut.setSide(PortSide.EAST);
                    
                    edge.setTarget(dummyIn);
                    
                    LEdge edge2 = new LEdge();
                    edge2.setProperty(InternalProperties.PARTITION_DUMMY, true);
                    edge2.setSource(dummyOut);
                    edge2.setTarget(targetPort);
                }
            }
        }
    }

}
