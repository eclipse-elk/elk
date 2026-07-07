/*******************************************************************************
 * Copyright (c) 2020 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.

 * This Source Code may also be made available under the following Secondary
 * Licenses when the conditions for such availability set forth in the Eclipse
 * Public License v. 2.0 are satisfied: GPL-3.0 which is available at
 * https://www.gnu.org/licenses/gpl-3.0-standalone.html.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.layered.intermediate;

import static org.junit.Assert.assertTrue;

import java.util.Collection;
import java.util.Deque;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.OptionalDouble;
import java.util.stream.Collectors;

import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.test.framework.LayoutTestRunner;
import org.eclipse.elk.alg.test.framework.annotations.Algorithm;
import org.eclipse.elk.alg.test.framework.annotations.DefaultConfiguration;
import org.eclipse.elk.alg.test.framework.annotations.GraphResourceProvider;
import org.eclipse.elk.alg.test.framework.io.AbstractResourcePath;
import org.eclipse.elk.alg.test.framework.io.ModelResourcePath;
import org.eclipse.elk.graph.ElkNode;
import org.junit.Test;
import org.junit.runner.RunWith;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;

@RunWith(LayoutTestRunner.class)
@Algorithm(LayeredOptions.ALGORITHM_ID)
@DefaultConfiguration()
public class PartitionTest {

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Sources

    @GraphResourceProvider
    public List<AbstractResourcePath> testGraphs() {
        return Lists.newArrayList(
                new ModelResourcePath("tests/layered/partitioning/**/"),
                new ModelResourcePath("tickets/layered/656_brokenLayoutPartitions.elkt"));
    }


    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Tests

    /**
     * Each node needs to be placed to the left of nodes with higher partition IDs.
     */
    @Test
    public void testPartitionOrder(final ElkNode node) {
        Deque<ElkNode> graphQueue = new LinkedList<>(); 
        graphQueue.add(node);
        
        while (!graphQueue.isEmpty()) {
            ElkNode graph = graphQueue.poll();
            testPartitionOrderInGraph(graph);
            
            // Look for child nodes that also contain graphs
            graph.getChildren().stream()
                .filter(child -> node.isHierarchical())
                .forEach(child -> graphQueue.add(child));
        }
    }
    
    private void testPartitionOrderInGraph(final ElkNode graph) {
        // If partitioning is not even active, don't bother
        if (!graph.getProperty(LayeredOptions.PARTITIONING_ACTIVATE)) {
            return;
        }
        
        // Collect nodes which have a partition set
        Multimap<Integer, ElkNode> partitionToNodesMap = HashMultimap.create();
        
        // Go through all layerless nodes and collect all partition IDs in use
        graph.getChildren().stream()
            .filter(node -> node.hasProperty(LayeredOptions.PARTITIONING_PARTITION))
            .forEach(node -> partitionToNodesMap.put(node.getProperty(LayeredOptions.PARTITIONING_PARTITION), node));
        
        // If no nodes have partition IDs assigned to them, don't bother
        if (partitionToNodesMap.isEmpty()) {
            return;
        }
        
        // Collect a sorted list of partition IDs
        List<Integer> sortedPartitionIDs = partitionToNodesMap.keySet().stream()
                .sorted()
                .collect(Collectors.toList());
        
        // For each pair of consecutive partition IDs (a, b), connect all nodes from partition a to partition b
        Iterator<Integer> idIterator = sortedPartitionIDs.iterator();
        
        Integer firstId = idIterator.next();
        while (idIterator.hasNext()) {
            Integer secondId = idIterator.next();
            testPartitionCoordinates(partitionToNodesMap.get(firstId), partitionToNodesMap.get(secondId));
            firstId = secondId;
        }
    }

    private void testPartitionCoordinates(final Collection<ElkNode> firstPartition,
            final Collection<ElkNode> secondPartition) {
        
        // The rightmost coordinate of any node in the first partition...
        OptionalDouble rightmost = firstPartition.stream()
            .mapToDouble(node -> node.getX() + node.getWidth())
            .max();
        
        // ...must be to the left of all coordinates of the second partition
        OptionalDouble leftmost = secondPartition.stream()
                .mapToDouble(node -> node.getX())
                .min();
        
        if (rightmost.isPresent() && leftmost.isPresent()) {
            assertTrue(rightmost.getAsDouble() < leftmost.getAsDouble());
        }
    }

}
