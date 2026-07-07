/*******************************************************************************
 * Copyright (c) 2013, 2020 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.layered.intermediate;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;

import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LNode.NodeType;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.test.framework.LayoutTestRunner;
import org.eclipse.elk.alg.test.framework.annotations.Algorithm;
import org.eclipse.elk.alg.test.framework.annotations.DefaultConfiguration;
import org.eclipse.elk.alg.test.framework.annotations.GraphResourceProvider;
import org.eclipse.elk.alg.test.framework.annotations.TestAfterProcessor;
import org.eclipse.elk.alg.test.framework.io.AbstractResourcePath;
import org.eclipse.elk.alg.test.framework.io.ModelResourcePath;
import org.junit.runner.RunWith;

import com.google.common.collect.Lists;

@RunWith(LayoutTestRunner.class)
@Algorithm(LayeredOptions.ALGORITHM_ID)
@DefaultConfiguration()
public class PartitionPostProcessorTest {

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
     * After constraint processing no empty layers should be left behind.
     */
    @TestAfterProcessor(PartitionPostprocessor.class)
    public void noEmptyLayerTest(final Object graph) {
        LGraph lGraph = (LGraph) graph;
        
        for (Layer layer : lGraph) {
            assertFalse(layer.getNodes().isEmpty());
        }
    }

    /**
     * For each layer, every node in this layer has the same partition and and this partition is greater or equal to the
     * partition of the previous layer.
     */
    @TestAfterProcessor(PartitionPostprocessor.class)
    public void testPartitionOrder(final Object graph) {
        LGraph lGraph = (LGraph) graph;
        
        int lastPartition = -1;
        for (Layer layer : lGraph) {
            // Set to invalid partition
            int currentPartition = -1;
            
            for (LNode node : layer) {
                if (node.getType() == NodeType.NORMAL && node.hasProperty(LayeredOptions.PARTITIONING_PARTITION)) {
                    int nodePartition = node.getProperty(LayeredOptions.PARTITIONING_PARTITION);
                    if (currentPartition == -1) {
                        currentPartition = nodePartition;
                    } else {
                        assertEquals(currentPartition, nodePartition);
                    }
                }
            }
            
            assertTrue(lastPartition <= currentPartition);
            lastPartition = currentPartition;
        }
    }

}
