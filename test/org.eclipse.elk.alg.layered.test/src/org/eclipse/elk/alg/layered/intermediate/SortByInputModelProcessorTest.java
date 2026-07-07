/*******************************************************************************
 * Copyright (c) 2020 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 * 
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later 
 *******************************************************************************/
package org.eclipse.elk.alg.layered.intermediate;

import static org.junit.Assert.assertEquals;

import java.util.List;

import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LNode.NodeType;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.options.CrossingMinimizationStrategy;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.layered.options.OrderingStrategy;
import org.eclipse.elk.alg.layered.p3order.LayerSweepCrossingMinimizer;
import org.eclipse.elk.alg.test.framework.LayoutTestRunner;
import org.eclipse.elk.alg.test.framework.annotations.Algorithm;
import org.eclipse.elk.alg.test.framework.annotations.ConfiguratorProvider;
import org.eclipse.elk.alg.test.framework.annotations.DefaultConfiguration;
import org.eclipse.elk.alg.test.framework.annotations.GraphResourceProvider;
import org.eclipse.elk.alg.test.framework.annotations.TestAfterProcessor;
import org.eclipse.elk.alg.test.framework.io.AbstractResourcePath;
import org.eclipse.elk.alg.test.framework.io.ModelResourcePath;
import org.eclipse.elk.core.LayoutConfigurator;
import org.eclipse.elk.graph.ElkNode;
import org.junit.runner.RunWith;

import com.google.common.collect.Lists;

@RunWith(LayoutTestRunner.class)
@Algorithm(LayeredOptions.ALGORITHM_ID)
@DefaultConfiguration()
public class SortByInputModelProcessorTest {

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Sources

    @GraphResourceProvider
    public List<AbstractResourcePath> testGraphs() {
        return Lists.newArrayList(new ModelResourcePath("tests/layered/preserve_order/**"));
    }
    

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Configuration

    @ConfiguratorProvider
    public LayoutConfigurator preserveOrderConfigurator() {
        LayoutConfigurator config = new LayoutConfigurator();
        config.configure(ElkNode.class).setProperty(
                LayeredOptions.CROSSING_MINIMIZATION_STRATEGY,
                CrossingMinimizationStrategy.LAYER_SWEEP);
        config.configure(ElkNode.class).setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY, OrderingStrategy.NODES_AND_EDGES);
        return config;
    }
    

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Tests

    @TestAfterProcessor(LayerSweepCrossingMinimizer.class)
    public void preserveNodeAndEdgeOrder(final Object graph) {
        LGraph lgraph = (LGraph) graph;
        int layerIndex = 0;
        for (Layer layer : lgraph) {
            int nodeIndex = 0;
            for (LNode node : layer.getNodes()) {
                if (node.getType() == NodeType.NORMAL) {
                    assertEquals("A node is on an unexpected position.",
                            "n_l" + layerIndex + "p" + nodeIndex, node.toString());
                }
                nodeIndex++;
            }
            layerIndex++;
        }
    }
}
