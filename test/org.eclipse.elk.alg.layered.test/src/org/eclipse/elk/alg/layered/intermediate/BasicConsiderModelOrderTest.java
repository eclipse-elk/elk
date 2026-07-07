/*******************************************************************************
 * Copyright (c) 2021 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 * 
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later 
 *******************************************************************************/
package org.eclipse.elk.alg.layered.intermediate;

import java.util.List;

import org.eclipse.elk.alg.layered.options.CrossingMinimizationStrategy;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.layered.options.OrderingStrategy;
import org.eclipse.elk.alg.test.framework.LayoutTestRunner;
import org.eclipse.elk.alg.test.framework.annotations.Algorithm;
import org.eclipse.elk.alg.test.framework.annotations.ConfiguratorProvider;
import org.eclipse.elk.alg.test.framework.annotations.DefaultConfiguration;
import org.eclipse.elk.alg.test.framework.annotations.GraphResourceProvider;
import org.eclipse.elk.alg.test.framework.annotations.TestAfterProcessor;
import org.eclipse.elk.alg.test.framework.io.AbstractResourcePath;
import org.eclipse.elk.alg.test.framework.io.FileExtensionFilter;
import org.eclipse.elk.alg.test.framework.io.ModelResourcePath;
import org.eclipse.elk.core.LayoutConfigurator;
import org.eclipse.elk.graph.ElkNode;
import org.junit.runner.RunWith;

import com.google.common.collect.Lists;

@RunWith(LayoutTestRunner.class)
@Algorithm(LayeredOptions.ALGORITHM_ID)
@DefaultConfiguration()
public class BasicConsiderModelOrderTest {

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Sources

    @GraphResourceProvider
    public List<AbstractResourcePath> testGraphs() {
        return Lists.newArrayList(
                new ModelResourcePath("realworld/ptolemy/**/").withFilter(new FileExtensionFilter("elkg")));
    }
    

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Configuration

    @ConfiguratorProvider
    public LayoutConfigurator nodesAndEdgesConfigurator() {
        LayoutConfigurator config = new LayoutConfigurator();
        config.configure(ElkNode.class).setProperty(
                LayeredOptions.CROSSING_MINIMIZATION_STRATEGY,
                CrossingMinimizationStrategy.LAYER_SWEEP);
        config.configure(ElkNode.class).setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY,
                OrderingStrategy.NODES_AND_EDGES);
        return config;
    }
    
    @ConfiguratorProvider
    public LayoutConfigurator nodesAndEdgesWeightedConfigurator() {
        LayoutConfigurator config = new LayoutConfigurator();
        config.configure(ElkNode.class).setProperty(
                LayeredOptions.CROSSING_MINIMIZATION_STRATEGY,
                CrossingMinimizationStrategy.LAYER_SWEEP);
        config.configure(ElkNode.class).setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY,
                OrderingStrategy.NODES_AND_EDGES);
        config.configure(ElkNode.class).setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_CROSSING_COUNTER_NODE_INFLUENCE,
                0.001);
        config.configure(ElkNode.class).setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_CROSSING_COUNTER_PORT_INFLUENCE,
                0.001);
        return config;
    }

    @ConfiguratorProvider
    public LayoutConfigurator preferEdgesConfigurator() {
        LayoutConfigurator config = new LayoutConfigurator();
        config.configure(ElkNode.class).setProperty(
                LayeredOptions.CROSSING_MINIMIZATION_STRATEGY,
                CrossingMinimizationStrategy.LAYER_SWEEP);
        config.configure(ElkNode.class).setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY,
                OrderingStrategy.PREFER_EDGES);
        return config;
    }

    @ConfiguratorProvider
    public LayoutConfigurator preferEdgesWeightedConfigurator() {
        LayoutConfigurator config = new LayoutConfigurator();
        config.configure(ElkNode.class).setProperty(
                LayeredOptions.CROSSING_MINIMIZATION_STRATEGY,
                CrossingMinimizationStrategy.LAYER_SWEEP);
        config.configure(ElkNode.class).setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY,
                OrderingStrategy.PREFER_EDGES);
        config.configure(ElkNode.class).setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_CROSSING_COUNTER_NODE_INFLUENCE,
                0.001);
        config.configure(ElkNode.class).setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_CROSSING_COUNTER_PORT_INFLUENCE,
                0.001);
        return config;
    }

    @ConfiguratorProvider
    public LayoutConfigurator preferNodesConfigurator() {
        LayoutConfigurator config = new LayoutConfigurator();
        config.configure(ElkNode.class).setProperty(
                LayeredOptions.CROSSING_MINIMIZATION_STRATEGY,
                CrossingMinimizationStrategy.LAYER_SWEEP);
        config.configure(ElkNode.class).setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY,
                OrderingStrategy.PREFER_NODES);
        return config;
    }

    @ConfiguratorProvider
    public LayoutConfigurator preferNodesWeightedConfigurator() {
        LayoutConfigurator config = new LayoutConfigurator();
        config.configure(ElkNode.class).setProperty(
                LayeredOptions.CROSSING_MINIMIZATION_STRATEGY,
                CrossingMinimizationStrategy.LAYER_SWEEP);
        config.configure(ElkNode.class).setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY,
                OrderingStrategy.PREFER_NODES);
        config.configure(ElkNode.class).setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_CROSSING_COUNTER_NODE_INFLUENCE,
                0.001);
        config.configure(ElkNode.class).setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_CROSSING_COUNTER_PORT_INFLUENCE,
                0.001);
        return config;
    }

    @ConfiguratorProvider
    public LayoutConfigurator portModelOrderPreferEdges() {
        LayoutConfigurator config = new LayoutConfigurator();
        config.configure(ElkNode.class).setProperty(
                LayeredOptions.CROSSING_MINIMIZATION_STRATEGY,
                CrossingMinimizationStrategy.LAYER_SWEEP);
        config.configure(ElkNode.class).setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY,
                OrderingStrategy.PREFER_EDGES);
        config.configure(ElkNode.class).setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_PORT_MODEL_ORDER,
                true);
        return config;
    }

    @ConfiguratorProvider
    public LayoutConfigurator portModelOrderNodesAndEdges() {
        LayoutConfigurator config = new LayoutConfigurator();
        config.configure(ElkNode.class).setProperty(
                LayeredOptions.CROSSING_MINIMIZATION_STRATEGY,
                CrossingMinimizationStrategy.LAYER_SWEEP);
        config.configure(ElkNode.class).setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY,
                OrderingStrategy.NODES_AND_EDGES);
        config.configure(ElkNode.class).setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_PORT_MODEL_ORDER,
                true);
        return config;
    }

    @ConfiguratorProvider
    public LayoutConfigurator portModelOrderPreferNodes() {
        LayoutConfigurator config = new LayoutConfigurator();
        config.configure(ElkNode.class).setProperty(
                LayeredOptions.CROSSING_MINIMIZATION_STRATEGY,
                CrossingMinimizationStrategy.LAYER_SWEEP);
        config.configure(ElkNode.class).setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY,
                OrderingStrategy.PREFER_NODES);
        config.configure(ElkNode.class).setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_PORT_MODEL_ORDER,
                true);
        return config;
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Tests
    
    // Just check for errors that might occur
    
    @TestAfterProcessor(SortByInputModelProcessor.class)
    public void test(final Object graph) {
        assert(true);
    }
}
