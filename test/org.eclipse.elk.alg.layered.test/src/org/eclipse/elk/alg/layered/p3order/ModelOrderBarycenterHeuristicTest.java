/*******************************************************************************
 * Copyright (c) 2022 Kiel University and others.
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
package org.eclipse.elk.alg.layered.p3order;

import static org.junit.Assert.assertTrue;

import java.util.List;

import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.options.CrossingMinimizationStrategy;
import org.eclipse.elk.alg.layered.options.GreedySwitchType;
import org.eclipse.elk.alg.layered.options.InternalProperties;
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
public class ModelOrderBarycenterHeuristicTest {

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
    public LayoutConfigurator preferEdgesWeightedConfigurator() {
        LayoutConfigurator config = new LayoutConfigurator();
        config.configure(ElkNode.class).setProperty(
                LayeredOptions.CROSSING_MINIMIZATION_STRATEGY,
                CrossingMinimizationStrategy.LAYER_SWEEP);
        config.configure(ElkNode.class).setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY,
                OrderingStrategy.NODES_AND_EDGES);
        config.configure(ElkNode.class).setProperty(LayeredOptions.CROSSING_MINIMIZATION_FORCE_NODE_MODEL_ORDER,
                true);

        config.configure(ElkNode.class).setProperty(LayeredOptions.CROSSING_MINIMIZATION_GREEDY_SWITCH_TYPE,
                GreedySwitchType.OFF);
        return config;
    }    

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Tests
    
    // Just check for errors that might occur

    @TestAfterProcessor(LayerSweepCrossingMinimizer.class)
    public void testModelOrderRespected(final Object graph) {
        for (Layer layer : (LGraph) graph) {
            // We iterate over the layer's nodes check whether the real nodes (the one with a model order)
            // are correctly ordered.
            int modelOrder = -1;
            for (LNode node : layer) {
                if (node.hasProperty(InternalProperties.MODEL_ORDER)) {
                    int newModelOrder = node.getProperty(InternalProperties.MODEL_ORDER);
                    assertTrue("Node " + node + " has model order " + newModelOrder
                                    + ", which is smaller than the previous model order of " + modelOrder,
                            newModelOrder > modelOrder);
                    modelOrder = newModelOrder;
                }
            }
        }
    }
}