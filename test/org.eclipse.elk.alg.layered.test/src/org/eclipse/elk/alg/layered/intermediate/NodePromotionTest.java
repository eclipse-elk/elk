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

import java.util.List;

import org.eclipse.elk.alg.layered.options.CycleBreakingStrategy;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.layered.options.LayeringStrategy;
import org.eclipse.elk.alg.layered.options.NodePromotionStrategy;
import org.eclipse.elk.alg.layered.options.OrderingStrategy;
import org.eclipse.elk.alg.layered.p2layers.BasicLayerAssignmentTest;
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


/**
 * Basic tests for the node promotion heuristic. The test uses the {@link BasicLayerAssignmentTest} methods.
 */
@RunWith(LayoutTestRunner.class)
@Algorithm(LayeredOptions.ALGORITHM_ID)
@DefaultConfiguration()
public class NodePromotionTest {
    
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
    public LayoutConfigurator nikolovConfigurator() {
        return configuratorFor(NodePromotionStrategy.NIKOLOV);
    }
    
    @ConfiguratorProvider
    public LayoutConfigurator nikolovImprovedConfigurator() {
        return configuratorFor(NodePromotionStrategy.NIKOLOV_IMPROVED);
    }
    
    @ConfiguratorProvider
    public LayoutConfigurator nikolovPixelConfigurator() {
        return configuratorFor(NodePromotionStrategy.NIKOLOV_PIXEL);
    }
    
    @ConfiguratorProvider
    public LayoutConfigurator nikolovImprovedPixelConfigurator() {
        return configuratorFor(NodePromotionStrategy.NIKOLOV_IMPROVED_PIXEL);
    }
    
    @ConfiguratorProvider
    public LayoutConfigurator NoBoundaryConfigurator() {
        return configuratorFor(NodePromotionStrategy.NO_BOUNDARY);
    }
    
    @ConfiguratorProvider
    public LayoutConfigurator modelOrderLeftToRight() {
        LayoutConfigurator config = configuratorFor(NodePromotionStrategy.MODEL_ORDER_LEFT_TO_RIGHT);
        config.configure(ElkNode.class).setProperty(LayeredOptions.CYCLE_BREAKING_STRATEGY,
                CycleBreakingStrategy.MODEL_ORDER);
        config.configure(ElkNode.class).setProperty(LayeredOptions.LAYERING_STRATEGY,
                LayeringStrategy.LONGEST_PATH_SOURCE);
        config.configure(ElkNode.class).setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY,
                OrderingStrategy.PREFER_EDGES);
        return config;
    }
    
    @ConfiguratorProvider
    public LayoutConfigurator modelOrderRightToLeft() {
        LayoutConfigurator config = configuratorFor(NodePromotionStrategy.MODEL_ORDER_RIGHT_TO_LEFT);
        config.configure(ElkNode.class).setProperty(LayeredOptions.CYCLE_BREAKING_STRATEGY,
                CycleBreakingStrategy.MODEL_ORDER);
        config.configure(ElkNode.class).setProperty(LayeredOptions.LAYERING_STRATEGY,
                LayeringStrategy.LONGEST_PATH);
        config.configure(ElkNode.class).setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY,
                OrderingStrategy.PREFER_EDGES);
        return config;
    }
    
    private LayoutConfigurator configuratorFor(final NodePromotionStrategy strategy) {
        LayoutConfigurator config = new LayoutConfigurator();
        config.configure(ElkNode.class).setProperty(LayeredOptions.LAYERING_NODE_PROMOTION_STRATEGY, strategy);
        return config;
    }
    

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Tests
    
    @TestAfterProcessor(NodePromotion.class)
    public void testNoLayerlessNodes(final Object graph) {
        new BasicLayerAssignmentTest().testNoLayerlessNodes(graph);
    }

    @TestAfterProcessor(NodePromotion.class)
    public void testNoEmptyLayers(final Object graph) {
        new BasicLayerAssignmentTest().testNoEmptyLayers(graph);
    }

    @TestAfterProcessor(NodePromotion.class)
    public void testEdgesPointTowardsNextLayers(final Object graph) {
        new BasicLayerAssignmentTest().testEdgesPointTowardsNextLayers(graph);
    }
    
}