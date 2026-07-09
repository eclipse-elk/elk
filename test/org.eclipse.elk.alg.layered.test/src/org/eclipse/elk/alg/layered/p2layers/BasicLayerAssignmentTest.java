/*******************************************************************************
 * Copyright (c) 2013, 2020 Kiel University and others.
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
package org.eclipse.elk.alg.layered.p2layers;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.Deque;
import java.util.LinkedList;
import java.util.List;
import java.util.Random;

import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.layered.options.LayeringStrategy;
import org.eclipse.elk.alg.test.framework.LayoutTestRunner;
import org.eclipse.elk.alg.test.framework.annotations.Algorithm;
import org.eclipse.elk.alg.test.framework.annotations.Configurator;
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
public class BasicLayerAssignmentTest {

    private static final int INTERACTIVE_RANDOM_SEED = 0;
    private static final int INTERACTIVE_MAX_POS = 10_000;
    
    
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
    public LayoutConfigurator coffmanGrahamConfigurator() {
        return configuratorFor(LayeringStrategy.COFFMAN_GRAHAM);
    }
    
    @ConfiguratorProvider
    public LayoutConfigurator longestPathConfigurator() {
        return configuratorFor(LayeringStrategy.LONGEST_PATH);
    }
    
    @ConfiguratorProvider
    public LayoutConfigurator longestPathSourceConfigurator() {
        return configuratorFor(LayeringStrategy.LONGEST_PATH_SOURCE);
    }
    
    @ConfiguratorProvider
    public LayoutConfigurator minWidthConfigurator() {
        return configuratorFor(LayeringStrategy.MIN_WIDTH);
    }
    
    @ConfiguratorProvider
    public LayoutConfigurator networkSimplexConfigurator() {
        return configuratorFor(LayeringStrategy.NETWORK_SIMPLEX);
    }
    
    @ConfiguratorProvider
    public LayoutConfigurator stretchWidthConfigurator() {
        return configuratorFor(LayeringStrategy.STRETCH_WIDTH);
    }
    
    private LayoutConfigurator configuratorFor(final LayeringStrategy strategy) {
        LayoutConfigurator config = new LayoutConfigurator();
        config.configure(ElkNode.class).setProperty(LayeredOptions.LAYERING_STRATEGY, strategy);
        return config;
    }
    
    @Configurator
    public void interactiveConfigurator(final ElkNode elkGraph) {
        // This is a bit more involved. We not only have to set the proper cycle breaking strategy, but we also
        // have to assign bogus coordinates to the nodes
        Deque<ElkNode> nodeQueue = new LinkedList<>();
        nodeQueue.add(elkGraph);
        
        // Fixed seed to keep tests repeatable
        Random rand = new Random(INTERACTIVE_RANDOM_SEED);
        
        while (!nodeQueue.isEmpty()) {
            ElkNode node = nodeQueue.poll();
            nodeQueue.addAll(node.getChildren());
            
            node.setProperty(LayeredOptions.LAYERING_STRATEGY, LayeringStrategy.INTERACTIVE);
            
            if (node.getX() == 0) {
                node.setX(rand.nextInt(INTERACTIVE_MAX_POS));
            }
    
            if (node.getY() == 0) {
                node.setY(rand.nextInt(INTERACTIVE_MAX_POS));
            }
        }
    }
    

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Tests
    
    @TestAfterProcessor(CoffmanGrahamLayerer.class)
    @TestAfterProcessor(InteractiveLayerer.class)
    @TestAfterProcessor(LongestPathLayerer.class)
    @TestAfterProcessor(MinWidthLayerer.class)
    @TestAfterProcessor(NetworkSimplexLayerer.class)
    @TestAfterProcessor(StretchWidthLayerer.class)
    public void testNoLayerlessNodes(final Object graph) {
        LGraph lGraph = (LGraph) graph;
        assertTrue(lGraph.getLayerlessNodes().isEmpty());
    }
    
    @TestAfterProcessor(CoffmanGrahamLayerer.class)
    @TestAfterProcessor(InteractiveLayerer.class)
    @TestAfterProcessor(LongestPathLayerer.class)
    @TestAfterProcessor(MinWidthLayerer.class)
    @TestAfterProcessor(NetworkSimplexLayerer.class)
    @TestAfterProcessor(StretchWidthLayerer.class)
    public void testNoEmptyLayers(final Object graph) {
        LGraph lGraph = (LGraph) graph;
        
        for (Layer layer : lGraph) {
            assertFalse(layer.getNodes().isEmpty());
        }
    }
    
    @TestAfterProcessor(CoffmanGrahamLayerer.class)
    @TestAfterProcessor(InteractiveLayerer.class)
    @TestAfterProcessor(LongestPathLayerer.class)
    @TestAfterProcessor(MinWidthLayerer.class)
    @TestAfterProcessor(NetworkSimplexLayerer.class)
    @TestAfterProcessor(StretchWidthLayerer.class)
    public void testEdgesPointTowardsNextLayers(final Object graph) {
        LGraph lGraph = (LGraph) graph;
        
        // Assign increasing IDs to the layers
        int nextLayerId = 0;
        for (Layer layer : lGraph) {
            layer.id = nextLayerId++;
        }
        
        for (Layer layer : lGraph) {
            for (LNode node : layer) {
                int nodeLayerId = node.getLayer().id;
                
                for (LEdge edge : node.getOutgoingEdges()) {
                    assertTrue(nodeLayerId < edge.getTarget().getNode().getLayer().id);
                }
            }
        }
    }
    
}
