/*******************************************************************************
 * Copyright (c) 2013, 2020 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.layered.p1cycles;

import static org.junit.Assert.fail;

import java.util.Deque;
import java.util.LinkedList;
import java.util.List;
import java.util.Random;
import java.util.stream.Collectors;

import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.intermediate.greedyswitch.TestGraphCreator;
import org.eclipse.elk.alg.layered.options.CycleBreakingStrategy;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.layered.options.OrderingStrategy;
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
public class BasicCycleBreakerTest extends TestGraphCreator {

    private static final int INTERACTIVE_RANDOM_SEED = 0;
    private static final int INTERACTIVE_MAX_POS = 100;
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Sources

    @GraphResourceProvider
    public List<AbstractResourcePath> testGraphs() {
        return Lists.newArrayList(
                new ModelResourcePath("realworld/ptolemy/**/").withFilter(new FileExtensionFilter("elkg")),
                new ModelResourcePath("tickets/layered/600_outgoingEdgeInLastSeparateNode.elkt"));
    }
    

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Configuration
    
    @ConfiguratorProvider
    public LayoutConfigurator depthFirstConfigurator() {
        return configuratorFor(CycleBreakingStrategy.DEPTH_FIRST);
    }
    
    @ConfiguratorProvider
    public LayoutConfigurator greedyConfigurator() {
        return configuratorFor(CycleBreakingStrategy.GREEDY);
    }
    
    @ConfiguratorProvider
    public LayoutConfigurator greedyModelOrderConfigurator() {
        return configuratorFor(CycleBreakingStrategy.GREEDY_MODEL_ORDER);
    }
    
    @ConfiguratorProvider
    public LayoutConfigurator modelOrderPreferEdgesConfigurator() {
        LayoutConfigurator config = configuratorFor(CycleBreakingStrategy.MODEL_ORDER);
        config.configure(ElkNode.class).setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY,
                OrderingStrategy.PREFER_EDGES);
        return config;
    }
    
    @ConfiguratorProvider
    public LayoutConfigurator modelOrderNodesAndEdgesConfigurator() {
        LayoutConfigurator config = configuratorFor(CycleBreakingStrategy.MODEL_ORDER);
        config.configure(ElkNode.class).setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY,
                OrderingStrategy.NODES_AND_EDGES);
        return config;
    }
    
    private LayoutConfigurator configuratorFor(final CycleBreakingStrategy strategy) {
        LayoutConfigurator config = new LayoutConfigurator();
        config.configure(ElkNode.class).setProperty(LayeredOptions.CYCLE_BREAKING_STRATEGY, strategy);
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
            
            node.setProperty(LayeredOptions.CYCLE_BREAKING_STRATEGY, CycleBreakingStrategy.INTERACTIVE);
            
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
    
    private static final int DFS_UNVISITED = 0;
    private int nextDFSNumber = 1;
    
    @TestAfterProcessor(DepthFirstCycleBreaker.class)
    @TestAfterProcessor(GreedyCycleBreaker.class)
    @TestAfterProcessor(InteractiveCycleBreaker.class)
    @TestAfterProcessor(ModelOrderCycleBreaker.class)
    public void testIsAcyclic(final Object graph) {
        LGraph lGraph = (LGraph) graph;
        
        // Find all nodes that don't have incoming edges
        List<LNode> sources = lGraph.getLayerlessNodes().stream()
                .filter(node -> !node.getIncomingEdges().iterator().hasNext())
                .collect(Collectors.toList());
        
        // Reset all node IDs
        lGraph.getLayerlessNodes().stream()
                .forEach(node -> node.id = DFS_UNVISITED);
        
        // Start a DFS at each source
        for (LNode source : sources) {
            if (source.id != DFS_UNVISITED) {
                fail("Should never happen");
            }
            
            dfs(source);
        }
    }

    private void dfs(final LNode node) {
        if (node.id != DFS_UNVISITED) {
            return;
        }
        
        node.id = nextDFSNumber++;
        putOnCurrentStack(node);
        
        for (LEdge lEdge : node.getOutgoingEdges()) {
            LNode target = lEdge.getTarget().getNode();
            
            if (target.id == DFS_UNVISITED) {
                dfs(target);
            } else if (isOnCurrentStack(target)) {
                fail("Cycle detected!");
            }
        }
        
        removeFromCurrentStack(node);
    }
    
    private boolean isOnCurrentStack(final LNode node) {
        return node.id < 0;
    }
    
    private void putOnCurrentStack(final LNode node) {
        node.id = -1 * Math.abs(node.id);
    }
    
    private void removeFromCurrentStack(final LNode node) {
        node.id = Math.abs(node.id);
    }
    
}
