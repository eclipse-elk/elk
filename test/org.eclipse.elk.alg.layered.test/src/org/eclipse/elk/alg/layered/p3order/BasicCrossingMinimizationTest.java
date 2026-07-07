/*******************************************************************************
 * Copyright (c) 2013, 2020 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.layered.p3order;

import static org.junit.Assert.assertTrue;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.intermediate.greedyswitch.TestGraphCreator;
import org.eclipse.elk.alg.layered.options.CrossingMinimizationStrategy;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.test.framework.LayoutTestRunner;
import org.eclipse.elk.alg.test.framework.annotations.Algorithm;
import org.eclipse.elk.alg.test.framework.annotations.ConfiguratorProvider;
import org.eclipse.elk.alg.test.framework.annotations.DefaultConfiguration;
import org.eclipse.elk.alg.test.framework.annotations.GraphResourceProvider;
import org.eclipse.elk.alg.test.framework.annotations.TestAfterProcessor;
import org.eclipse.elk.alg.test.framework.annotations.TestBeforeProcessor;
import org.eclipse.elk.alg.test.framework.io.AbstractResourcePath;
import org.eclipse.elk.alg.test.framework.io.ModelResourcePath;
import org.eclipse.elk.core.LayoutConfigurator;
import org.eclipse.elk.graph.ElkNode;
import org.junit.runner.RunWith;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;

@RunWith(LayoutTestRunner.class)
@Algorithm(LayeredOptions.ALGORITHM_ID)
@DefaultConfiguration()
public class BasicCrossingMinimizationTest extends TestGraphCreator {
    
    private Map<LGraph, Multimap<Layer, LNode>> layerAssignments = new HashMap<>();

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Sources

    @GraphResourceProvider
    public List<AbstractResourcePath> testGraphs() {
        return Lists.newArrayList(new ModelResourcePath("tests/**/"));
    }
    

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Configuration

    @ConfiguratorProvider
    public LayoutConfigurator layerSweepConfigurator() {
        LayoutConfigurator config = new LayoutConfigurator();
        config.configure(ElkNode.class).setProperty(
                LayeredOptions.CROSSING_MINIMIZATION_STRATEGY,
                CrossingMinimizationStrategy.LAYER_SWEEP);
        return config;
    }
    

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Tests

    @TestBeforeProcessor(LayerSweepCrossingMinimizer.class)
    @TestAfterProcessor(LayerSweepCrossingMinimizer.class)
    public void ensureLayerAssignmentUnchanged(final Object graph) {
        LGraph lgraph = (LGraph) graph;
        
        if (layerAssignments.containsKey(lgraph)) {
            verifyLayerAssignment(lgraph);
        } else {
            recordLayerAssignment(lgraph);
        }
    }

    private void verifyLayerAssignment(LGraph lgraph) {
        Multimap<Layer, LNode> layerAssignment = layerAssignments.get(lgraph);
        layerAssignments.remove(lgraph);
        
        for (Layer layer : lgraph) {
            for (LNode node : layer) {
                // node in the same layer?
                assertTrue(layerAssignment.get(layer).contains(node));
            }
        }
    }

    private void recordLayerAssignment(LGraph lgraph) {
        Multimap<Layer, LNode> layerAssignment = HashMultimap.create();
        
        for (Layer layer : lgraph) {
            for (LNode node : layer) {
                layerAssignment.put(layer, node);
            }
        }
        
        layerAssignments.put(lgraph, layerAssignment);
    }

}
