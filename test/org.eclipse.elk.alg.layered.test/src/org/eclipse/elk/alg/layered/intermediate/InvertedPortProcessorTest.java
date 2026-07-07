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

import static org.junit.Assert.assertTrue;

import java.util.List;

import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LNode.NodeType;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.intermediate.greedyswitch.TestGraphCreator;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
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
import org.eclipse.elk.core.options.EdgeRouting;
import org.eclipse.elk.core.options.PortSide;
import org.eclipse.elk.graph.ElkNode;
import org.junit.runner.RunWith;

import com.google.common.collect.Lists;

@RunWith(LayoutTestRunner.class)
@Algorithm(LayeredOptions.ALGORITHM_ID)
@DefaultConfiguration()
public class InvertedPortProcessorTest extends TestGraphCreator {

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
    public LayoutConfigurator configurator() {
        LayoutConfigurator config = new LayoutConfigurator();
        config.configure(ElkNode.class).setProperty(LayeredOptions.EDGE_ROUTING, EdgeRouting.ORTHOGONAL);
        return config;
    }
    

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Tests
    
    @TestAfterProcessor(InvertedPortProcessor.class)
    public void testValidNodeOrder(final Object graph) {
        LGraph lGraph = (LGraph) graph;
        
        for (Layer layer : lGraph) {
            for (LNode node : layer) {
                // outgoing edges left
                for (LEdge edge : node.getOutgoingEdges()) {
                    if (edge.getSource().getSide() == PortSide.WEST) {
                        LNode target = edge.getTarget().getNode();
                        
                        // dummy in same layer
                        assertTrue(target.getType() == NodeType.LONG_EDGE);
                        assertTrue(target.getLayer().getIndex() == layer.getIndex());
                    }
                }

                // incoming edges right
                for (LEdge edge : node.getIncomingEdges()) {
                    if (edge.getTarget().getSide() == PortSide.EAST) {
                        LNode source = edge.getSource().getNode();
                        
                        // dummy in same layer
                        assertTrue(source.getType() == NodeType.LONG_EDGE);
                        assertTrue(source.getLayer().getIndex() == layer.getIndex());
                    }
                }
            }
        }
    }
    
}
