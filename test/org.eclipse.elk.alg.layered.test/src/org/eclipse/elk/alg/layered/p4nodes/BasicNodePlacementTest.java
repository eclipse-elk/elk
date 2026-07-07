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
package org.eclipse.elk.alg.layered.p4nodes;

import static org.junit.Assert.assertTrue;

import java.util.Iterator;
import java.util.List;

import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.intermediate.greedyswitch.TestGraphCreator;
import org.eclipse.elk.alg.layered.options.EdgeStraighteningStrategy;
import org.eclipse.elk.alg.layered.options.FixedAlignment;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.layered.options.NodePlacementStrategy;
import org.eclipse.elk.alg.layered.p4nodes.bk.BKNodePlacer;
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
public class BasicNodePlacementTest extends TestGraphCreator {

    private static final int INTERACTIVE_RANDOM_SEED = 0;
    private static final int INTERACTIVE_MAX_POS = 10_000;
    
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Sources

    @GraphResourceProvider
    public List<AbstractResourcePath> testGraphs() {
        return Lists.newArrayList(
                new ModelResourcePath("realworld/ptolemy/**/").withFilter(new FileExtensionFilter("elkg")),
                new ModelResourcePath("tests/layered/node_placement/**/"));
    }
    

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Configuration
    
    // This would be a lot nicer with parameterized tests, but JUnit 4 requires those to be run with the Parameterized
    // test runner. We already use a custom runner, so that's out of the question. :(
    
    @ConfiguratorProvider
    public LayoutConfigurator linearSegmentsConfigurator() {
        return configuratorFor(NodePlacementStrategy.LINEAR_SEGMENTS, null, null);
    }
    
    @ConfiguratorProvider
    public LayoutConfigurator networkSimplexConfigurator() {
        return configuratorFor(NodePlacementStrategy.NETWORK_SIMPLEX, null, null);
    }
    
    @ConfiguratorProvider
    public LayoutConfigurator simpleConfigurator() {
        return configuratorFor(NodePlacementStrategy.SIMPLE, null, null);
    }
    
    @ConfiguratorProvider
    public LayoutConfigurator bkNoneNoneConfigurator() {
        return configuratorFor(
                NodePlacementStrategy.BRANDES_KOEPF,
                EdgeStraighteningStrategy.NONE,
                FixedAlignment.NONE);
    }
    
    @ConfiguratorProvider
    public LayoutConfigurator bkNoneLeftUpConfigurator() {
        return configuratorFor(
                NodePlacementStrategy.BRANDES_KOEPF,
                EdgeStraighteningStrategy.NONE,
                FixedAlignment.LEFTUP);
    }
    
    @ConfiguratorProvider
    public LayoutConfigurator bkNoneLeftDownConfigurator() {
        return configuratorFor(
                NodePlacementStrategy.BRANDES_KOEPF,
                EdgeStraighteningStrategy.NONE,
                FixedAlignment.LEFTDOWN);
    }
    
    @ConfiguratorProvider
    public LayoutConfigurator bkNoneRightUpConfigurator() {
        return configuratorFor(
                NodePlacementStrategy.BRANDES_KOEPF,
                EdgeStraighteningStrategy.NONE,
                FixedAlignment.RIGHTUP);
    }
    
    @ConfiguratorProvider
    public LayoutConfigurator bkNoneRightDownConfigurator() {
        return configuratorFor(
                NodePlacementStrategy.BRANDES_KOEPF,
                EdgeStraighteningStrategy.NONE,
                FixedAlignment.RIGHTDOWN);
    }
    
    @ConfiguratorProvider
    public LayoutConfigurator bkNoneBalancedConfigurator() {
        return configuratorFor(
                NodePlacementStrategy.BRANDES_KOEPF,
                EdgeStraighteningStrategy.NONE,
                FixedAlignment.BALANCED);
    }
    
    @ConfiguratorProvider
    public LayoutConfigurator bkImproveStraightnessLeftUpConfigurator() {
        return configuratorFor(
                NodePlacementStrategy.BRANDES_KOEPF,
                EdgeStraighteningStrategy.IMPROVE_STRAIGHTNESS,
                FixedAlignment.LEFTUP);
    }
    
    @ConfiguratorProvider
    public LayoutConfigurator bkImproveStraightnessLeftDownConfigurator() {
        return configuratorFor(
                NodePlacementStrategy.BRANDES_KOEPF,
                EdgeStraighteningStrategy.IMPROVE_STRAIGHTNESS,
                FixedAlignment.LEFTDOWN);
    }
    
    @ConfiguratorProvider
    public LayoutConfigurator bkImproveStraightnessRightUpConfigurator() {
        return configuratorFor(
                NodePlacementStrategy.BRANDES_KOEPF,
                EdgeStraighteningStrategy.IMPROVE_STRAIGHTNESS,
                FixedAlignment.RIGHTUP);
    }
    
    @ConfiguratorProvider
    public LayoutConfigurator bkImproveStraightnessRightDownConfigurator() {
        return configuratorFor(
                NodePlacementStrategy.BRANDES_KOEPF,
                EdgeStraighteningStrategy.IMPROVE_STRAIGHTNESS,
                FixedAlignment.RIGHTDOWN);
    }
    
    @ConfiguratorProvider
    public LayoutConfigurator bkImproveStraightnessBalancedConfigurator() {
        return configuratorFor(
                NodePlacementStrategy.BRANDES_KOEPF,
                EdgeStraighteningStrategy.IMPROVE_STRAIGHTNESS,
                FixedAlignment.BALANCED);
    }
    
    private LayoutConfigurator configuratorFor(final NodePlacementStrategy strategy,
            final EdgeStraighteningStrategy bkStraighteningStrategy, final FixedAlignment bkFixedAlignment) {
        
        LayoutConfigurator config = new LayoutConfigurator();
        config.configure(ElkNode.class)
            .setProperty(LayeredOptions.NODE_PLACEMENT_STRATEGY, strategy)
            .setProperty(LayeredOptions.NODE_PLACEMENT_BK_FIXED_ALIGNMENT, bkFixedAlignment)
            .setProperty(LayeredOptions.NODE_PLACEMENT_BK_EDGE_STRAIGHTENING, bkStraighteningStrategy);
        return config;
    }
    

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Tests
    
    @TestAfterProcessor(BKNodePlacer.class)
    @TestAfterProcessor(LinearSegmentsNodePlacer.class)
    @TestAfterProcessor(NetworkSimplexPlacer.class)
    @TestAfterProcessor(SimpleNodePlacer.class)
    public void testProperCoordinates(final Object graph) {
        LGraph lGraph = (LGraph) graph;
        
        for (Layer layer : lGraph) {
            Iterator<LNode> nodeIter = layer.getNodes().iterator();
            LNode node = nodeIter.next();
            double lastBottomY = node.getPosition().y + node.getSize().y;
            
            while (nodeIter.hasNext()) {
                node = nodeIter.next();
                assertTrue(lastBottomY < node.getPosition().y);
                
                lastBottomY = node.getPosition().y + node.getSize().y;
            }
        }
    }
    
}
