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

import java.util.Deque;
import java.util.LinkedList;
import java.util.List;
import java.util.Random;

import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.intermediate.greedyswitch.TestGraphCreator;
import org.eclipse.elk.alg.layered.options.InLayerConstraint;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.test.framework.LayoutTestRunner;
import org.eclipse.elk.alg.test.framework.annotations.Algorithm;
import org.eclipse.elk.alg.test.framework.annotations.Configurator;
import org.eclipse.elk.alg.test.framework.annotations.DefaultConfiguration;
import org.eclipse.elk.alg.test.framework.annotations.GraphResourceProvider;
import org.eclipse.elk.alg.test.framework.annotations.TestAfterProcessor;
import org.eclipse.elk.alg.test.framework.io.AbstractResourcePath;
import org.eclipse.elk.alg.test.framework.io.FileExtensionFilter;
import org.eclipse.elk.alg.test.framework.io.ModelResourcePath;
import org.eclipse.elk.graph.ElkNode;
import org.junit.runner.RunWith;

import com.google.common.collect.Lists;

@RunWith(LayoutTestRunner.class)
@Algorithm(LayeredOptions.ALGORITHM_ID)
@DefaultConfiguration()
public class InLayerConstraintProcessorTest extends TestGraphCreator {

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Sources

    @GraphResourceProvider
    public List<AbstractResourcePath> testGraphs() {
        return Lists.newArrayList(
                new ModelResourcePath("realworld/ptolemy/**/").withFilter(new FileExtensionFilter("elkg")));
    }
    

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Configuration
    
    /** Fixed seed to keep tests reproducible. */
    private static final int RANDOM_SEED = 0;
    
    @Configurator
    public void configureRandomInLayerConstraints(final ElkNode graph) {
        Deque<ElkNode> nodeQueue = new LinkedList<>(graph.getChildren());
        Random rand = new Random(RANDOM_SEED);
        
        while (!nodeQueue.isEmpty()) {
            ElkNode node = nodeQueue.poll();
            nodeQueue.addAll(node.getChildren());
            
            double r = rand.nextDouble();
            if (r < 0.2) {
                node.setProperty(InternalProperties.IN_LAYER_CONSTRAINT, InLayerConstraint.TOP);
            } else if (r < 0.8) {
                node.setProperty(InternalProperties.IN_LAYER_CONSTRAINT, InLayerConstraint.NONE);
            } else {
                node.setProperty(InternalProperties.IN_LAYER_CONSTRAINT, InLayerConstraint.BOTTOM);
            }
        }
    }
    

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Tests
    
    @TestAfterProcessor(InLayerConstraintProcessor.class)
    public void testValidNodeOrder(final Object graph) {
        LGraph lGraph = (LGraph) graph;
        
        for (Layer layer : lGraph) {
            InLayerConstraint lastConstraint = null;
            for (LNode node : layer) {
                InLayerConstraint currentConstraint = node.getProperty(InternalProperties.IN_LAYER_CONSTRAINT);

                if (lastConstraint != null && currentConstraint != lastConstraint) {
                    // if the value changes check valid transitions
                    String error = "Invalid constraint transition: " + lastConstraint.name()
                        + " -> " + currentConstraint.name();
                    
                    if (currentConstraint == InLayerConstraint.NONE) {
                        assertTrue(error, lastConstraint == InLayerConstraint.TOP);
                        
                    } else if (currentConstraint == InLayerConstraint.BOTTOM) {
                        assertTrue(error, lastConstraint == InLayerConstraint.TOP
                                || lastConstraint == InLayerConstraint.NONE);
                        
                    }

                    lastConstraint = currentConstraint;
                }
                
                lastConstraint = currentConstraint;
            }
        }
    }
    
}
