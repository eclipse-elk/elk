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

import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.test.framework.LayoutTestRunner;
import org.eclipse.elk.alg.test.framework.annotations.Algorithm;
import org.eclipse.elk.alg.test.framework.annotations.DefaultConfiguration;
import org.eclipse.elk.alg.test.framework.annotations.GraphResourceProvider;
import org.eclipse.elk.alg.test.framework.annotations.TestAfterProcessor;
import org.eclipse.elk.alg.test.framework.io.AbstractResourcePath;
import org.eclipse.elk.alg.test.framework.io.FileExtensionFilter;
import org.eclipse.elk.alg.test.framework.io.ModelResourcePath;
import org.junit.Assert;
import org.junit.runner.RunWith;

import com.google.common.collect.Lists;

@RunWith(LayoutTestRunner.class)
@Algorithm(LayeredOptions.ALGORITHM_ID)
@DefaultConfiguration()
public class LongEdgeSplitterTest {

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Sources

    @GraphResourceProvider
    public List<AbstractResourcePath> testGraphs() {
        return Lists.newArrayList(
                new ModelResourcePath("realworld/ptolemy/**/").withFilter(new FileExtensionFilter("elkg")));
    }


    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Tests

    /**
     * All edges connect directly adjacent layer where the source layer's index is smaller than the target layer's
     * index.
     */
    @TestAfterProcessor(LongEdgeSplitter.class)
    public void testEdgesBetweenAdjacentLayers(final Object graph) {
        LGraph lGraph = (LGraph) graph;
        
        // Assign layer IDs
        int nextLayerId = 0;
        for (Layer layer : lGraph) {
            layer.id = nextLayerId++;
        }
        
        for (Layer layer : lGraph) {
            for (LNode node : layer) {
                for (LEdge edge : node.getOutgoingEdges()) {
                    // increasing layer ID
                    Assert.assertTrue(layer.id < edge.getTarget().getNode().getLayer().id);
                    
                    // directly adjacent index
                    Assert.assertTrue(layer.id + 1 == edge.getTarget().getNode().getLayer().id);
                }
            }
        }
    }

}
