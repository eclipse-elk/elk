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
import org.junit.runner.RunWith;

import com.google.common.collect.Lists;

@RunWith(LayoutTestRunner.class)
@Algorithm(LayeredOptions.ALGORITHM_ID)
@DefaultConfiguration()
public class LayerSizeAndGraphHeigthCalculatorTest {

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Sources

    @GraphResourceProvider
    public List<AbstractResourcePath> testGraphs() {
        return Lists.newArrayList(
                new ModelResourcePath("realworld/ptolemy/**/").withFilter(new FileExtensionFilter("elkg")));
    }


    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Tests
    
    private static final double TOLERANCE = 0.1;
    
    @TestAfterProcessor(LayerSizeAndGraphHeightCalculator.class)
    public void testAllNodesContained(final Object graph) {
        LGraph lGraph = (LGraph) graph;
        
        double top = 0.0 - lGraph.getOffset().y;
        double bottom = lGraph.getSize().y - lGraph.getOffset().y;
        
        for (Layer layer : lGraph) {
            for (LNode node : layer) {
                double nodeTop = node.getPosition().y - node.getMargin().top;
                assertTrue(nodeTop > top || Math.abs(nodeTop - top) < TOLERANCE);

                double nodeBottom = node.getPosition().y + node.getSize().y + node.getMargin().bottom;
                assertTrue(nodeBottom < bottom || Math.abs(nodeBottom - bottom) < TOLERANCE);
            }
        }
    }

}
