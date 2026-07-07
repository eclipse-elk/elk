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


/**
 * Tests node margin calculation, which is spread out over several processors.
 */
@RunWith(LayoutTestRunner.class)
@Algorithm(LayeredOptions.ALGORITHM_ID)
@DefaultConfiguration()
public class NodeMarginCalculatorTest {

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Sources

    @GraphResourceProvider
    public List<AbstractResourcePath> testGraphs() {
        return Lists.newArrayList(
                new ModelResourcePath("realworld/ptolemy/**/").withFilter(new FileExtensionFilter("elkg")));
    }


    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Tests
    
    // There is one method for each intermediate processor that computes node margins. The idea is to make it
    // clear which processor failed the test. Another way to solve this would be to extend the test framework
    // to list a method triggered after several processors once for each such processor.
    
    @TestAfterProcessor(InnermostNodeMarginCalculator.class)
    public void testNodeMarginsAfterInnermostNodeMarginCalculator(final Object graph) {
        testNodeMargins(graph);
    }
    
    @TestAfterProcessor(SelfLoopRouter.class)
    public void testNodeMarginsAfterSelfLoopRouter(final Object graph) {
        testNodeMargins(graph);
    }
    
    @TestAfterProcessor(CommentNodeMarginCalculator.class)
    public void testNodeMarginsAfterCommentNodeMarginCalculator(final Object graph) {
        testNodeMargins(graph);
    }
    
    @TestAfterProcessor(EndLabelPreprocessor.class)
    public void testNodeMarginsAfterEndLabelPreprocessor(final Object graph) {
        testNodeMargins(graph);
    }
    
    private void testNodeMargins(final Object graph) {
        LGraph lGraph = (LGraph) graph;
        
        for (Layer layer : lGraph) {
            for (LNode node : layer) {
                assertTrue(node.getMargin() != null);
                assertTrue(node.getMargin().top >= 0);
                assertTrue(node.getMargin().right >= 0);
                assertTrue(node.getMargin().bottom >= 0);
                assertTrue(node.getMargin().left >= 0);
            }
        }
    }

}
