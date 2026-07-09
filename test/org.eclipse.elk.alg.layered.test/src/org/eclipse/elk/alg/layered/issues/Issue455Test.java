/*******************************************************************************
 * Copyright (c) 2020 Kiel University and others.
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
package org.eclipse.elk.alg.layered.issues;

import static org.junit.Assert.assertTrue;

import java.util.DoubleSummaryStatistics;
import java.util.List;

import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.test.framework.LayoutTestRunner;
import org.eclipse.elk.alg.test.framework.annotations.Algorithm;
import org.eclipse.elk.alg.test.framework.annotations.DefaultConfiguration;
import org.eclipse.elk.alg.test.framework.annotations.GraphResourceProvider;
import org.eclipse.elk.alg.test.framework.io.AbstractResourcePath;
import org.eclipse.elk.alg.test.framework.io.ModelResourcePath;
import org.eclipse.elk.graph.ElkNode;
import org.junit.Test;
import org.junit.runner.RunWith;

import com.google.common.collect.Lists;

/**
 * Test for issue 455.
 */
@RunWith(LayoutTestRunner.class)
@Algorithm(LayeredOptions.ALGORITHM_ID)
@DefaultConfiguration
public class Issue455Test {
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Sources

    @GraphResourceProvider
    public List<AbstractResourcePath> testGraphs() {
        return Lists.newArrayList(
                new ModelResourcePath("tickets/layered/455_longHierarchicalEdges.elkt"),
                new ModelResourcePath("tickets/layered/455_shortHierarchicalEdges.elkt"));
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Tests
    
    /** Allowed delta when comparing coordinates. */
    private final double COORDINATE_FUZZYNESS = 0.5;

    @Test
    public void testChildNodesAreInSameLayer(final ElkNode graph) {
        // We want to inspect the graph's compound node
        graph.getChildren().stream()
            .filter(child -> child.isHierarchical())
            .forEach(node -> doTest(node));
    }
    
    private void doTest(final ElkNode node) {
        DoubleSummaryStatistics stats = node.getChildren().stream()
            .mapToDouble(child -> child.getX())
            .summaryStatistics();
        
        assertTrue("Child nodes are not placed in the same layer.",
                stats.getMax() - stats.getMin() <= COORDINATE_FUZZYNESS);
    }

}
