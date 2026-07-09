/*******************************************************************************
 * Copyright (c) 2019 Kiel University and others.
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

import static org.junit.Assert.assertEquals;

import java.util.List;

import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.test.framework.LayoutTestRunner;
import org.eclipse.elk.alg.test.framework.annotations.Algorithm;
import org.eclipse.elk.alg.test.framework.annotations.DefaultConfiguration;
import org.eclipse.elk.alg.test.framework.annotations.GraphResourceProvider;
import org.eclipse.elk.alg.test.framework.io.AbstractResourcePath;
import org.eclipse.elk.alg.test.framework.io.ModelResourcePath;
import org.eclipse.elk.graph.ElkEdge;
import org.eclipse.elk.graph.ElkLabel;
import org.eclipse.elk.graph.ElkNode;
import org.junit.Test;
import org.junit.runner.RunWith;

import com.google.common.collect.Lists;

/**
 * Test for issue 471.
 */
@RunWith(LayoutTestRunner.class)
@Algorithm(LayeredOptions.ALGORITHM_ID)
@DefaultConfiguration
public class Issue471Test {
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Sources

    @GraphResourceProvider
    public List<AbstractResourcePath> testGraphs() {
        return Lists.newArrayList(new ModelResourcePath("tickets/layered/471_verticalInlineLabelAlignment.elkt"));
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Tests
    
    /** Allowed delta when comparing coordinates. */
    private final double COORDINATE_FUZZYNESS = 0.5;

    @Test
    public void testEdgeLabelAlignment(final ElkNode graph) {
        graph.getContainedEdges().stream()
            .forEach(edge -> testLabels(edge));
    }
    
    private void testLabels(final ElkEdge edge) {
        double minX = Double.MAX_VALUE;
        double maxWidth = 0;
        
        for (ElkLabel label : edge.getLabels()) {
            minX = Math.min(minX, label.getX());
            maxWidth = Math.max(maxWidth, label.getWidth());
        }
        
        // Each label must be centered
        for (ElkLabel label : edge.getLabels()) {
            assertEquals(minX + (maxWidth - label.getWidth()) / 2, label.getX(), COORDINATE_FUZZYNESS);
        }
    }

}
