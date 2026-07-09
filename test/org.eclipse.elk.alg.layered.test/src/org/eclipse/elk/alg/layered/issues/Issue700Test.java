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

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;

import org.eclipse.elk.alg.test.GraphTestUtils;
import org.eclipse.elk.alg.test.framework.LayoutTestRunner;
import org.eclipse.elk.alg.test.framework.annotations.DefaultConfiguration;
import org.eclipse.elk.alg.test.framework.annotations.GraphResourceProvider;
import org.eclipse.elk.alg.test.framework.io.AbstractResourcePath;
import org.eclipse.elk.alg.test.framework.io.FileNameFilter;
import org.eclipse.elk.alg.test.framework.io.ModelResourcePath;
import org.eclipse.elk.graph.ElkEdge;
import org.eclipse.elk.graph.ElkEdgeSection;
import org.eclipse.elk.graph.ElkNode;
import org.junit.Test;
import org.junit.runner.RunWith;

import com.google.common.collect.Lists;

/**
 * Test for issue 700.
 */
@RunWith(LayoutTestRunner.class)
@DefaultConfiguration()
public class Issue700Test {

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Sources

    @GraphResourceProvider
    public List<AbstractResourcePath> testGraphs() {
        // The bug was introduced while fixing #665, so ensure that the fix for #700 does not return the favour...
        return Lists.newArrayList(
                new ModelResourcePath("tickets/layered/**").withFilter(new FileNameFilter("665.+\\.elkt")),
                new ModelResourcePath("tickets/layered/**").withFilter(new FileNameFilter("700.+\\.elkt")),
                new ModelResourcePath("examples/hierarchy/**"));
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Tests
    
    @Test
    public void testEdgesRouted(final ElkNode graph) {
        GraphTestUtils.allNodes(graph, true).stream()
            .flatMap(node -> node.getContainedEdges().stream())
            .forEach(this::testEdgeRouted);
    }
    
    private void testEdgeRouted(final ElkEdge edge) {
        assertEquals("Edge has no sections.", 1, edge.getSections().size());
        
        // The start and end point may not be at the origin
        ElkEdgeSection section = edge.getSections().get(0);
        assertTrue("Edge has no proper start and end coordinates.",
                section.getStartX() != 0
                || section.getStartY() != 0
                || section.getEndX() != 0
                || section.getEndY() != 0);
    }

}
