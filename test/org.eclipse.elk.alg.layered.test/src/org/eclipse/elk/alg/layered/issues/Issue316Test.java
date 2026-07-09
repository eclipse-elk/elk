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

import static org.junit.Assert.fail;

import java.util.List;

import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.test.GraphTestUtils;
import org.eclipse.elk.alg.test.framework.LayoutTestRunner;
import org.eclipse.elk.alg.test.framework.annotations.Algorithm;
import org.eclipse.elk.alg.test.framework.annotations.DefaultConfiguration;
import org.eclipse.elk.alg.test.framework.annotations.GraphResourceProvider;
import org.eclipse.elk.alg.test.framework.io.AbstractResourcePath;
import org.eclipse.elk.alg.test.framework.io.ModelResourcePath;
import org.eclipse.elk.graph.ElkLabel;
import org.eclipse.elk.graph.ElkNode;
import org.junit.Test;
import org.junit.runner.RunWith;

import com.google.common.collect.Lists;

/**
 * Tests whether node size calculation is properly applied to hierarchical nodes if hierarchy handling is set to
 * {@code INCLUDE_CHILDREN}.
 */
@RunWith(LayoutTestRunner.class)
@Algorithm(LayeredOptions.ALGORITHM_ID)
@DefaultConfiguration()
public class Issue316Test {

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Sources

    @GraphResourceProvider
    public List<AbstractResourcePath> testGraphs() {
        return Lists.newArrayList(
                new ModelResourcePath("tickets/layered/316_wrongGraphSizeWithIncludeChildren.elkt"),
                new ModelResourcePath("tickets/layered/316_wrongGraphSizeWithIncludeChildren2.elkt"));
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Tests

    @Test
    public void ensurePortsDontOverlap(final ElkNode graph) {
        boolean overlaps = GraphTestUtils.allNodes(graph, false).stream()
            .map(node -> node.getPorts())
            .anyMatch(GraphTestUtils::haveOverlaps);
        
        if (overlaps) {
            fail("Port overlaps detected!");
        }
    }

    @Test
    public void ensureNodesBigEnoughForLabels(final ElkNode graph) {
        for (ElkNode node : GraphTestUtils.allNodes(graph, false)) {
            // We only care about hierarchical nodes
            if (!node.getChildren().isEmpty()) {
                for (ElkLabel label : node.getLabels()) {
                    if (label.getX() < 0
                            || label.getY() < 0
                            || label.getX() + label.getWidth() > node.getWidth()
                            || label.getY() + label.getHeight() > node.getHeight()) {
                        
                        fail("Node label not fully contained in node.");
                    }
                }
            }
        }
    }

}
