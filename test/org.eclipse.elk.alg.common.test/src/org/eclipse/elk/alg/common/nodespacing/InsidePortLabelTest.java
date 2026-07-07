/*******************************************************************************
 * Copyright (c) 2019 Kiel University and others.
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
package org.eclipse.elk.alg.common.nodespacing;

import static org.junit.Assert.fail;

import java.util.ArrayList;
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
 * Tests whether node size computation handles inside port labels properly. Since the test graph contains a node
 * without a node label, we need to deactivate the default configuration for nodes.
 */
@RunWith(LayoutTestRunner.class)
@Algorithm(LayeredOptions.ALGORITHM_ID)
@DefaultConfiguration(edges = false, nodes = false, ports = true)
public class InsidePortLabelTest {

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Sources

    @GraphResourceProvider
    public List<AbstractResourcePath> testGraphs() {
        return Lists.newArrayList(new ModelResourcePath("tests/core/node_size/inside_port_labels.elkt"));
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Tests

    /**
     * Check for each node that no two labels overlap.
     */
    @Test
    public void testNoLabelOverlaps(final ElkNode graph) {
        boolean overlaps = graph.getChildren().stream()
            .map(node -> assembleLabels(node))
            .anyMatch(GraphTestUtils::haveOverlaps);
        
        if (overlaps) {
            fail("Overlaps between labels detected!");
        }
    }
    
    private List<ElkLabel> assembleLabels(final ElkNode node) {
        List<ElkLabel> labels = new ArrayList<>();
        
        // Add all node labels
        labels.addAll(node.getLabels());
        node.getPorts().stream()
            .flatMap(port -> port.getLabels().stream())
            .forEach(label -> labels.add(label));
        
        return labels;
    }

}
