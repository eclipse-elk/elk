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
package org.eclipse.elk.alg.common.issues;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;

import org.eclipse.elk.alg.test.framework.LayoutTestRunner;
import org.eclipse.elk.alg.test.framework.annotations.Algorithm;
import org.eclipse.elk.alg.test.framework.annotations.Configurator;
import org.eclipse.elk.alg.test.framework.annotations.DefaultConfiguration;
import org.eclipse.elk.alg.test.framework.annotations.GraphResourceProvider;
import org.eclipse.elk.alg.test.framework.io.AbstractResourcePath;
import org.eclipse.elk.alg.test.framework.io.ModelResourcePath;
import org.eclipse.elk.core.options.FixedLayouterOptions;
import org.eclipse.elk.graph.ElkNode;
import org.junit.Test;
import org.junit.runner.RunWith;

import com.google.common.collect.Lists;

/**
 * Test for issue 475.
 */
@RunWith(LayoutTestRunner.class)
@Algorithm(FixedLayouterOptions.ALGORITHM_ID)
@DefaultConfiguration
public class Issue475Test {

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Sources

    @GraphResourceProvider
    public List<AbstractResourcePath> testGraphs() {
        return Lists.newArrayList(new ModelResourcePath("tickets/core/475_fixedLayoutWithLongLabels.elkt"));
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Configuration
    
    /** Default compound node size. */
    private static final double SIZE = 20;

    @Configurator
    public void configureFixedGraphSize(final ElkNode graph) {
        for (ElkNode child : graph.getChildren()) {
            // Make sure we know the size the node shall have
            child.setProperty(FixedLayouterOptions.NODE_SIZE_FIXED_GRAPH_SIZE, true);
            child.setHeight(SIZE);
            child.setWidth(SIZE);
        }
    }

    @Configurator
    public void configureNoFixedGraphSize(final ElkNode graph) {
        graph.getChildren().stream()
            .forEach(node -> node.setProperty(FixedLayouterOptions.NODE_SIZE_FIXED_GRAPH_SIZE, false));
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Tests
    
    @Test
    public void testCompoundNodeSize(final ElkNode graph) {
        graph.getChildren().stream().forEach(node -> doTestCompoundNodeSize(node));
    }

    public void doTestCompoundNodeSize(final ElkNode node) {
        if (node.getProperty(FixedLayouterOptions.NODE_SIZE_FIXED_GRAPH_SIZE)) {
            // The node's size should not have changed
            assertEquals(SIZE, node.getWidth(), 0);
            assertEquals(SIZE, node.getHeight(), 0);
            
        } else {
            double maxX = 0;
            double maxY = 0;
            
            for (ElkNode child : node.getChildren()) {
                maxX = Math.max(maxX, child.getX() + child.getWidth());
                maxY = Math.max(maxY, child.getY() + child.getHeight());
            }
            
            assertTrue(node.getWidth() >= maxX);
            assertTrue(node.getHeight() >= maxY);
        }
    }

}
