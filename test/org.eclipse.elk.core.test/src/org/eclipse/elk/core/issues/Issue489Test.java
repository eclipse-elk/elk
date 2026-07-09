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
package org.eclipse.elk.core.issues;

import static org.junit.Assert.assertEquals;

import java.util.List;

import org.eclipse.elk.alg.test.framework.LayoutTestRunner;
import org.eclipse.elk.alg.test.framework.annotations.Algorithm;
import org.eclipse.elk.alg.test.framework.annotations.Configurator;
import org.eclipse.elk.alg.test.framework.annotations.DefaultConfiguration;
import org.eclipse.elk.alg.test.framework.annotations.GraphResourceProvider;
import org.eclipse.elk.alg.test.framework.io.AbstractResourcePath;
import org.eclipse.elk.alg.test.framework.io.ModelResourcePath;
import org.eclipse.elk.core.math.ElkPadding;
import org.eclipse.elk.core.options.BoxLayouterOptions;
import org.eclipse.elk.graph.ElkNode;
import org.junit.Test;
import org.junit.runner.RunWith;

import com.google.common.collect.Lists;

/**
 * Test for issue 489.
 */
@RunWith(LayoutTestRunner.class)
@Algorithm(BoxLayouterOptions.ALGORITHM_ID)
@DefaultConfiguration
public class Issue489Test {

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Sources

    @GraphResourceProvider
    public List<AbstractResourcePath> testGraphs() {
        return Lists.newArrayList(
                new ModelResourcePath("tickets/core/489_paddingAndMinSize_Small.elkt"),
                new ModelResourcePath("tickets/core/489_paddingAndMinSize_Large.elkt"));
    }
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Configuration
    
    private static final double PADDING = 5;

    @Configurator
    public void configurePadding(final ElkNode graph) {
        // Set the padding of the nodes on the top level
        for (ElkNode node : graph.getChildren()) {
            node.setProperty(BoxLayouterOptions.PADDING, new ElkPadding(PADDING));
        }
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Tests
    
    private static final double TOLERANCE = 0.1;

    @Test
    public void testProperChildPlacement(final ElkNode graph) {
        // Ensure that the nodes have the proper size
        for (ElkNode node : graph.getChildren()) {
            // We expect only one child to be in the test graph!
            assertEquals(1, node.getChildren().size());
            
            ElkNode child = node.getChildren().get(0);
            
            assertEquals(PADDING, child.getX(), TOLERANCE);
            assertEquals(PADDING, child.getY(), TOLERANCE);
            
            assertEquals(node.getWidth() - 2 * PADDING, child.getWidth(), TOLERANCE);
            assertEquals(node.getHeight() - 2 * PADDING, child.getHeight(), TOLERANCE);
        }
    }

}
