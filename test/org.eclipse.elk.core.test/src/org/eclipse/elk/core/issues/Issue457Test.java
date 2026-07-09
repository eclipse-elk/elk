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
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.options.BoxLayouterOptions;
import org.eclipse.elk.core.options.SizeConstraint;
import org.eclipse.elk.graph.ElkNode;
import org.junit.Test;
import org.junit.runner.RunWith;

import com.google.common.collect.Lists;

/**
 * Test for issue 457.
 */
@RunWith(LayoutTestRunner.class)
@Algorithm(BoxLayouterOptions.ALGORITHM_ID)
@DefaultConfiguration
public class Issue457Test {

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Sources

    @GraphResourceProvider
    public List<AbstractResourcePath> testGraphs() {
        return Lists.newArrayList(new ModelResourcePath("tickets/core/457_boxLayouterIgnoresMinimumSize.elkt"));
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Configuration
    
    private static final double MIN_SIZE = 60;

    @Configurator
    public void configureMinSize(final ElkNode graph) {
        // Set the minimum size of the nodes on the top level
        for (ElkNode node : graph.getChildren()) {
            node.setProperty(BoxLayouterOptions.NODE_SIZE_CONSTRAINTS, SizeConstraint.minimumSize());
            node.setProperty(BoxLayouterOptions.NODE_SIZE_MINIMUM, new KVector(MIN_SIZE, MIN_SIZE));
            
            // Make sure the nodes in there are small
            for (ElkNode child : node.getChildren()) {
                child.setDimensions(5, 5);
            }
        }
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Tests

    @Test
    public void testProperLabelCoordinates(final ElkNode graph) {
        // Ensure that the nodes have the proper size
        for (ElkNode node : graph.getChildren()) {
            assertEquals(MIN_SIZE, node.getWidth(), 1);
            assertEquals(MIN_SIZE, node.getHeight(), 1);
        }
    }

}
