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

import static org.junit.Assert.assertTrue;

import java.util.List;

import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.test.framework.LayoutTestRunner;
import org.eclipse.elk.alg.test.framework.annotations.Algorithm;
import org.eclipse.elk.alg.test.framework.annotations.Configurator;
import org.eclipse.elk.alg.test.framework.annotations.DefaultConfiguration;
import org.eclipse.elk.alg.test.framework.annotations.GraphResourceProvider;
import org.eclipse.elk.alg.test.framework.io.AbstractResourcePath;
import org.eclipse.elk.alg.test.framework.io.ModelResourcePath;
import org.eclipse.elk.core.math.ElkPadding;
import org.eclipse.elk.graph.ElkNode;
import org.junit.Test;
import org.junit.runner.RunWith;

import com.google.common.collect.Lists;

/**
 * Test for issue 445.
 */
@RunWith(LayoutTestRunner.class)
@Algorithm(LayeredOptions.ALGORITHM_ID)
@DefaultConfiguration
public class Issue445Test {

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Sources

    @GraphResourceProvider
    public List<AbstractResourcePath> testGraphs() {
        return Lists.newArrayList(new ModelResourcePath("tickets/layered/445_feedbackEdgeTailLabel.elkt"));
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Configuration

    @Configurator
    public void configurePadding(final ElkNode graph) {
        // Set the graph's padding to zero. Due to this bug, the label will end up at (0,0) as well, then
        graph.setProperty(LayeredOptions.PADDING, new ElkPadding(0));
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Tests

    /**
     * Check that the property was set on the graph.
     */
    @Test
    public void testProperLabelCoordinates(final ElkNode graph) {
        long misplacedLabels = graph.getContainedEdges().stream()
                .flatMap(edge -> edge.getLabels().stream())
                .filter(label -> label.getX() == 0 && label.getY() == 0)
                .count();
        assertTrue("Found " + misplacedLabels + " label(s) placed at (0, 0) instead of a proper position.",
                misplacedLabels == 0);
    }

}
