/*******************************************************************************
 * Copyright (c) 2020 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.layered.issues;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;

import org.eclipse.elk.alg.layered.options.GreedySwitchType;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.test.framework.LayoutTestRunner;
import org.eclipse.elk.alg.test.framework.annotations.Algorithm;
import org.eclipse.elk.alg.test.framework.annotations.Configurator;
import org.eclipse.elk.alg.test.framework.annotations.DefaultConfiguration;
import org.eclipse.elk.alg.test.framework.annotations.GraphResourceProvider;
import org.eclipse.elk.alg.test.framework.io.AbstractResourcePath;
import org.eclipse.elk.alg.test.framework.io.FileNameFilter;
import org.eclipse.elk.alg.test.framework.io.ModelResourcePath;
import org.eclipse.elk.graph.ElkNode;
import org.junit.Test;
import org.junit.runner.RunWith;

import com.google.common.collect.Lists;
import com.google.common.math.DoubleMath;

/**
 * Test for issue 628.
 */
@RunWith(LayoutTestRunner.class)
@Algorithm(LayeredOptions.ALGORITHM_ID)
@DefaultConfiguration(nodes = false)
public class Issue628Test {

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Sources

    @GraphResourceProvider
    public List<AbstractResourcePath> testGraphs() {
        return Lists.newArrayList(
                new ModelResourcePath("tickets/layered/**").withFilter(new FileNameFilter("628.+\\.elkt")));
    }

    @Configurator
    public void configureGreedySwitchHierarchical(final ElkNode graph) {
        graph.setProperty(LayeredOptions.CROSSING_MINIMIZATION_GREEDY_SWITCH_HIERARCHICAL_TYPE,
                GreedySwitchType.TWO_SIDED);
    }

    @Configurator
    public void unconfigureGreedySwitchHierarchical(final ElkNode graph) {
        graph.setProperty(LayeredOptions.CROSSING_MINIMIZATION_GREEDY_SWITCH_HIERARCHICAL_TYPE, null);
    }

    @Configurator
    public void configureGreedySwitch(final ElkNode graph) {
        graph.setProperty(LayeredOptions.CROSSING_MINIMIZATION_GREEDY_SWITCH_TYPE, GreedySwitchType.TWO_SIDED);
    }

    @Configurator
    public void unconfigureGreedySwitch(final ElkNode graph) {
        graph.setProperty(LayeredOptions.CROSSING_MINIMIZATION_GREEDY_SWITCH_TYPE, null);
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Tests

    @Test
    public void testNoLabelNodeOverlaps(final ElkNode graph) {
        ElkNode l1 = graph.getChildren().get(0);
        ElkNode l2_1 = l1.getChildren().stream().filter(n -> n.getIdentifier().equals("l2_1")).findFirst().get();
        ElkNode l2_2 = l1.getChildren().stream().filter(n -> n.getIdentifier().equals("l2_2")).findFirst().get();

        assertTrue("Y-coordinates not equal.", DoubleMath.fuzzyEquals(l2_1.getY(), l2_2.getY(), 1e-5));
        assertEquals("Spacing between the nodes not as expected.", 10.0, l2_2.getX() - (l2_1.getX() + l2_1.getWidth()),
                1e-5);
    }

}
