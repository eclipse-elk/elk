/*******************************************************************************
 * Copyright (c) 2020 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.test.tests;

import static org.junit.Assert.assertTrue;

import org.eclipse.elk.alg.disco.options.DisCoOptions;
import org.eclipse.elk.alg.force.options.ForceOptions;
import org.eclipse.elk.alg.force.options.StressOptions;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.mrtree.options.MrTreeOptions;
import org.eclipse.elk.alg.radial.options.RadialOptions;
import org.eclipse.elk.alg.rectpacking.options.RectPackingOptions;
import org.eclipse.elk.alg.spore.options.SporeCompactionOptions;
import org.eclipse.elk.alg.spore.options.SporeOverlapRemovalOptions;
import org.eclipse.elk.alg.test.framework.LayoutTestRunner;
import org.eclipse.elk.alg.test.framework.annotations.Algorithm;
import org.eclipse.elk.alg.test.framework.annotations.GraphProvider;
import org.eclipse.elk.core.options.BoxLayouterOptions;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.options.FixedLayouterOptions;
import org.eclipse.elk.core.options.RandomLayouterOptions;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.util.ElkGraphUtil;
import org.junit.Test;
import org.junit.runner.RunWith;

/**
 * Test that checks whether all expected algorithms are registered properly. If an algorithm is not registered properly,
 * one of the tests here will fail.
 */
@RunWith(LayoutTestRunner.class)
@Algorithm(BoxLayouterOptions.ALGORITHM_ID)
@Algorithm(DisCoOptions.ALGORITHM_ID)
@Algorithm(FixedLayouterOptions.ALGORITHM_ID)
@Algorithm(ForceOptions.ALGORITHM_ID)
@Algorithm(LayeredOptions.ALGORITHM_ID)
@Algorithm(MrTreeOptions.ALGORITHM_ID)
@Algorithm(RadialOptions.ALGORITHM_ID)
@Algorithm(RandomLayouterOptions.ALGORITHM_ID)
@Algorithm(RectPackingOptions.ALGORITHM_ID)
@Algorithm(SporeCompactionOptions.ALGORITHM_ID)
@Algorithm(SporeOverlapRemovalOptions.ALGORITHM_ID)
@Algorithm(StressOptions.ALGORITHM_ID)
public class PlainJavaInitializationTest {

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Sources

    /**
     * Supply a basic graph.
     */
    @GraphProvider
    public ElkNode basicGraph() {
        return ElkGraphUtil.createGraph();
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Tests

    @Test
    public void testAlgorithmPresent(final ElkNode graph) {
        // A simple tummy test
        assertTrue(graph.hasProperty(CoreOptions.ALGORITHM));
    }

}
