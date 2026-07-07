/*******************************************************************************
 * Copyright (c) 2018, 2019 Kiel University and others.
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
package org.eclipse.elk.alg.test.tests;

import static org.junit.Assert.assertTrue;

import org.eclipse.elk.alg.force.options.ForceOptions;
import org.eclipse.elk.alg.test.framework.LayoutTestRunner;
import org.eclipse.elk.alg.test.framework.annotations.Algorithm;
import org.eclipse.elk.alg.test.framework.annotations.GraphProvider;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.util.ElkGraphUtil;
import org.junit.Ignore;
import org.junit.Test;
import org.junit.runner.RunWith;

/**
 * An example for ignored tests.
 */
@RunWith(LayoutTestRunner.class)
@Algorithm(ForceOptions.ALGORITHM_ID)
public class IgnoredTest {

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Sources

    /**
     * Supply a basic graph.
     */
    @GraphProvider
    public ElkNode firstGraph() {
        return ElkGraphUtil.createGraph();
    }
    
    /**
     * Supply a basic graph.
     */
    @GraphProvider
    public ElkNode secondGraph() {
        return ElkGraphUtil.createGraph();
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Configuration

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Tests

    /**
     * A simple test that always succeeds.
     */
    @Test
    public void success(final ElkNode graph) {
        assertTrue(true);
    }
    
    /**
     * A simple test that always fails.
     */
    @Test
    @Ignore("Test the ignore facilities")
    public void fail(final ElkNode graph) {
        assertTrue(false);
    }

}
