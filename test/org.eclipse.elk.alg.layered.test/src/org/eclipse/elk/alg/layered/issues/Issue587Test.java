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

import java.util.List;

import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.test.framework.LayoutTestRunner;
import org.eclipse.elk.alg.test.framework.annotations.Configurator;
import org.eclipse.elk.alg.test.framework.annotations.DefaultConfiguration;
import org.eclipse.elk.alg.test.framework.annotations.GraphResourceProvider;
import org.eclipse.elk.alg.test.framework.io.AbstractResourcePath;
import org.eclipse.elk.alg.test.framework.io.FileNameFilter;
import org.eclipse.elk.alg.test.framework.io.ModelResourcePath;
import org.eclipse.elk.core.options.EdgeRouting;
import org.eclipse.elk.graph.ElkNode;
import org.junit.Test;
import org.junit.runner.RunWith;

import com.google.common.collect.Lists;

/**
 * Test for issue 587.
 */
@RunWith(LayoutTestRunner.class)
@DefaultConfiguration()
public class Issue587Test {

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Sources

    @GraphResourceProvider
    public List<AbstractResourcePath> testGraphs() {
        return Lists.newArrayList(
                new ModelResourcePath("tickets/layered/**").withFilter(new FileNameFilter("587.+\\.elkt")));
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Configurations

    @Configurator
    public void configureOrthogonalEdgeRouting(final ElkNode graph) {
        graph.setProperty(LayeredOptions.EDGE_ROUTING, EdgeRouting.ORTHOGONAL);
    }

    @Configurator
    public void configurePolylineEdgeRouting(final ElkNode graph) {
        graph.setProperty(LayeredOptions.EDGE_ROUTING, EdgeRouting.POLYLINE);
    }

    @Configurator
    public void configureSplineEdgeRouting(final ElkNode graph) {
        graph.setProperty(LayeredOptions.EDGE_ROUTING, EdgeRouting.SPLINES);
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Tests
    
    @Test
    public void testNoException(final ElkNode graph) {
        // The bug caused an exception; nothing particular to test here
    }

}
