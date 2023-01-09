/*******************************************************************************
 * Copyright (c) 2020 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.shared.test;

import static org.junit.Assert.assertTrue;

import java.util.List;

import org.eclipse.elk.alg.test.framework.LayoutTestRunner;
import org.eclipse.elk.alg.test.framework.annotations.DefaultConfiguration;
import org.eclipse.elk.alg.test.framework.annotations.GraphResourceProvider;
import org.eclipse.elk.alg.test.framework.io.AbstractResourcePath;
import org.eclipse.elk.alg.test.framework.io.FileExtensionFilter;
import org.eclipse.elk.alg.test.framework.io.ModelResourcePath;
import org.eclipse.elk.graph.ElkNode;
import org.junit.Test;
import org.junit.runner.RunWith;

import com.google.common.collect.Lists;

/**
 * Test ensuring that the example models are laid out properly.
 */
@RunWith(LayoutTestRunner.class)
@DefaultConfiguration
public class ElkLiveExamplesTest {

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Sources

    @GraphResourceProvider
    public List<AbstractResourcePath> testGraphs() {
        return Lists.newArrayList(new ModelResourcePath("examples/**/").withFilter(new FileExtensionFilter("elkt")));
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Configuration

    @Test
    public void noException(final ElkNode graph) {
        // Nothing to test in particular, simply ensure that no exception is raised.
        assertTrue(true);
    }

}
