/*******************************************************************************
 * Copyright (c) 2022 Kiel University and others.
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

import java.util.List;

import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.test.framework.LayoutTestRunner;
import org.eclipse.elk.alg.test.framework.annotations.Algorithm;
import org.eclipse.elk.alg.test.framework.annotations.GraphResourceProvider;
import org.eclipse.elk.alg.test.framework.io.AbstractResourcePath;
import org.eclipse.elk.alg.test.framework.io.ModelResourcePath;
import org.eclipse.elk.graph.ElkNode;
import org.junit.Test;
import org.junit.runner.RunWith;

import com.google.common.collect.Lists;

/**
 * Tests whether the network simplex node placer correctly handles
 * a specific edge configuration.
 *
 */
@RunWith(LayoutTestRunner.class)
@Algorithm(LayeredOptions.ALGORITHM_ID)
public class Issue870Test {
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Sources
    
    @GraphResourceProvider
    public List<AbstractResourcePath> testGraphs() {
        return Lists.newArrayList(
               new ModelResourcePath("tickets/layered/870_networkSimplexNodePlacer.elkt"));
    }
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Tests
    
    @Test
    public void ensureLayoutDoesNotFail(final ElkNode graph) {
    }


}
