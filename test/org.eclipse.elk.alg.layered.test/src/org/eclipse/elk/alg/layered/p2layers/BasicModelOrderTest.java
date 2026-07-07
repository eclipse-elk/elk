/*******************************************************************************
 * Copyright (c) 2023 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 * 
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later 
 *******************************************************************************/
package org.eclipse.elk.alg.layered.p2layers;

import java.util.List;

import org.eclipse.elk.alg.layered.options.CycleBreakingStrategy;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.layered.options.LayeringStrategy;
import org.eclipse.elk.alg.test.framework.LayoutTestRunner;
import org.eclipse.elk.alg.test.framework.annotations.Algorithm;
import org.eclipse.elk.alg.test.framework.annotations.ConfiguratorProvider;
import org.eclipse.elk.alg.test.framework.annotations.DefaultConfiguration;
import org.eclipse.elk.alg.test.framework.annotations.GraphResourceProvider;
import org.eclipse.elk.alg.test.framework.io.AbstractResourcePath;
import org.eclipse.elk.alg.test.framework.io.FileExtensionFilter;
import org.eclipse.elk.alg.test.framework.io.ModelResourcePath;
import org.eclipse.elk.core.LayoutConfigurator;
import org.eclipse.elk.graph.ElkNode;
import org.junit.Test;
import org.junit.runner.RunWith;

import com.google.common.collect.Lists;

@RunWith(LayoutTestRunner.class)
@Algorithm(LayeredOptions.ALGORITHM_ID)
@DefaultConfiguration()
public class BasicModelOrderTest {

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Sources

    @GraphResourceProvider
    public List<AbstractResourcePath> testGraphs() {
        return Lists.newArrayList(
                new ModelResourcePath("realworld/ptolemy/**/").withFilter(new FileExtensionFilter("elkg")));
    }
    

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Configuration

    @ConfiguratorProvider
    public LayoutConfigurator bfConfiguration() {
        LayoutConfigurator config = new LayoutConfigurator();
        config.configure(ElkNode.class).setProperty(LayeredOptions.CYCLE_BREAKING_STRATEGY,
                CycleBreakingStrategy.MODEL_ORDER);
        config.configure(ElkNode.class).setProperty(LayeredOptions.LAYERING_STRATEGY,
                LayeringStrategy.BF_MODEL_ORDER);
        return config;
    }

    @ConfiguratorProvider
    public LayoutConfigurator dfConfiguration() {
        LayoutConfigurator config = new LayoutConfigurator();
        config.configure(ElkNode.class).setProperty(LayeredOptions.CYCLE_BREAKING_STRATEGY,
                CycleBreakingStrategy.MODEL_ORDER);
        config.configure(ElkNode.class).setProperty(LayeredOptions.LAYERING_STRATEGY,
                LayeringStrategy.DF_MODEL_ORDER);
        return config;
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Tests
    
    // Just check for errors that might occur
    
    @Test
    public void test(final Object graph) {
        assert(true);
    }
}
