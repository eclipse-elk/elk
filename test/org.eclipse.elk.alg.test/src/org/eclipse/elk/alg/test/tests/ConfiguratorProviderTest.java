/*******************************************************************************
 * Copyright (c) 2018, 2019 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.test.tests;

import static org.junit.Assert.assertEquals;

import org.eclipse.elk.alg.force.options.ForceOptions;
import org.eclipse.elk.alg.test.framework.LayoutTestRunner;
import org.eclipse.elk.alg.test.framework.annotations.Algorithm;
import org.eclipse.elk.alg.test.framework.annotations.ConfiguratorProvider;
import org.eclipse.elk.alg.test.framework.annotations.GraphProvider;
import org.eclipse.elk.core.LayoutConfigurator;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.properties.IProperty;
import org.eclipse.elk.graph.properties.Property;
import org.eclipse.elk.graph.util.ElkGraphUtil;
import org.junit.Test;
import org.junit.runner.RunWith;

/**
 * An example for using a layout configurator provider to configure the layout graph.
 */
@RunWith(LayoutTestRunner.class)
@Algorithm(ForceOptions.ALGORITHM_ID)
public class ConfiguratorProviderTest {

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
    // Configuration

    // A property to test with
    private static final IProperty<String> TEST_PROPERTY = new Property<>("org.eclipse.elk.alg.common.test", "");
    private static final String TEST_PROPERTY_VALUE = "Test";

    /**
     * Supplies a configurator that sets our test property on nodes.
     */
    @ConfiguratorProvider
    public LayoutConfigurator configuratorProvider() {
        LayoutConfigurator layoutConfig = new LayoutConfigurator();
        layoutConfig.configure(ElkNode.class).setProperty(TEST_PROPERTY, TEST_PROPERTY_VALUE);
        return layoutConfig;
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Tests

    /**
     * Check that the property was set on the graph.
     */
    @Test
    public void testConfiguratorProvider(final ElkNode graph) {
        assertEquals("Property is not set properly.", TEST_PROPERTY_VALUE, graph.getProperty(TEST_PROPERTY));
    }

}
