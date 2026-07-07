/*******************************************************************************
 * Copyright (c) 2019 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.test.framework.config;

import java.util.ArrayList;
import java.util.List;

import org.eclipse.elk.alg.test.framework.annotations.ConfiguratorProvider;
import org.eclipse.elk.alg.test.framework.util.TestUtil;
import org.eclipse.elk.core.LayoutConfigurator;
import org.eclipse.elk.core.util.ElkUtil;
import org.eclipse.elk.graph.ElkNode;
import org.junit.runners.model.FrameworkMethod;
import org.junit.runners.model.TestClass;

/**
 * A configuration which applies a {@link LayoutConfigurator} to a graph. The configurator is supplied by a method.
 * Instances of this class are generated whenever the {@link ConfiguratorProvider} annotation is encountered.
 */
public final class ConfiguratorTestConfiguration extends TestConfiguration {

    /** The method that supplies the layout configurator. */
    private final FrameworkMethod configuratorMethod;
    /** The layout configurator obtained from the method. Cached after first call. */
    private LayoutConfigurator configurator = null;

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // Creation

    private ConfiguratorTestConfiguration(final TestClass testClass, final FrameworkMethod configuratorMethod,
            final LayoutConfigurator configurator) {

        super(testClass);
        this.configuratorMethod = configuratorMethod;
        this.configurator = configurator;
    }

    /**
     * Loads the configurator supplier methods defined in the test class.
     * 
     * @param testClass
     *            the test class.
     * @param test
     *            instance of the test class to call methods.
     * @param errors
     *            list of error conditions encountered.
     * @return a list of configurations loaded from the test class.
     */
    public static List<ConfiguratorTestConfiguration> fromTestClass(final TestClass testClass, final Object test,
            final List<Throwable> errors) {

        List<ConfiguratorTestConfiguration> configurations = new ArrayList<>();

        for (FrameworkMethod method : testClass.getAnnotatedMethods(ConfiguratorProvider.class)) {
            TestUtil.ensurePublic(method, errors);
            TestUtil.ensureReturnsType(method, errors, LayoutConfigurator.class);
            TestUtil.ensureParameters(method, errors);

            // Invoke the method to obtain the layout configurator
            try {
                configurations.add(new ConfiguratorTestConfiguration(testClass, method,
                        (LayoutConfigurator) method.invokeExplosively(test)));
            } catch (Throwable e) {
                errors.add(e);
            }
        }

        return configurations;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // TestConfiguration

    @Override
    protected void applyBaseConfiguration(final Object testClass, final ElkNode graph) throws Throwable {
        ElkUtil.applyVisitors(graph, configurator);
    }

    @Override
    public String toString() {
        return super.toString() + " layoutConfigurator(" + configuratorMethod.getName() + ")";
    }

}
