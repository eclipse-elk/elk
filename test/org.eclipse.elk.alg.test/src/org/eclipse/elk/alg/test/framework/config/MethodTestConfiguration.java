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

import org.eclipse.elk.alg.test.framework.annotations.Configurator;
import org.eclipse.elk.alg.test.framework.util.TestUtil;
import org.eclipse.elk.graph.ElkNode;
import org.junit.runners.model.FrameworkMethod;
import org.junit.runners.model.TestClass;

/**
 * A configuration which calls a method to configure to a graph. Instances of this class are generated whenever the
 * {@link Configurator} annotation is encountered.
 */
public final class MethodTestConfiguration extends TestConfiguration {

    /** The method that configures the graph. */
    private final FrameworkMethod configuratorMethod;

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // Creation

    private MethodTestConfiguration(final TestClass testClass, final FrameworkMethod configuratorMethod) {
        super(testClass);
        this.configuratorMethod = configuratorMethod;
    }

    /**
     * Loads the configuration methods defined in the test class.
     * 
     * @param testClass
     *            the test class.
     * @param errors
     *            list of error conditions encountered.
     * @return a list of configurations loaded from the test class.
     */
    public static List<MethodTestConfiguration> fromTestClass(final TestClass testClass, final List<Throwable> errors) {
        List<MethodTestConfiguration> configurations = new ArrayList<>();

        for (FrameworkMethod method : testClass.getAnnotatedMethods(Configurator.class)) {
            TestUtil.ensurePublic(method, errors);
            TestUtil.ensureReturnsType(method, errors, Void.TYPE);
            TestUtil.ensureParameters(method, errors, ElkNode.class);

            configurations.add(new MethodTestConfiguration(testClass, method));
        }

        return configurations;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // TestConfiguration

    @Override
    protected void applyBaseConfiguration(final Object test, final ElkNode graph) throws Throwable {
        configuratorMethod.invokeExplosively(test, graph);
    }

    @Override
    public String toString() {
        return super.toString() + " configMethod(" + configuratorMethod.getName() + ")";
    }

}
