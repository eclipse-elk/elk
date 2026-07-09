/*******************************************************************************
 * Copyright (c) 2019 Kiel University and others.
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
package org.eclipse.elk.alg.test.framework.graph;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.eclipse.elk.alg.test.framework.annotations.RandomGeneratorOptions;
import org.eclipse.elk.alg.test.framework.util.TestUtil;
import org.eclipse.elk.core.debug.grandom.generators.GeneratorOptions;
import org.eclipse.elk.core.debug.grandom.generators.RandomGraphGenerator;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.emf.ecore.util.EcoreUtil;
import org.junit.runners.model.FrameworkMethod;
import org.junit.runners.model.TestClass;

/**
 * A test graph that generates input graphs using the {@link RandomGraphGenerator}. The configuration for the generator
 * is supplied by a configuration method. Instances of this class are generated whenever the
 * {@link RandomGeneratorOptions} annotation is encountered.
 */
public final class RandomGraphFromMethod extends TestGraph {

    /** Seed for the random graph generator to make tests deterministic. */
    private static final long RANDOM_SEED = 1337;

    /** The method that supplied the random graph configuration. */
    private final FrameworkMethod configuratorMethod;
    /** Since a file can generate several random graphs, this is the one this instance will supply. */
    private final ElkNode graph;

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // Creation

    private RandomGraphFromMethod(final FrameworkMethod configuratorMethod, final ElkNode graph) {
        this.configuratorMethod = configuratorMethod;
        this.graph = graph;
    }

    /**
     * Loads random graphs from methods that specify random graph configurations.
     * 
     * @param testClass
     *            the test class.
     * @param test
     *            instance of the test class to call methods.
     * @param errors
     *            list of error conditions encountered.
     * @return a list of configurations loaded from the test class.
     */
    public static List<RandomGraphFromMethod> fromTestClass(final TestClass testClass, final Object test,
            final List<Throwable> errors) {

        List<RandomGraphFromMethod> graphs = new ArrayList<>();

        // We use a fixed seed to keep tests deterministic
        RandomGraphGenerator generator = new RandomGraphGenerator(new Random(RANDOM_SEED));

        for (FrameworkMethod method : testClass.getAnnotatedMethods(RandomGeneratorOptions.class)) {
            TestUtil.ensurePublic(method, errors);
            TestUtil.ensureReturnsType(method, errors, GeneratorOptions.class);
            TestUtil.ensureParameters(method, errors);

            // Invoke the method to obtain random graphs
            try {
                GeneratorOptions options = (GeneratorOptions) method.invokeExplosively(test);
                graphs.add(new RandomGraphFromMethod(method, generator.generate(options)));
            } catch (Throwable e) {
                errors.add(e);
            }
        }

        return graphs;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // TestGraph

    @Override
    public ElkNode provideGraph(final Object test) {
        // Simply return a copy of the graph
        return EcoreUtil.copy(graph);
    }

    @Override
    public String toString() {
        return "randomGraphMethod(" + configuratorMethod.getName() + ")";
    }

}
