/*******************************************************************************
 * Copyright (c) 2019 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.test.framework;

import java.util.ArrayList;
import java.util.List;

import org.eclipse.elk.alg.test.framework.algorithm.TestAlgorithm;
import org.eclipse.elk.alg.test.framework.config.TestConfiguration;
import org.eclipse.elk.alg.test.framework.graph.TestGraph;
import org.eclipse.elk.graph.ElkNode;

/**
 * An experimental object is a concrete tuple of a layout algorithm, a configuration, and a graph. The experimental
 * object does not keep specific graph instances. Instead, it just keeps everything necessary to know how to create a
 * graph. Clients call {@link #realize()} to obtain a concrete graph to run layout on.
 */
final class ExperimentalObject {

    /**
     * The layout algorithm strategy to be applied to our graph. This can be to either leave the graph's default
     * settings alone, or to set it to a particular layout algorithm.
     */
    private final TestAlgorithm layoutAlgorithm;

    /**
     * The configuration to apply to the graph. This can, again, be to leave the graph alone or something more
     * elaborate.
     */
    private final TestConfiguration configuration;

    /**
     * A way to get our hands at a graph.
     */
    private final TestGraph graph;

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // Construction

    /**
     * Creates a new instance for the given strategies.
     * 
     * @param layoutAlgorithm
     *            layout algorithm.
     * @param configuration
     *            configurator.
     * @param graph
     *            graph.
     */
    ExperimentalObject(final TestAlgorithm layoutAlgorithm, final TestConfiguration configuration,
            final TestGraph graph) {

        this.layoutAlgorithm = layoutAlgorithm;
        this.configuration = configuration;
        this.graph = graph;
    }

    /**
     * Creates experimental objects for all combinations of algorithms, configurators, and graphs.
     * 
     * @param algorithms
     *            list of algorithm strategies.
     * @param configurators
     *            list of configuration strategies.
     * @param graphs
     *            list of graphs.
     * @return list of experimental objects.
     */
    public static List<ExperimentalObject> inflate(final List<TestAlgorithm> algorithms,
            final List<TestConfiguration> configurators, final List<TestGraph> graphs) {

        // Create a list large enough
        List<ExperimentalObject> objects = new ArrayList<>(algorithms.size() * configurators.size() * graphs.size());

        for (TestAlgorithm algorithm : algorithms) {
            for (TestConfiguration configurator : configurators) {
                for (TestGraph graph : graphs) {
                    objects.add(new ExperimentalObject(algorithm, configurator, graph));
                }
            }
        }

        return objects;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // Realization

    /**
     * Produces a graph for this experimental object.
     * 
     * @param test
     *            instance of the test class that runs test methods.
     * @return a graph.
     * @throws Throwable
     *             if something goes wrong.
     */
    public ElkNode realize(final Object test) throws Throwable {
        ElkNode testGraph = graph.provideGraph(test);
        layoutAlgorithm.apply(test, testGraph);
        configuration.apply(test, testGraph);

        return testGraph;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // Getters

    /**
     * Returns the experimental object's layout algorithm strategy.
     */
    public TestAlgorithm getLayoutAlgorithm() {
        return layoutAlgorithm;
    }

    /**
     * Returns the experimental object's configuration strategy.
     */
    public TestConfiguration getConfiguration() {
        return configuration;
    }

    /**
     * Returns the experimental object's graph loading strategy.
     */
    public TestGraph getGraph() {
        return graph;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // Object

    @Override
    public String toString() {
        return layoutAlgorithm + " " + configuration + " " + graph;
    }

}
