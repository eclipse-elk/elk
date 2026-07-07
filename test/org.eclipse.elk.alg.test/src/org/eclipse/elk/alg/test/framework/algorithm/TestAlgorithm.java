/*******************************************************************************
 * Copyright (c) 2019 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.test.framework.algorithm;

import java.util.ArrayList;
import java.util.List;

import org.eclipse.elk.alg.test.framework.annotations.Algorithm;
import org.eclipse.elk.alg.test.framework.annotations.AllAlgorithms;
import org.eclipse.elk.core.LayoutConfigurator;
import org.eclipse.elk.core.data.LayoutAlgorithmData;
import org.eclipse.elk.core.data.LayoutMetaDataService;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.util.ElkUtil;
import org.eclipse.elk.graph.ElkNode;
import org.junit.runners.model.TestClass;

/**
 * Configures graphs to use a specific layout algorithm.
 */
public final class TestAlgorithm {
    
    /** The layout algorithm to apply to graphs. */
    private final LayoutAlgorithmData algorithm;
    /** The layout configurator that will be used for the actual application. */
    private final LayoutConfigurator algorithmConfigurator;

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // Creation
    
    /**
     * Creates a new instance that configures graphs to use the given layout algorithm.
     */
    private TestAlgorithm(final LayoutAlgorithmData algorithm) {
        this.algorithm = algorithm;
        
        if (algorithm != null) {
            algorithmConfigurator = new LayoutConfigurator();
            algorithmConfigurator.configure(ElkNode.class).setProperty(CoreOptions.ALGORITHM, algorithm.getId());
        } else {
            algorithmConfigurator = null;
        }
    }
    
    /**
     * Loads the test layout algorithms defined in the test class.
     * 
     * @param testClass
     *            the test class.
     * @param errors
     *            list of error conditions encountered while evaluating the layout algorithms.
     * @return a list of layout algorithms to use during the test.
     */
    public static List<TestAlgorithm> fromTestClass(final TestClass testClass, final List<Throwable> errors) {
        List<LayoutAlgorithmData> algorithmData = new ArrayList<>();

        // Obtain annotations
        Algorithm[] algorithmAnnotations = testClass.getJavaClass().getAnnotationsByType(Algorithm.class);
        AllAlgorithms allAlgorithmsAnnotation = testClass.getAnnotation(AllAlgorithms.class);

        LayoutMetaDataService service = LayoutMetaDataService.getInstance();

        if (algorithmAnnotations.length != 0 && allAlgorithmsAnnotation != null) {
            // Only one of the two may be specified
            errors.add(new Exception("If @AllAlgorithms is specified, @Algorithm cannot be used anymore."));

        } else if (allAlgorithmsAnnotation != null) {
            // Simply add all known algorithm IDs
            algorithmData.addAll(service.getAlgorithmData());

        } else {
            // Add all specifically supplied IDs as long as they refer to known layout algorithms
            for (Algorithm algAnnotation : algorithmAnnotations) {
                LayoutAlgorithmData algData = service.getAlgorithmData(algAnnotation.value());
                if (algData != null) {
                    algorithmData.add(algData);
                } else {
                    errors.add(new Exception("Unknown layout algorithm: " + algAnnotation.value()));
                }
            }
        }

        // Produce the list of algorithms
        List<TestAlgorithm> algorithms = new ArrayList<>();

        algorithmData.stream()
                .sorted((data1, data2) -> data1.getId().compareTo(data2.getId()))
                .map(algData -> new TestAlgorithm(algData))
                .forEach(strategy -> algorithms.add(strategy));

        return algorithms;
    }
    
    /**
     * Returns a test algorithm that simply leaves the input graphs unchanged.
     */
    public static TestAlgorithm identity() {
        return new TestAlgorithm(null);
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // Getters
    
    /**
     * Returns the layout algorithm or {@code null} if this instance doesn't change layout algorithms.
     */
    public LayoutAlgorithmData getAlgorithmData() {
        return algorithm;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // TestAlgorithm
    
    /**
     * Applies the strategy to the given graph.
     * 
     * @param test
     *            instance of the test class that runs test methods.
     * @param graph
     *            the graph to apply the strategy to..
     */
    public void apply(final Object test, final ElkNode graph) {
        if (algorithmConfigurator != null) {
            ElkUtil.applyVisitors(graph, algorithmConfigurator);
        }
    }

    @Override
    public String toString() {
        if (algorithm == null) {
            return "algorithm(as in test graph)";
        } else {
            return "algorithm(" + algorithm.getId() + ")";
        }
    }

}
