/*******************************************************************************
 * Copyright (c) 2018, 2019 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.test.framework;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.eclipse.elk.alg.test.PlainJavaInitialization;
import org.eclipse.elk.alg.test.framework.algorithm.TestAlgorithm;
import org.eclipse.elk.alg.test.framework.annotations.OnlyOnRootNode;
import org.eclipse.elk.alg.test.framework.annotations.TestAfterProcessor;
import org.eclipse.elk.alg.test.framework.annotations.TestAfterProcessors;
import org.eclipse.elk.alg.test.framework.annotations.TestBeforeProcessor;
import org.eclipse.elk.alg.test.framework.annotations.TestBeforeProcessors;
import org.eclipse.elk.alg.test.framework.config.ConfiguratorTestConfiguration;
import org.eclipse.elk.alg.test.framework.config.MethodTestConfiguration;
import org.eclipse.elk.alg.test.framework.config.TestConfiguration;
import org.eclipse.elk.alg.test.framework.graph.GraphFromFile;
import org.eclipse.elk.alg.test.framework.graph.GraphFromMethod;
import org.eclipse.elk.alg.test.framework.graph.RandomGraphFromFile;
import org.eclipse.elk.alg.test.framework.graph.RandomGraphFromMethod;
import org.eclipse.elk.alg.test.framework.graph.TestGraph;
import org.eclipse.elk.alg.test.framework.util.TestUtil;
import org.eclipse.elk.core.AbstractLayoutProvider;
import org.eclipse.elk.core.alg.ILayoutProcessor;
import org.eclipse.elk.core.data.LayoutAlgorithmData;
import org.eclipse.elk.core.testing.IWhiteBoxTestable;
import org.eclipse.elk.graph.ElkNode;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import org.junit.runner.Description;
import org.junit.runner.RunWith;
import org.junit.runner.notification.RunNotifier;
import org.junit.runners.ParentRunner;
import org.junit.runners.model.FrameworkMethod;
import org.junit.runners.model.InitializationError;
import org.junit.runners.model.TestClass;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;

/**
 * The runner that manages the execution of all tests on the layout algorithms. This test runner has to be referenced in
 * the {@link RunWith} annotation of a test class. This runner is responsible for a few things:
 * 
 * <ul>
 * <li>It validates that the annotations make sense.</li>
 * <li>It loads all test conditions and produces the experimental objects (combination of layout algorithm,
 * configuration, and graph).</li>
 * <li>It collects all test methods and sorts them into blackbox and whitebox tests. These are collected once for the
 * child test runners to use.</li>
 * </ul>
 */
public class LayoutTestRunner extends ParentRunner<ExperimentRunner> {

    /** List of layout algorithms loaded from the test class. */
    private List<TestAlgorithm> testAlgorithms;
    /** List of layout configurations loaded from the test class. */
    private List<TestConfiguration> testConfigurations;
    /** List of graphs loaded from the test class. */
    private List<TestGraph> testGraphs;

    /** List of blackbox test methods. */
    private List<FrameworkMethod> blackboxTests;
    /** List of whitebox test methods. */
    private List<FrameworkMethod> whiteboxTests;
    /** Whitebox tests to be executed before a given layout processor. */
    private Multimap<Class<? extends ILayoutProcessor<?>>, FrameworkMethod> whiteboxBeforeTests;
    /** Whitebox tests to be executed after a given layout processor. */
    private Multimap<Class<? extends ILayoutProcessor<?>>, FrameworkMethod> whiteboxAfterTests;
    /** Set of whitebox test methods that only want to be executed on an input graph's root node. */
    private Set<FrameworkMethod> whiteboxOnlyOnRoot;

    /** All of the test runners we'll be running. */
    private final List<ExperimentRunner> childRunners = new ArrayList<>();

    /**
     * Creates a new instance for the given test class.
     * 
     * @param testClass
     *            the test class to be executed with this test runner.
     * @throws InitializationError
     *             represents an error that occurred during initialization
     */
    public LayoutTestRunner(final Class<?> testClass) throws InitializationError {
        super(testClass);

        // At this point, all configuration data have been loaded by collectInitializationErrors(...) and we have made
        // an effort to detect faulty configurations. Build our child runners
        List<ExperimentalObject> experimentalObjects =
                ExperimentalObject.inflate(testAlgorithms, testConfigurations, testGraphs);

        for (ExperimentalObject eo : experimentalObjects) {
            childRunners.add(new ExperimentRunner(this, eo));
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // Getters

    /**
     * Returns the list of blackbox tests.
     */
    public List<FrameworkMethod> getBlackboxTests() {
        return blackboxTests;
    }

    /**
     * Returns a list of all whitebox tests.
     */
    public List<FrameworkMethod> getWhiteboxTests() {
        return whiteboxTests;
    }

    /**
     * Returns a map that maps layout processors to all test methods that want to be executed before that processor.
     */
    public Multimap<Class<? extends ILayoutProcessor<?>>, FrameworkMethod> getWhiteboxBeforeTests() {
        return whiteboxBeforeTests;
    }

    /**
     * Returns a map that maps layout processors to all test methods that want to be executed after that processor.
     */
    public Multimap<Class<? extends ILayoutProcessor<?>>, FrameworkMethod> getWhiteboxAfterTests() {
        return whiteboxAfterTests;
    }

    /**
     * Returns the set of whitebox test that only want to be run on a root graph.
     */
    public Set<FrameworkMethod> getWhiteboxOnlyOnRoot() {
        return whiteboxOnlyOnRoot;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // ParentRunner Initialization

    @Override
    protected void collectInitializationErrors(final List<Throwable> errors) {
        super.collectInitializationErrors(errors);

        TestClass testClass = getTestClass();

        // We're in plain Java mode, so ensure everything's initialized properly
        PlainJavaInitialization.initializePlainJavaLayout();

        initializeFields();

        // Initialize our test methods
        initializeBlackboxTests(errors);
        initializeWhiteboxTests(errors);

        // Obtain an instance of our test class for initialization code to call
        Object test = null;
        if (!TestUtil.ensureSinglePublicNoParameterConstructor(testClass, errors)) {
            return;
        }

        try {
            test = TestUtil.createTestClassInstance(testClass);
        } catch (Throwable e) {
            errors.add(e);
            return;
        }

        // Initialize test data
        initializeGraphs(testClass, test, errors);
        initializeAlgorithms(testClass, errors);
        initializeConfigurations(testClass, test, errors);
    }

    /**
     * We need to initialize some of our fields here because the initialization methods are called during construction
     * of the superclass, which means that our fields wouldn't have been initialized yet.
     */
    private void initializeFields() {
        testAlgorithms = new ArrayList<>();
        testConfigurations = new ArrayList<>();
        testGraphs = new ArrayList<>();

        blackboxTests = new ArrayList<>();

        whiteboxTests = new ArrayList<>();
        whiteboxBeforeTests = HashMultimap.create();
        whiteboxAfterTests = HashMultimap.create();
        whiteboxOnlyOnRoot = new HashSet<>();
    }

    /**
     * Collects and validates all blackbox test methods.
     */
    private void initializeBlackboxTests(final List<Throwable> errors) {
        // We need to initialize our data structures here because this method is called during construction of the
        // superclass, which means that our fields wouldn't have been initialized yet

        for (FrameworkMethod test : getTestClass().getAnnotatedMethods(Test.class)) {
            TestUtil.ensurePublic(test, errors);
            TestUtil.ensureParameters(test, errors, ElkNode.class);
            TestUtil.ensureNotAnnotatedWith(test, errors, TestBeforeProcessor.class, TestAfterProcessor.class,
                    BeforeClass.class, Before.class, After.class, AfterClass.class, OnlyOnRootNode.class);

            blackboxTests.add(test);
        }
    }

    /**
     * Collects and validates all whitebox test methods.
     */
    private void initializeWhiteboxTests(final List<Throwable> errors) {
        Set<FrameworkMethod> encounteredTestMethods = new HashSet<>();
        
        // Tests that want to be executed before a given layout processor
        for (FrameworkMethod test : getTestClass().getAnnotatedMethods(TestBeforeProcessors.class)) {
            initializeWhiteboxTest(errors, test, encounteredTestMethods);

            // There are multiple annotations
            for (TestBeforeProcessor annotation : test.getMethod().getAnnotationsByType(TestBeforeProcessor.class)) {
                whiteboxBeforeTests.put(annotation.value(), test);
            }
        }
        
        for (FrameworkMethod test : getTestClass().getAnnotatedMethods(TestBeforeProcessor.class)) {
            initializeWhiteboxTest(errors, test, encounteredTestMethods);
            whiteboxBeforeTests.put(test.getMethod().getAnnotationsByType(TestBeforeProcessor.class)[0].value(), test);
        }

        // Tests that want to be executed after a given layout processor
        for (FrameworkMethod test : getTestClass().getAnnotatedMethods(TestAfterProcessors.class)) {
            initializeWhiteboxTest(errors, test, encounteredTestMethods);

            // There are multiple annotations
            for (TestAfterProcessor annotation : test.getMethod().getAnnotationsByType(TestAfterProcessor.class)) {
                whiteboxAfterTests.put(annotation.value(), test);
            }
        }
        
        for (FrameworkMethod test : getTestClass().getAnnotatedMethods(TestAfterProcessor.class)) {
            initializeWhiteboxTest(errors, test, encounteredTestMethods);
            whiteboxAfterTests.put(test.getMethod().getAnnotationsByType(TestAfterProcessor.class)[0].value(), test);
        }
    }

    /**
     * Performs some general initializations for whitebox tests.
     */
    private void initializeWhiteboxTest(final List<Throwable> errors, final FrameworkMethod testMethod,
            final Set<FrameworkMethod> encounteredTestMethods) {
        
        TestUtil.ensurePublic(testMethod, errors);
        TestUtil.ensureParameters(testMethod, errors, Object.class);
        TestUtil.ensureNotAnnotatedWith(testMethod, errors, Test.class, TestAfterProcessor.class, BeforeClass.class,
                Before.class, After.class, AfterClass.class);

        // Check if the test only wants to be run on the root node
        if (testMethod.getAnnotation(OnlyOnRootNode.class) != null) {
            whiteboxOnlyOnRoot.add(testMethod);
        }
        
        // Add to our list of whitebox tests if it isn't already there
        if (encounteredTestMethods.add(testMethod)) {
            whiteboxTests.add(testMethod);
        }
    }

    /**
     * Loads all graphs configured in the test class.
     */
    private void initializeGraphs(final TestClass testClass, final Object test, final List<Throwable> errors) {
        testGraphs.addAll(GraphFromFile.fromTestClass(testClass, test, errors));
        testGraphs.addAll(GraphFromMethod.fromTestClass(testClass, test, errors));
        testGraphs.addAll(RandomGraphFromFile.fromTestClass(testClass, test, errors));
        testGraphs.addAll(RandomGraphFromMethod.fromTestClass(testClass, test, errors));

        if (testGraphs.isEmpty()) {
            // There's no need to go any further if there are no graphs
            errors.add(new Exception("No test graphs specified."));
        }
    }

    /**
     * Loads all algorithms configured in the test class.
     */
    private void initializeAlgorithms(final TestClass testClass, final List<Throwable> errors) {
        testAlgorithms.addAll(TestAlgorithm.fromTestClass(testClass, errors));

        // If we have whitebox tests, we need to check stuff
        if (!whiteboxTests.isEmpty()) {
            // Ensure that there is at least one specific layout algorithm
            if (testAlgorithms.isEmpty() || testAlgorithms.get(0).getAlgorithmData() == null) {
                errors.add(new Exception("Whitebox tests require explicit @Algorithm annotations."));

            } else {
                // Ensure that all algorithms are whitebox testable
                for (TestAlgorithm algorithm : testAlgorithms) {
                    LayoutAlgorithmData algorithmData = algorithm.getAlgorithmData();

                    AbstractLayoutProvider layoutProvider = algorithmData.getInstancePool().fetch();

                    if (!(layoutProvider instanceof IWhiteBoxTestable)) {
                        errors.add(new Exception("Algorithm " + algorithmData.getId() + " is not whitebox testable."));
                    }

                    algorithmData.getInstancePool().release(layoutProvider);
                }
            }
        }

        if (testAlgorithms.isEmpty()) {
            testAlgorithms.add(TestAlgorithm.identity());
        }
    }

    /**
     * Loads all configurations configured in the test class.
     */
    private void initializeConfigurations(final TestClass testClass, final Object test, final List<Throwable> errors) {
        testConfigurations.addAll(ConfiguratorTestConfiguration.fromTestClass(testClass, test, errors));
        testConfigurations.addAll(MethodTestConfiguration.fromTestClass(testClass, errors));

        if (testConfigurations.isEmpty()) {
            testConfigurations.add(TestConfiguration.fromTestClass(testClass));
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // ParentRunner Child Management

    @Override
    protected List<ExperimentRunner> getChildren() {
        return childRunners;
    }

    @Override
    protected Description describeChild(final ExperimentRunner child) {
        return child.getDescription();
    }

    @Override
    protected void runChild(final ExperimentRunner child, final RunNotifier notifier) {
        child.run(notifier);
    }

}
