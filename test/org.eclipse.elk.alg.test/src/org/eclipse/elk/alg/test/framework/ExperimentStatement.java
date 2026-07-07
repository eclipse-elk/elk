/*******************************************************************************
 * Copyright (c) 2019 le-cds and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.test.framework;

import java.util.Collection;

import org.eclipse.elk.alg.test.framework.algorithm.TestAlgorithm;
import org.eclipse.elk.alg.test.framework.annotations.FailIfNotExecuted;
import org.eclipse.elk.alg.test.framework.util.TestUtil;
import org.eclipse.elk.core.RecursiveGraphLayoutEngine;
import org.eclipse.elk.core.alg.ILayoutProcessor;
import org.eclipse.elk.core.testing.TestController;
import org.eclipse.elk.core.testing.TestController.ILayoutExecutionListener;
import org.eclipse.elk.core.util.BasicProgressMonitor;
import org.eclipse.elk.graph.ElkNode;
import org.junit.runner.Description;
import org.junit.runner.notification.Failure;
import org.junit.runner.notification.RunNotifier;
import org.junit.runners.model.FrameworkMethod;
import org.junit.runners.model.Statement;

import com.google.common.collect.Multimap;

/**
 * Statement in charge of running the actual experiment. This statement works in several steps:
 * 
 * <ol>
 * <li>If there are whitebox tests, create a {@link TestController} to be installed on the layout algorithm.</li>
 * <li>Run automatic layout through the {@link RecursiveGraphLayoutEngine}. If there are whitebox tests, this will
 * cause the algorithm to tell us whenever there is an opportunity for whitebox tests to run.</li>
 * <li>Whitebox tests are run, but instead of running them like blackbox tests, we record any results without
 * notifying the run notifiers just yet. This is because the notifiers should be notified for each test method
 * exactly once. The recording happens through {@link WhiteboxRecorder}, a custom run notifier.</li>
 * <li>After layout, we run all blackbox tests through the usual channels.</li>
 * <li>Now that we have all the data, we can simulate whitebox test runs for the notifiers. Based on whether the
 * tests were executed (and with what results) we will simply call the appropriate notifier methods, possibly with
 * the failures we have since collected.</li>
 * </ol>
 */
final class ExperimentStatement extends Statement implements ILayoutExecutionListener {

    /** The runner that created us. */
    private final ExperimentRunner experimentRunner;
    /** The notifier we'll be notifying. */
    private final RunNotifier notifier;
    /** A special notifier that records results of whitebox tests. */
    private final WhiteboxRecorder whiteboxRecorder = new WhiteboxRecorder();
    /** The concrete graph for the experimental object. Created when this runner is started. */
    private ElkNode testGraph;

    /**
     * Creates a new statement for the given runner that will notify the given notifier.
     */
    ExperimentStatement(final ExperimentRunner experimentRunner, final RunNotifier notifier) {
        this.experimentRunner = experimentRunner;
        this.notifier = notifier;
    }

    @Override
    public void evaluate() throws Throwable {
        createTestGraph();
        
        // Create a test controller if we have whitebox tests and the layout algorithm
        TestController testController = null;

        if (!this.experimentRunner.getParentRunner().getWhiteboxTests().isEmpty()) {
            TestAlgorithm algorithm = experimentRunner.getExperimentalObject().getLayoutAlgorithm();
            testController = new TestController(algorithm.getAlgorithmData().getId());
            testController.addLayoutExecutionListener(this);
        }

        // Let the recursive graph layout engine run the layout algorithm (if we have a test controller, that will
        // automatically cause white box tests to run)
        RecursiveGraphLayoutEngine layoutEngine = new RecursiveGraphLayoutEngine();
        layoutEngine.layout(testGraph, testController, new BasicProgressMonitor());

        // Stop listening to the test controller
        if (testController != null) {
            testController.removeLayoutExecutionListener(this);
        }

        // Run blackbox tests now that the layout algorithm has finished
        runBlackboxTests();

        // Cleanly handle all whitebox tests that did not execute because their processors didn't run
        simulateWhiteboxTestsForNotifiers();
        
        testGraph = null;
    }

    /**
     * Turns the experimental object into an actual test graph.
     */
    private void createTestGraph() throws Throwable {
        Object testInstance = TestUtil.createTestClassInstance(this.experimentRunner.getTestClass());
        testGraph = this.experimentRunner.getExperimentalObject().realize(testInstance);
    }

    /**
     * Runs all blackbox tests after the layout algorithm has finished. By this time, all whitebox tests will have
     * finished.
     */
    private void runBlackboxTests() {
        for (FrameworkMethod bbTest : experimentRunner.getParentRunner().getBlackboxTests()) {
            runTest(bbTest, testGraph, notifier);
        }
    }

    /**
     * Runs whitebox tests on the given processor.
     * 
     * @param wbTests
     *            the collection of whitebox tests. This will be whitebox tests to be run either before or after
     *            processors.
     * @param processor
     *            the processor about to be executed or just executed.
     * @param graph
     *            the graph the processor ran on.
     * @param isRoot
     *            whether the graph is a root graph.
     */
    private void runWhiteboxTests(final Multimap<Class<? extends ILayoutProcessor<?>>, FrameworkMethod> wbTests,
            final ILayoutProcessor<?> processor, final Object graph, final boolean isRoot) {

        @SuppressWarnings("unchecked")
        Collection<FrameworkMethod> wbTestsForProcessor =
                wbTests.get((Class<? extends ILayoutProcessor<?>>) processor.getClass());

        for (FrameworkMethod wbTest : wbTestsForProcessor) {
            // Only run the test if this is either the root graph (all whitebox tests want to be run on the root
            // graph) or if the test is not restricted to run only on the root
            if (isRoot || !experimentRunner.getParentRunner().getWhiteboxOnlyOnRoot().contains(wbTest)) {
                // Run the test, but not "properly"... we shall only record the test results
                runTest(wbTest, graph, whiteboxRecorder);
            }
        }
    }

    /**
     * Runs a test method and notifies the given notifier of any results.
     * 
     * @param test
     *            the test to run.
     * @param graph
     *            the graph to be passed to the test method.
     * @param actualNotifier
     *            the notifier to notify of any results. This will be our usual notifier for blackbox test and a
     *            special notifier for whitebox test that records the results.
     */
    private void runTest(final FrameworkMethod test, final Object graph, final RunNotifier actualNotifier) {
        Description testDescription = this.experimentRunner.describeChild(test);
        
        // Only run the test if it's not currently ignored
        if (this.experimentRunner.isIgnored(test)) {
            actualNotifier.fireTestIgnored(testDescription);
        } else {
            Statement bbTestStatement = new Statement() {
                @Override
                public void evaluate() throws Throwable {
                    Object testInstance = TestUtil.createTestClassInstance(experimentRunner.getTestClass());
                    test.invokeExplosively(testInstance, graph);
                }
            };

            experimentRunner.doRunLeaf(bbTestStatement, testDescription, actualNotifier);
        }
    }

    /**
     * Runs through each whitebox test and simulates a test run for the notifiers based on the results recorded.
     */
    private void simulateWhiteboxTestsForNotifiers() {
        for (FrameworkMethod wbTest : experimentRunner.getParentRunner().getWhiteboxTests()) {
            Description wbTestDescription = this.experimentRunner.describeChild(wbTest);
            
            if (this.experimentRunner.isIgnored(wbTest)) {
                notifier.fireTestIgnored(wbTestDescription);
                
            } else {
                // We definitely have to simulate a test run
                notifier.fireTestStarted(wbTestDescription);
                
                if (whiteboxRecorder.wasExecuted(wbTestDescription)) {
                    // The test ran at least once. Fail if we have recorded any failures
                    for (Failure failure : whiteboxRecorder.getFailures(wbTestDescription)) {
                        notifier.fireTestFailure(failure);
                    }
                    
                } else {
                    // The test didn't run. Fire a failure if it wanted to run
                    if (wbTest.getAnnotation(FailIfNotExecuted.class) != null) {
                        notifier.fireTestFailure(new Failure(wbTestDescription,
                                new Exception("Test was not executed, but was supposed to.")));
                    }
                }
                
                notifier.fireTestFinished(wbTestDescription);
            }
        }
    }

    @Override
    public void layoutProcessorReady(final ILayoutProcessor<?> processor, final Object graph,
            final boolean isRoot) {

        runWhiteboxTests(experimentRunner.getParentRunner().getWhiteboxBeforeTests(), processor, graph, isRoot);
    }

    @Override
    public void layoutProcessorFinished(final ILayoutProcessor<?> processor, final Object graph,
            final boolean isRoot) {

        runWhiteboxTests(experimentRunner.getParentRunner().getWhiteboxAfterTests(), processor, graph, isRoot);
    }

}