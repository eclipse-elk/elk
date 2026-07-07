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
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Stream;

import org.junit.Ignore;
import org.junit.runner.Description;
import org.junit.runner.notification.RunNotifier;
import org.junit.runners.ParentRunner;
import org.junit.runners.model.FrameworkMethod;
import org.junit.runners.model.InitializationError;
import org.junit.runners.model.Statement;

/**
 * A test runner for a given experimental object. This test runner will perform layout on the experimental object and
 * call the appropriate test methods during and after the layout algorithm's execution. The actual experiment execution
 * works as follows:
 * 
 * <ol>
 * <li>{@link ParentRunner} calls {@link #childrenInvoker(RunNotifier)} to have us generate a {@link Statement} that
 * runs all of the tests. The statement is in charge of executing the layout tests.</li>
 * <li>A {@link ExperimentStatement} is then used to run the actual tests.</li>
 * </ol>
 * 
 * For more details, see the documentation of {@link ExperimentStatement}. These complications are mainly due to the
 * fact that our whitebox tests are not executed the way usual test methods are.
 */
final class ExperimentRunner extends ParentRunner<FrameworkMethod> {

    /** ID of the next {@link ExperimentRunner} to be created. */
    private static int nextRunnerId = 0;
    /** ID of this runner, used to disambiguate descriptions of test methods. */
    private int runnerId;

    /** The runner that created us. We'll use that to enumerate our tests. */
    private final LayoutTestRunner parentRunner;
    /** The experimental object we'll be performing tests on. */
    private final ExperimentalObject experimentalObject;
    /** Descriptions for our test methods. */
    private Map<FrameworkMethod, Description> testDescriptions = new HashMap<>();

    /**
     * Creates a new test runner that will perform tests on the given experimental object.
     * 
     * @throws InitializationError
     *             this shouldn't happen since any initialization problems should already have been caught by
     *             {@link LayoutTestRunner}.
     */
    ExperimentRunner(final LayoutTestRunner parentRunner, final ExperimentalObject experimentalObject)
            throws InitializationError {

        super(parentRunner.getTestClass().getJavaClass());

        this.runnerId = ++nextRunnerId;

        this.parentRunner = parentRunner;
        this.experimentalObject = experimentalObject;

        // Our parent makes all tests available to use. Now we just have to create descriptions for them
        Stream.concat(parentRunner.getBlackboxTests().stream(), parentRunner.getWhiteboxTests().stream())
                .forEach(test -> describeChild(test));
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // Getters

    /**
     * Returns the {@link LayoutTestRunner} that created us.
     */
    public LayoutTestRunner getParentRunner() {
        return parentRunner;
    }
    
    /**
     * Returns the experimental object for this runner.
     */
    public ExperimentalObject getExperimentalObject() {
        return experimentalObject;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // ParentRunner
    
    @Override
    protected String getName() {
        return experimentalObject.toString();
    }

    @Override
    protected boolean isIgnored(final FrameworkMethod child) {
        return child.getAnnotation(Ignore.class) != null;
    }

    @Override
    protected List<FrameworkMethod> getChildren() {
        return new ArrayList<>(testDescriptions.keySet());
    }

    @Override
    protected Description describeChild(final FrameworkMethod child) {
        if (!testDescriptions.containsKey(child)) {
            // Since other experiment runners will build descriptions for the same methods that we build descriptions
            // for, we need to disambiguate them.
            Description childDescription = Description.createTestDescription(getTestClass().getJavaClass(),
                    child.getName() + " (" + runnerId + ")", child.getAnnotations());
            testDescriptions.put(child, childDescription);
        }

        return testDescriptions.get(child);
    }

    @Override
    protected Statement childrenInvoker(final RunNotifier notifier) {
        // Instead of invoking our children directly, we use this statement to run tests
        return new ExperimentStatement(this, notifier);
    }

    @Override
    protected void runChild(final FrameworkMethod child, final RunNotifier notifier) {
        // We don't do anything here because this method is not called due to the way our implementation works
    }

    /**
     * Calls {@link #runLeaf(Statement, Description, RunNotifier)}, whose visibility cannot be changed to public.
     */
    void doRunLeaf(final Statement statement, final Description description, final RunNotifier notifier) {
        runLeaf(statement, description, notifier);
    }

}
