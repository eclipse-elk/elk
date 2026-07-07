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
import java.util.HashSet;
import java.util.Set;

import org.eclipse.elk.alg.test.framework.annotations.FailIfNotExecuted;
import org.junit.runner.Description;
import org.junit.runner.Result;
import org.junit.runner.notification.Failure;
import org.junit.runner.notification.RunNotifier;
import org.junit.runner.notification.StoppedByUserException;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;

/**
 * A notifier which, instead of passing events on, simply records them.
 */
final class WhiteboxRecorder extends RunNotifier {

    /** Set of whitebox tests that have been executed. Important for {@link FailIfNotExecuted} support. */
    private final Set<Description> executedWhiteboxTests = new HashSet<>();
    /** Record of failures that occurred while running whitebox tests. */
    private final Multimap<Description, Failure> whiteboxTestFailures = HashMultimap.create();

    /** The test that is currently running. */
    private Description currentlyRunningTest = null;

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // Accessors

    /**
     * Checks whether the test with the given description was executed at all.
     */
    public boolean wasExecuted(final Description description) {
        return executedWhiteboxTests.contains(description);
    }

    /**
     * Returns the failures accumulated over the given test's test runs. The returned collection will never be
     * {@code null}.
     */
    public Collection<Failure> getFailures(final Description description) {
        return whiteboxTestFailures.get(description);
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // RunNotifier

    @Override
    public void fireTestRunStarted(final Description description) {
        // Do nothing
    }

    @Override
    public void fireTestRunFinished(final Result result) {
        // Do nothing
    }

    @Override
    public void fireTestStarted(final Description description) throws StoppedByUserException {
        currentlyRunningTest = description;
        executedWhiteboxTests.add(description);
    }

    @Override
    public void fireTestFailure(final Failure failure) {
        whiteboxTestFailures.put(currentlyRunningTest, failure);
    }

    @Override
    public void fireTestAssumptionFailed(final Failure failure) {
        whiteboxTestFailures.put(currentlyRunningTest, failure);
    }

    @Override
    public void fireTestIgnored(final Description description) {
        // Do nothing
    }

    @Override
    public void fireTestFinished(final Description description) {
        currentlyRunningTest = null;
    }

}
