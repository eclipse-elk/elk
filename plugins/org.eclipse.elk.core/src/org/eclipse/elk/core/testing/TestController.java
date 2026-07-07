/*******************************************************************************
 * Copyright (c) 2018, 2019 Kiel University and others.
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
package org.eclipse.elk.core.testing;

import java.util.LinkedHashSet;
import java.util.Set;

import org.eclipse.elk.core.AbstractLayoutProvider;
import org.eclipse.elk.core.alg.ILayoutProcessor;
import org.eclipse.elk.core.data.LayoutAlgorithmData;

/**
 * A test controller is used to provide listeners with a way to track the execution of a layout algorithm structured
 * into {@link ILayoutProcessor}s. The main reason for the existence of this class is that our test framework supports
 * white box tests, that is, tests that examine the internal state of a layout algorithm before or after certain
 * processors are run. Thus, the test framework uses test controllers to be notified whenever a layout processor is
 * about to be executed or has just finished executing.
 * 
 * <p>A test controller needs to be registered with the layout provider by calling the
 * {@link #install(AbstractLayoutProvider)} method. Once the layout algorithm has finished executing, it is a very good
 * idea to 
 */
public class TestController {

    /** Identifier of the layout algorithm that should be tested with this controller. */
    private final String layoutAlgorithmId;
    /** The layout provider the controller is currently installed on, if any. */
    private IWhiteBoxTestable currentLayoutProvider;
    /** Listeners that want to be notified during the execution of a layout algorithm. */
    private Set<ILayoutExecutionListener> listeners = new LinkedHashSet<>();

    /**
     * Creates a new instance for testing the layout algorithm with the given ID. Each controller is created for a
     * specific layout algorithm.
     */
    public TestController(final String layoutAlgorithmId) {
        this.layoutAlgorithmId = layoutAlgorithmId;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Setup

    /**
     * Returns the identifier of the layout algorithm that should be tested with this test controller.
     */
    public String getTargetAlgorithmId() {
        return layoutAlgorithmId;
    }

    /**
     * Checks whether this test controller targets layout algorithms described by the given data object. If so, this
     * test controller can be installed on the algorithm by calling {@link #install(AbstractLayoutProvider)}.
     */
    public boolean targets(final LayoutAlgorithmData algorithmData) {
        return algorithmData.getId().equals(layoutAlgorithmId);
    }

    /**
     * Installs this test controller on the given layout provider. Before calling this method, clients should have made
     * sure that the two fit by calling {@link #targets(LayoutAlgorithmData)}.
     * 
     * @throws IllegalArgumentException
     *             of the layout provider does not implement {@link IWhiteBoxTestable}.
     * @throws IllegalStateException
     *             if the test controller is currently installed on another layout provider.
     */
    public void install(final AbstractLayoutProvider layoutProvider) {
        if (currentLayoutProvider != null) {
            throw new IllegalStateException("Test controller may be installed on only one layout provider at a time");
        } else if (layoutProvider instanceof IWhiteBoxTestable) {
            currentLayoutProvider = (IWhiteBoxTestable) layoutProvider;
            currentLayoutProvider.setTestController(this);
        } else {
            throw new IllegalArgumentException(
                    "Test controllers can only be installed on white-box testable layout algorithms");
        }
    }
    
    /**
     * Uninstalls this test controller from its current layout provider, if any. If it is not currently installed on
     * one, this method does nothing.
     */
    public void uninstall() {
        if (currentLayoutProvider != null) {
            currentLayoutProvider.setTestController(null);
            currentLayoutProvider = null;
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Notification Methods

    /**
     * Notifies all listeners that a processor is about to be executed on a non-root graph.
     * 
     * @param graph
     *            the graph.
     * @param processor
     *            the processor.
     */
    public void notifyProcessorReady(final Object graph, final ILayoutProcessor<?> processor) {
        for (ILayoutExecutionListener listener : listeners) {
            listener.layoutProcessorReady(processor, graph, false);
        }
    }

    /**
     * Notifies all listeners that a processor just finished executing on a non-root graph.
     * 
     * @param graph
     *            the graph.
     * @param processor
     *            the processor.
     */
    public void notifyProcessorFinished(final Object graph, final ILayoutProcessor<?> processor) {
        for (ILayoutExecutionListener listener : listeners) {
            listener.layoutProcessorFinished(processor, graph, false);
        }
    }

    /**
     * Notifies all listeners that a processor is about to be executed on a root graph.
     * 
     * @param graph
     *            the graph.
     * @param processor
     *            the processor.
     */
    public void notifyRootProcessorReady(final Object graph, final ILayoutProcessor<?> processor) {
        for (ILayoutExecutionListener listener : listeners) {
            listener.layoutProcessorReady(processor, graph, true);
        }
    }

    /**
     * Notifies all listeners that a processor just finished executing on a root graph.
     * 
     * @param graph
     *            the graph.
     * @param processor
     *            the processor.
     */
    public void notifyRootProcessorFinished(final Object graph, final ILayoutProcessor<?> processor) {
        for (ILayoutExecutionListener listener : listeners) {
            listener.layoutProcessorFinished(processor, graph, true);
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Listener Management

    /**
     * Adds a listener to the test controller.
     * 
     * @param listener
     *            the listener to be added.
     */
    public void addLayoutExecutionListener(final ILayoutExecutionListener listener) {
        listeners.add(listener);
    }

    /**
     * Removes a listener from the test controller.
     * 
     * @param listener
     *            the listener to be added.
     */
    public void removeLayoutExecutionListener(final ILayoutExecutionListener listener) {
        listeners.add(listener);
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Listener Interfaces

    /**
     * Interface implemented by listeners interested in the execution of a layout algorithm.
     */
    public interface ILayoutExecutionListener {

        /**
         * Fired when a layout processor is about to be executed.
         * 
         * @param processor
         *            the layout processor about to be executed.
         * @param graph
         *            the graph the processor will be run on.
         * @param isRoot
         *            {@code true} if the graph is the root of the graph hierarchy laid out by the layout algorithm.
         *            This does not have to be the root of the entire graph!
         */
        void layoutProcessorReady(ILayoutProcessor<?> processor, Object graph, boolean isRoot);

        /**
         * Fired when a layout processor finished executing.
         * 
         * @param processor
         *            the layout processor about to be executed.
         * @param graph
         *            the graph the processor will be run on.
         * @param isRoot
         *            {@code true} if the graph is the root of the graph hierarchy laid out by the layout algorithm.
         *            This does not have to be the root of the entire graph!
         */
        void layoutProcessorFinished(ILayoutProcessor<?> processor, Object graph, boolean isRoot);

    }

}
