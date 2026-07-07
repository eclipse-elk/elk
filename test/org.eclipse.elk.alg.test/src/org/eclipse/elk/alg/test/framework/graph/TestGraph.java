/*******************************************************************************
 * Copyright (c) 2019 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.test.framework.graph;

import org.eclipse.elk.graph.ElkNode;

/**
 * Test graphs provide input graphs for layout tests. How those graphs are created is up to the subclasses.
 */
public abstract class TestGraph {

    /**
     * Applies the strategy to load a graph.
     * 
     * @param test
     *            instance of the test class that runs test methods.
     * @return the test graph.
     * @throws Throwable
     *             if anything goes wrong while trying to provide the graph.
     */
    public abstract ElkNode provideGraph(Object test) throws Throwable;

    // Force sub classes to provide proper descriptions
    @Override
    public abstract String toString();

}
