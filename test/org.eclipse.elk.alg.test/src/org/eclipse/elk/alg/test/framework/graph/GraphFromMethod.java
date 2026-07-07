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

import java.util.ArrayList;
import java.util.List;

import org.eclipse.elk.alg.test.framework.annotations.GraphProvider;
import org.eclipse.elk.alg.test.framework.util.TestUtil;
import org.eclipse.elk.graph.ElkNode;
import org.junit.runners.model.FrameworkMethod;
import org.junit.runners.model.TestClass;

/**
 * A test graph that gets an input graph from a method. Instances of this class are generated whenever the
 * {@link GraphProvider} annotation is encountered.
 */
public final class GraphFromMethod extends TestGraph {
    
    /** The method we'll invoke to obtain our graph. */
    private final FrameworkMethod method;

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // Creation

    private GraphFromMethod(final FrameworkMethod method) {
        this.method = method;
    }

    /**
     * Loads graphs from files as specified by configuration methods in the test class.
     * 
     * @param testClass
     *            the test class.
     * @param test
     *            instance of the test class to call methods.
     * @param errors
     *            list of error conditions encountered.
     * @return a list of graphs.
     */
    public static List<GraphFromMethod> fromTestClass(final TestClass testClass, final Object test,
            final List<Throwable> errors) {
        
        List<GraphFromMethod> result = new ArrayList<>();

        for (FrameworkMethod method : testClass.getAnnotatedMethods(GraphProvider.class)) {
            TestUtil.ensurePublic(method, errors);
            TestUtil.ensureReturnsType(method, errors, ElkNode.class);
            TestUtil.ensureParameters(method, errors);

            result.add(new GraphFromMethod(method));
        }
        
        return result;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // TestGraph

    @Override
    public ElkNode provideGraph(final Object test) throws Throwable {
        return (ElkNode) method.invokeExplosively(test);
    }

    @Override
    public String toString() {
        return "graphMethod(" + method.getName() + ")";
    }

}
