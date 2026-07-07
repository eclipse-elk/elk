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

import java.io.FileFilter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import org.eclipse.elk.alg.test.framework.annotations.GraphResourceProvider;
import org.eclipse.elk.alg.test.framework.io.AbsoluteResourcePath;
import org.eclipse.elk.alg.test.framework.io.AbstractResourcePath;
import org.eclipse.elk.alg.test.framework.io.FileExtensionFilter;
import org.eclipse.elk.alg.test.framework.util.TestUtil;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.emf.common.util.URI;
import org.eclipse.emf.ecore.EObject;
import org.eclipse.emf.ecore.resource.Resource;
import org.eclipse.emf.ecore.resource.ResourceSet;
import org.eclipse.emf.ecore.resource.impl.ResourceSetImpl;
import org.junit.runners.model.FrameworkMethod;
import org.junit.runners.model.TestClass;

/**
 * A test graph that gets an input graph from a file supplied by a method. Instances of this class are generated
 * whenever the {@link GraphResourceProvider} annotation is encountered.
 */
public final class GraphFromFile extends TestGraph {

    /** A file filter that will accept standard ELK graph files. */
    public static final FileFilter GRAPH_FILE_FILTER = new FileExtensionFilter(".elkt", ".elkg");

    /** Path of the graph file to load. */
    private final AbsoluteResourcePath resourcePath;

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // Creation

    private GraphFromFile(final AbsoluteResourcePath resourcePath) {
        this.resourcePath = resourcePath;
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
    public static List<GraphFromFile> fromTestClass(final TestClass testClass, final Object test,
            final List<Throwable> errors) {

        List<AbstractResourcePath> resourcePaths = new ArrayList<>();

        for (FrameworkMethod method : testClass.getAnnotatedMethods(GraphResourceProvider.class)) {
            TestUtil.ensurePublic(method, errors);
            TestUtil.ensureReturnsType(method, errors, List.class);
            TestUtil.ensureParameters(method, errors);

            // Invoke the method to obtain the resource paths
            try {
                List<?> result = (List<?>) method.invokeExplosively(test);

                // Check that these are all resource paths
                if (!result.stream().allMatch(o -> o instanceof AbstractResourcePath)) {
                    errors.add(new Exception("Method " + method.getName() + " must return List<AbstractResourcePath"));
                } else {
                    result.stream().map(o -> (AbstractResourcePath) o).forEach(path -> resourcePaths.add(path));
                }
            } catch (Throwable e) {
                errors.add(e);
            }
        }

        // Set a default graph file filter in case no explicit filter has been specified
        resourcePaths.stream().filter(path -> path.getFilter() == null)
                .forEach(path -> path.setFilter(GRAPH_FILE_FILTER));

        // Turn the abstract paths into absolute paths
        return resourcePaths.stream().flatMap(path -> path.listResources().stream())
                .map(path -> new GraphFromFile(path)).collect(Collectors.toList());
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // TestGraph

    @Override
    public ElkNode provideGraph(final Object testClass) throws Throwable {
        // For this to work, we assume that PlainJavaInitialization was already triggered
        ResourceSet resourceSet = new ResourceSetImpl();
        Resource resource = resourceSet.getResource(URI.createFileURI(resourcePath.getFile().getAbsolutePath()), true);

        if (resource != null) {
            resource.load(Collections.emptyMap());
            EObject eo = resource.getContents().get(0);
            if (eo instanceof ElkNode) {
                return (ElkNode) eo;
            }
        }

        return null;
    }

    @Override
    public String toString() {
        return "graphFile(" + resourcePath.getFile().getPath() + ")";
    }

}
