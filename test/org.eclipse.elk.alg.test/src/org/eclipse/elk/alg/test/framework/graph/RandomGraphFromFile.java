/*******************************************************************************
 * Copyright (c) 2019, 2020 Kiel University and others.
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

import org.eclipse.elk.alg.test.framework.annotations.RandomGeneratorFile;
import org.eclipse.elk.alg.test.framework.io.AbsoluteResourcePath;
import org.eclipse.elk.alg.test.framework.io.AbstractResourcePath;
import org.eclipse.elk.alg.test.framework.io.FileExtensionFilter;
import org.eclipse.elk.alg.test.framework.util.TestUtil;
import org.eclipse.elk.core.debug.grandom.gRandom.Configuration;
import org.eclipse.elk.core.debug.grandom.gRandom.RandGraph;
import org.eclipse.elk.core.debug.grandom.generators.RandomGraphGenerator;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.emf.common.util.URI;
import org.eclipse.emf.ecore.EObject;
import org.eclipse.emf.ecore.resource.Resource;
import org.eclipse.emf.ecore.resource.ResourceSet;
import org.eclipse.emf.ecore.resource.impl.ResourceSetImpl;
import org.eclipse.emf.ecore.util.EcoreUtil;
import org.junit.runners.model.FrameworkMethod;
import org.junit.runners.model.TestClass;

/**
 * A test graph that generates input graphs using the {@link RandomGraphGenerator}. The configuration for the generator
 * is loaded from a file supplied by a configuration method. Instances of this class are generated whenever the
 * {@link RandomGeneratorFile} annotation is encountered.
 */
public final class RandomGraphFromFile extends TestGraph {
    
    /** A file filter that will accept standard ELK graph files. */
    public static final FileFilter RANDOM_GRAPH_FILE_FILTER = new FileExtensionFilter(".elkr");

    /** The file from which our random graph was created. */
    private final AbsoluteResourcePath resourcePath;
    /** Since a file can generate several random graphs, this is the one this instance will supply. */
    private final ElkNode graph;

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // Creation

    private RandomGraphFromFile(final AbsoluteResourcePath resourcePath, final ElkNode graph) {
        this.resourcePath = resourcePath;
        this.graph = graph;
    }

    /**
     * Loads random graphs from files as specified by configuration methods in the test class.
     * 
     * @param testClass
     *            the test class.
     * @param test
     *            instance of the test class to call methods.
     * @param errors
     *            list of error conditions encountered.
     * @return a list of graphs.
     */
    public static List<RandomGraphFromFile> fromTestClass(final TestClass testClass, final Object test,
            final List<Throwable> errors) {
        
        List<AbstractResourcePath> resourcePaths = new ArrayList<>();

        for (FrameworkMethod method : testClass.getAnnotatedMethods(RandomGeneratorFile.class)) {
            TestUtil.ensurePublic(method, errors);
            TestUtil.ensureReturnsType(method, errors, AbstractResourcePath.class);
            TestUtil.ensureParameters(method, errors);

            // Invoke the method to obtain the resource paths
            try {
                resourcePaths.add((AbstractResourcePath) method.invokeExplosively(test));
            } catch (Throwable e) {
                errors.add(e);
            }
        }
        
        // Set proper file filters and load graphs
        List<RandomGraphFromFile> graphs = new ArrayList<>();
        
        for (AbstractResourcePath abstractPath : resourcePaths) {
            abstractPath.setFilter(RANDOM_GRAPH_FILE_FILTER);
            
            for (AbsoluteResourcePath absolutePath : abstractPath.listResources()) {
                try {
                    for (ElkNode graph : randomGraphs(absolutePath)) {
                        graphs.add(new RandomGraphFromFile(absolutePath, graph));
                    }
                } catch (Throwable e) {
                    errors.add(new Exception("Unable to load random graphs from " + absolutePath));
                }
            }
        }
        
        return graphs;
    }
    
    /**
     * Creates random graphs from a random graph configuration file.
     */
    private static List<ElkNode> randomGraphs(final AbsoluteResourcePath path) throws Throwable {
        List<ElkNode> graphs = new ArrayList<>();

        // Load the random graph file
        ResourceSet resourceSet = new ResourceSetImpl();
        Resource resource = resourceSet.getResource(URI.createFileURI(path.getFile().getAbsolutePath()), true);

        if (resource != null) {
            resource.load(Collections.emptyMap());
            EObject eo = resource.getContents().get(0);
            if (eo instanceof RandGraph) {
                // Create random graphs
                RandGraph rand = (RandGraph) resource.getContents().get(0);

                RandomGraphGenerator graphGenerator = new RandomGraphGenerator(null);
                for (Configuration config : rand.getConfigs()) {
                    graphs.addAll(graphGenerator.generate(config));
                }
            }
        }

        return graphs;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // TestGraph
    
    @Override
    public ElkNode provideGraph(final Object testClass) {
        // Simply return a copy of the graph
        return EcoreUtil.copy(graph);
    }

    @Override
    public String toString() {
        return "randomGraphFile(" + resourcePath.getFile().getPath() + ")";
    }

}
