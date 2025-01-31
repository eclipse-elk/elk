/*******************************************************************************
 * Copyright (c) 2018 TypeFox GmbH and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.core;

import static org.junit.Assert.*;

import org.eclipse.elk.alg.test.PlainJavaInitialization;
import org.eclipse.elk.core.data.LayoutAlgorithmData;
import org.eclipse.elk.core.data.LayoutMetaDataService;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.util.BasicProgressMonitor;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.util.ElkGraphUtil;
import org.junit.BeforeClass;
import org.junit.Test;

/**
 * Tests for {@link RecursiveGraphLayoutEngine}.
 */
public class RecursiveGraphLayoutEngineTest {
    
    @BeforeClass
    public static void initPlainJavaLayout() {
        PlainJavaInitialization.initializePlainJavaLayout();
    }
    
    @Test
    public void testUnresolvedGraph() {
        Graph graph = new Graph();
        graph.root.setProperty(CoreOptions.ALGORITHM, "org.eclipse.elk.box");
        RecursiveGraphLayoutEngine engine = new RecursiveGraphLayoutEngine();
        engine.layout(graph.root, new BasicProgressMonitor());
        
        assertTrue(graph.root.getWidth() > 0);
        assertTrue(graph.root.getHeight() > 0);
    }
    
    @Test
    public void testResolvedGraph() {
        Graph graph = new Graph();
        LayoutAlgorithmData algorithmData = LayoutMetaDataService.getInstance().getAlgorithmData("org.eclipse.elk.box");
        graph.root.setProperty(CoreOptions.RESOLVED_ALGORITHM, algorithmData);
        RecursiveGraphLayoutEngine engine = new RecursiveGraphLayoutEngine();
        engine.layout(graph.root, new BasicProgressMonitor());
        
        assertTrue(graph.root.getWidth() > 0);
        assertTrue(graph.root.getHeight() > 0);
    }
    
    @Test(expected = UnsupportedConfigurationException.class)
    public void testUnknownAlgorithmId() {
        Graph graph = new Graph();
        graph.root.setProperty(CoreOptions.ALGORITHM, "foo.Bar");
        RecursiveGraphLayoutEngine engine = new RecursiveGraphLayoutEngine();
        engine.layout(graph.root, new BasicProgressMonitor());
    }
    
    @Test
    public void testEmptyAlgorithmId() {
        Graph graph = new Graph();
        RecursiveGraphLayoutEngine engine = new RecursiveGraphLayoutEngine();
        engine.layout(graph.root, new BasicProgressMonitor());
        
        assertEquals("org.eclipse.elk.layered", graph.root.getProperty(CoreOptions.RESOLVED_ALGORITHM).getId());
    }
    
    private class Graph {
        ElkNode root;
        private ElkNode n1;
        private ElkNode n2;
        public Graph() {
            root = ElkGraphUtil.createGraph();
            n1 = ElkGraphUtil.createNode(root);
            n1.setDimensions(10, 10);
            n2 = ElkGraphUtil.createNode(root);
            n2.setDimensions(10, 10);
        }
    }
    
}
