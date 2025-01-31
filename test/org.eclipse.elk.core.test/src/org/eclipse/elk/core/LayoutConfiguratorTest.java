/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.core;

import static org.junit.Assert.assertEquals;

import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.util.ElkUtil;
import org.eclipse.elk.graph.ElkGraphElement;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.util.ElkGraphUtil;
import org.junit.Test;

/**
 * Test the {@link LayoutConfigurator} class with regard to proper configuration of a passed graph. 
 */
public class LayoutConfiguratorTest {

    private static final double EPSILON = 0.0001d;
    
    @Test
    public void testElementOverridesClass() {

        Graph g = new Graph();
        
        LayoutConfigurator lc = new LayoutConfigurator();
        lc.configure(ElkNode.class).setProperty(CoreOptions.SPACING_NODE_NODE, 22d);
        lc.configure(g.n1).setProperty(CoreOptions.SPACING_NODE_NODE, 12d);
        lc.configure(g.n2).setProperty(CoreOptions.SPACING_NODE_NODE, 32d);

        ElkUtil.applyVisitors(g.root, lc);
        
        assertEquals(22d, g.root.getProperty(CoreOptions.SPACING_NODE_NODE), EPSILON);
        assertEquals(12d, g.n1.getProperty(CoreOptions.SPACING_NODE_NODE), EPSILON);
        assertEquals(32d, g.n2.getProperty(CoreOptions.SPACING_NODE_NODE), EPSILON);
    }

    @Test
    public void testCombineElementAndClass() {

        Graph g = new Graph();
        
        LayoutConfigurator lc = new LayoutConfigurator();
        lc.configure(ElkNode.class).setProperty(CoreOptions.SPACING_NODE_NODE, 42d);
        lc.configure(g.n1).setProperty(CoreOptions.SPACING_EDGE_NODE, 42d);

        ElkUtil.applyVisitors(g.root, lc);
        
        assertEquals(42d, g.n1.getProperty(CoreOptions.SPACING_NODE_NODE), EPSILON);
        assertEquals(42d, g.n1.getProperty(CoreOptions.SPACING_EDGE_NODE), EPSILON);
    }
    
    @Test
    public void testCombineMultipleClass() {

        Graph g = new Graph();
        
        LayoutConfigurator lc = new LayoutConfigurator();
        lc.configure(ElkNode.class).setProperty(CoreOptions.SPACING_NODE_NODE, 42d);
        lc.configure(ElkGraphElement.class).setProperty(CoreOptions.SPACING_EDGE_NODE, 42d);

        ElkUtil.applyVisitors(g.root, lc);
        
        assertEquals(42d, g.n1.getProperty(CoreOptions.SPACING_NODE_NODE), EPSILON);
        assertEquals(42d, g.n1.getProperty(CoreOptions.SPACING_EDGE_NODE), EPSILON);
    }
    
    @Test
    public void testOverridingExistingOptions() {
        
        Graph g = new Graph();
        g.root.setProperty(CoreOptions.SPACING_NODE_NODE, 13d);
        
        LayoutConfigurator lc = new LayoutConfigurator();
        lc.configure(ElkNode.class).setProperty(CoreOptions.SPACING_NODE_NODE, 42d);
        
        ElkUtil.applyVisitors(g.root, lc);
        
        assertEquals(42d, g.root.getProperty(CoreOptions.SPACING_NODE_NODE), EPSILON);
    }
    
    @Test
    public void testPreventOverridingExistingOptions() {
        
        Graph g = new Graph();
        g.root.setProperty(CoreOptions.SPACING_NODE_NODE, 13d);
        
        LayoutConfigurator lc = new LayoutConfigurator();
        lc.configure(ElkNode.class).setProperty(CoreOptions.SPACING_NODE_NODE, 42d);
        lc.addFilter(LayoutConfigurator.NO_OVERWRITE);
        
        ElkUtil.applyVisitors(g.root, lc);
        
        assertEquals(13d, g.root.getProperty(CoreOptions.SPACING_NODE_NODE), EPSILON);
    }
    
    private class Graph {
        ElkNode root;
        ElkNode n1;
        ElkNode n2;
        public Graph() {
            root = ElkGraphUtil.createGraph();
            n1 = ElkGraphUtil.createNode(root);
            n2 = ElkGraphUtil.createNode(root);
        }
    }
    
}
