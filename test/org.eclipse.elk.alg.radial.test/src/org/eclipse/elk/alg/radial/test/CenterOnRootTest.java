/*******************************************************************************
 * Copyright (c) 2023 Kiel University and others.
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
package org.eclipse.elk.alg.radial.test;

import static org.junit.Assert.assertEquals;

import org.eclipse.elk.alg.radial.RadialLayoutProvider;
import org.eclipse.elk.alg.radial.options.RadialOptions;
import org.eclipse.elk.alg.test.PlainJavaInitialization;
import org.eclipse.elk.core.math.ElkMargin;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.util.BasicProgressMonitor;
import org.eclipse.elk.graph.ElkEdge;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.util.ElkGraphUtil;
import org.junit.BeforeClass;
import org.junit.Test;

/**
 * Test for the center on root option.
 *
 */
public class CenterOnRootTest {
    
    @BeforeClass
    public static void init() {
        PlainJavaInitialization.initializePlainJavaLayout();
    }
    
    /**
     * Layout a tree of height 1 and center it. The center point of the root
     * node should be placed in the middle of the graph as defined by its width
     * and height.
     */
    @Test
    public void testSimpleCentering() {
        ElkNode parent = ElkGraphUtil.createGraph();
        ElkNode root = ElkGraphUtil.createNode(parent);
        ElkNode n1 = ElkGraphUtil.createNode(parent);
        ElkEdge e1 = ElkGraphUtil.createSimpleEdge(root, n1);
        ElkNode n2 = ElkGraphUtil.createNode(parent);
        ElkEdge e2 = ElkGraphUtil.createSimpleEdge(root, n2);
        ElkNode n3 = ElkGraphUtil.createNode(parent);
        ElkEdge e3 = ElkGraphUtil.createSimpleEdge(root, n3);
        
        parent.setProperty(CoreOptions.ALGORITHM, RadialOptions.ALGORITHM_ID);
        parent.setProperty(RadialOptions.CENTER_ON_ROOT, true);
        
        RadialLayoutProvider layoutProvider = new RadialLayoutProvider();
        layoutProvider.layout(parent, new BasicProgressMonitor());
        
        ElkMargin margins = root.getProperty(CoreOptions.MARGINS);
        assertEquals("Horizontal centering", parent.getWidth()/2, root.getX() + margins.left + root.getWidth()/2, 0.1);
        assertEquals("Vertical centering", parent.getHeight()/2, root.getY() + margins.top + root.getHeight()/2, 0.1);
    }
    
    /**
     * Layout of a tree with a height larger than 1.
     */
    @Test
    public void testLargerGraphCentering() {
        ElkNode parent = ElkGraphUtil.createGraph();
        ElkNode root = ElkGraphUtil.createNode(parent);
        
        ElkNode n1 = ElkGraphUtil.createNode(parent);
        ElkEdge e1 = ElkGraphUtil.createSimpleEdge(root, n1);
        ElkNode n2 = ElkGraphUtil.createNode(parent);
        ElkEdge e2 = ElkGraphUtil.createSimpleEdge(root, n2);
        ElkNode n3 = ElkGraphUtil.createNode(parent);
        ElkEdge e3 = ElkGraphUtil.createSimpleEdge(root, n3);
        
        ElkNode n11 = ElkGraphUtil.createNode(parent);
        ElkEdge e11 = ElkGraphUtil.createSimpleEdge(n1, n11);
        ElkNode n12 = ElkGraphUtil.createNode(parent);
        ElkEdge e12 = ElkGraphUtil.createSimpleEdge(n1, n12);
        ElkNode n13 = ElkGraphUtil.createNode(parent);
        ElkEdge e13 = ElkGraphUtil.createSimpleEdge(n1, n13);
        
        parent.setProperty(CoreOptions.ALGORITHM, RadialOptions.ALGORITHM_ID);
        parent.setProperty(RadialOptions.CENTER_ON_ROOT, true);
        
        RadialLayoutProvider layoutProvider = new RadialLayoutProvider();
        layoutProvider.layout(parent, new BasicProgressMonitor());
        
        ElkMargin margins = root.getProperty(CoreOptions.MARGINS);
        assertEquals("Horizontal centering", parent.getWidth()/2, root.getX() + margins.left + root.getWidth()/2, 0.1);
        assertEquals("Vertical centering", parent.getHeight()/2, root.getY() + margins.top + root.getHeight()/2, 0.1);
    }

}
