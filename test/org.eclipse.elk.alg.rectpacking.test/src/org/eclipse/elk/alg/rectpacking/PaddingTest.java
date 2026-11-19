package org.eclipse.elk.alg.rectpacking;
/*******************************************************************************
 * Copyright (c) 2020 Kiel University and others.
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


import static org.junit.Assert.assertEquals;

import org.eclipse.elk.alg.rectpacking.options.RectPackingOptions;
import org.eclipse.elk.alg.test.PlainJavaInitialization;
import org.eclipse.elk.core.math.ElkPadding;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.util.BasicProgressMonitor;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.util.ElkGraphUtil;
import org.junit.BeforeClass;
import org.junit.Test;

/**
 * Test padding for rectpacking graphs.
 * 
 * Padding is difficult for rectpacking, since it should influence the aspect ratio.
 *
 */
public class PaddingTest {

    @BeforeClass
    public static void init() {
        PlainJavaInitialization.initializePlainJavaLayout();
    }
    
    @Test
    public void testTopPadding() {
        ElkNode parent = ElkGraphUtil.createGraph();
        ElkNode n1 = ElkGraphUtil.createNode(parent);
        ElkNode n2 = ElkGraphUtil.createNode(parent);
        ElkNode n3 = ElkGraphUtil.createNode(parent);
        ElkNode n4 = ElkGraphUtil.createNode(parent);
        parent.setProperty(CoreOptions.ALGORITHM, RectPackingOptions.ALGORITHM_ID);
        parent.setProperty(RectPackingOptions.ASPECT_RATIO, 1.3);
        parent.setProperty(CoreOptions.SPACING_NODE_NODE, 10.0);
        parent.setProperty(CoreOptions.PADDING, new ElkPadding(1000, 0, 0, 0));
        n1.setDimensions(30, 30);
        n2.setDimensions(30, 30);
        n3.setDimensions(30, 30);
        n4.setDimensions(30, 30);
        //
        //
        //
        // __  __  __  __
        //|__||__||__||__|
        
        RectPackingLayoutProvider layoutProvider = new RectPackingLayoutProvider();
        layoutProvider.layout(parent, new BasicProgressMonitor());
        assertEquals("", 150.0, parent.getWidth(), 1);
        assertEquals("", 1030.0, parent.getHeight(), 1);
        assertEquals("", 0.0, n1.getX(), 1);
        assertEquals("", 1000.0, n1.getY(), 1);
        assertEquals("", 40.0, n2.getX(), 1);
        assertEquals("", 1000.0, n2.getY(), 1);
        assertEquals("", 80.0, n3.getX(), 1);
        assertEquals("", 1000.0, n3.getY(), 1);
        assertEquals("", 120.0, n4.getX(), 1);
        assertEquals("", 1000.0, n4.getY(), 1);      
    }
    
    @Test
    public void testLeftPadding() {
        ElkNode parent = ElkGraphUtil.createGraph();
        ElkNode n1 = ElkGraphUtil.createNode(parent);
        ElkNode n2 = ElkGraphUtil.createNode(parent);
        ElkNode n3 = ElkGraphUtil.createNode(parent);
        ElkNode n4 = ElkGraphUtil.createNode(parent);
        parent.setProperty(CoreOptions.ALGORITHM, RectPackingOptions.ALGORITHM_ID);
        parent.setProperty(RectPackingOptions.ASPECT_RATIO, 1.3);
        parent.setProperty(CoreOptions.SPACING_NODE_NODE, 10.0);
        parent.setProperty(CoreOptions.PADDING, new ElkPadding(0, 0, 0, 1000));
        n1.setDimensions(30, 30);
        n2.setDimensions(30, 30);
        n3.setDimensions(30, 30);
        n4.setDimensions(30, 30);
        //                    __
        //                   |__|
        //                   |__| 
        //                   |__|
        //                   |__|
        
        RectPackingLayoutProvider layoutProvider = new RectPackingLayoutProvider();
        layoutProvider.layout(parent, new BasicProgressMonitor());
        assertEquals("", 1030.0, parent.getWidth(), 1);
        assertEquals("", 150.0, parent.getHeight(), 1);
        assertEquals("", 1000.0, n1.getX(), 1);
        assertEquals("", 0.0, n1.getY(), 1);
        assertEquals("", 1000.0, n2.getX(), 1);
        assertEquals("", 40.0, n2.getY(), 1);
        assertEquals("", 1000.0, n3.getX(), 1);
        assertEquals("", 80.0, n3.getY(), 1);
        assertEquals("", 1000.0, n4.getX(), 1);
        assertEquals("", 120.0, n4.getY(), 1);      
    }
    
    @Test
    public void testBottomPadding() {
        ElkNode parent = ElkGraphUtil.createGraph();
        ElkNode n1 = ElkGraphUtil.createNode(parent);
        ElkNode n2 = ElkGraphUtil.createNode(parent);
        ElkNode n3 = ElkGraphUtil.createNode(parent);
        ElkNode n4 = ElkGraphUtil.createNode(parent);
        parent.setProperty(CoreOptions.ALGORITHM, RectPackingOptions.ALGORITHM_ID);
        parent.setProperty(RectPackingOptions.ASPECT_RATIO, 1.3);
        parent.setProperty(CoreOptions.SPACING_NODE_NODE, 10.0);
        parent.setProperty(CoreOptions.PADDING, new ElkPadding(0, 0, 1000, 0));
        n1.setDimensions(30, 30);
        n2.setDimensions(30, 30);
        n3.setDimensions(30, 30);
        n4.setDimensions(30, 30);
        // __  __  __  __
        //|__||__||__||__|
        //
        //
        //
        //
        
        RectPackingLayoutProvider layoutProvider = new RectPackingLayoutProvider();
        layoutProvider.layout(parent, new BasicProgressMonitor());
        assertEquals("", 150.0, parent.getWidth(), 1);
        assertEquals("", 1030.0, parent.getHeight(), 1);
        assertEquals("", 0.0, n1.getX(), 1);
        assertEquals("", 0.0, n1.getY(), 1);
        assertEquals("", 40.0, n2.getX(), 1);
        assertEquals("", 0.0, n2.getY(), 1);
        assertEquals("", 80.0, n3.getX(), 1);
        assertEquals("", 0.0, n3.getY(), 1);
        assertEquals("", 120.0, n4.getX(), 1);
        assertEquals("", 0.0, n4.getY(), 1);      
    }
    
    @Test
    public void testRightPadding() {
        ElkNode parent = ElkGraphUtil.createGraph();
        ElkNode n1 = ElkGraphUtil.createNode(parent);
        ElkNode n2 = ElkGraphUtil.createNode(parent);
        ElkNode n3 = ElkGraphUtil.createNode(parent);
        ElkNode n4 = ElkGraphUtil.createNode(parent);
        parent.setProperty(CoreOptions.ALGORITHM, RectPackingOptions.ALGORITHM_ID);
        parent.setProperty(RectPackingOptions.ASPECT_RATIO, 1.3);
        parent.setProperty(CoreOptions.SPACING_NODE_NODE, 10.0);
        parent.setProperty(CoreOptions.PADDING, new ElkPadding(0, 1000, 0, 0));
        n1.setDimensions(30, 30);
        n2.setDimensions(30, 30);
        n3.setDimensions(30, 30);
        n4.setDimensions(30, 30);
        // __
        //|__|
        //|__| 
        //|__|
        //|__|
        
        RectPackingLayoutProvider layoutProvider = new RectPackingLayoutProvider();
        layoutProvider.layout(parent, new BasicProgressMonitor());
        assertEquals("", 1030.0, parent.getWidth(), 1);
        assertEquals("", 150.0, parent.getHeight(), 1);
        assertEquals("", 0.0, n1.getX(), 1);
        assertEquals("", 0.0, n1.getY(), 1);
        assertEquals("", 0.0, n2.getX(), 1);
        assertEquals("", 40.0, n2.getY(), 1);
        assertEquals("", 0.0, n3.getX(), 1);
        assertEquals("", 80.0, n3.getY(), 1);
        assertEquals("", 0.0, n4.getX(), 1);
        assertEquals("", 120.0, n4.getY(), 1);      
    }

}
