/*******************************************************************************
 * Copyright (c) 2024 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 * 
 * SPDX-License-Identifier: EPL-2.0 
 *******************************************************************************/
package org.eclipse.elk.alg.mrtree.test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;

import java.util.Arrays;
import java.util.Collection;

import org.eclipse.elk.alg.mrtree.options.MrTreeOptions;
import org.eclipse.elk.alg.test.PlainJavaInitialization;
import org.eclipse.elk.core.LayoutConfigurator;
import org.eclipse.elk.core.RecursiveGraphLayoutEngine;
import org.eclipse.elk.core.UnsupportedGraphException;
import org.eclipse.elk.core.data.LayoutAlgorithmResolver;
import org.eclipse.elk.core.math.ElkPadding;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.util.BasicProgressMonitor;
import org.eclipse.elk.core.util.ElkUtil;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.util.ElkGraphUtil;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameters;

/**
 * Tests that the graph size calculated by MrTree is correct.
 *
 */
@RunWith(Parameterized.class)
public class MrTreeGraphSizeTest {

    private static double DOUBLE_EQ_EPSILON = 10e-5;
    
    
    @Parameters
    public static Collection<Object[]> testValues() {
        return Arrays.asList(new Object[][] {
            {20, 20, 10, 10, 10, 10, 10}, {15, 30, 7, 9, 11, 13, 15}, {25, 15, 0, 0, 0, 0, 20}
        });
    }
    
    private final double nodeWidth;
    private final double nodeHeight;
    private final double paddingLeft;
    private final double paddingRight;
    private final double paddingTop;
    private final double paddingBottom;
    private final double nodeNodeSpacing;
    
    public MrTreeGraphSizeTest(double nodeWidth, double nodeHeight, double paddingLeft, double paddingRight, double paddingTop, double paddingBottom, double nodeNodeSpacing) {
        this.nodeWidth = nodeWidth;
        this.nodeHeight = nodeHeight;
        this.paddingLeft = paddingLeft;
        this.paddingRight = paddingRight;
        this.paddingTop = paddingTop;
        this.paddingBottom = paddingBottom;
        this.nodeNodeSpacing = nodeNodeSpacing;
    }
    
    /**
     * Tests that the graph size is correct.
     */
    @Test
    public void GraphSizeCalculationTest() {
        PlainJavaInitialization.initializePlainJavaLayout();
        ElkNode graph = ElkGraphUtil.createGraph();
        graph.setProperty(CoreOptions.ALGORITHM, MrTreeOptions.ALGORITHM_ID);
        graph.setProperty(CoreOptions.PADDING, new ElkPadding(paddingTop, paddingRight, paddingBottom, paddingLeft));
        graph.setProperty(CoreOptions.SPACING_NODE_NODE, nodeNodeSpacing);
        
        ElkNode n1 = ElkGraphUtil.createNode(graph);
        n1.setDimensions(nodeWidth, nodeHeight);
        ElkNode n2 = ElkGraphUtil.createNode(graph);
        n2.setDimensions(nodeWidth, nodeHeight);
        ElkNode n3 = ElkGraphUtil.createNode(graph);
        n3.setDimensions(nodeWidth, nodeHeight);
        
        ElkGraphUtil.createSimpleEdge(n1, n2);
        ElkGraphUtil.createSimpleEdge(n1, n3);
        
        // prepare layout engine
        LayoutConfigurator config = new LayoutConfigurator();
        ElkUtil.applyVisitors(graph, config, new LayoutAlgorithmResolver());
        // call layout with layout engine
        try {
            new RecursiveGraphLayoutEngine().layout(graph, new BasicProgressMonitor());
        } catch (UnsupportedGraphException exception) {
            fail(exception.toString());
        }
        
        assertEquals(paddingLeft + nodeWidth + nodeNodeSpacing + nodeWidth + paddingRight, graph.getWidth(), DOUBLE_EQ_EPSILON);
        assertEquals(paddingTop + nodeHeight + nodeNodeSpacing + nodeHeight + paddingBottom, graph.getHeight(), DOUBLE_EQ_EPSILON);
    }
}
