/*******************************************************************************
 * Copyright (c) 2020 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 * 
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later  *******************************************************************************/
package org.eclipse.elk.core.util;

import static org.junit.Assert.assertEquals;

import java.util.EnumSet;

import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.options.ContentAlignment;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.util.ElkGraphUtil;
import org.junit.Test;

/**
 * Class for testing methods in ElkUtil
 */
public class ElkUtilTest {

    @Test
    public void translateWithContentAlignmentTest() {
        // Imagine a graph with one node that has inner behavior.
        // If this node is resized, the contents of said node should be alignment
        // based on the ContentAlignment specific by the parent of the node.
        
        // Test no alignment set
        {
            ElkNode parent = createContentAlignmentTestGraph(null);
            ElkNode node = parent.getChildren().get(0);
            ElkNode innerBehavior = node.getChildren().get(0);
            ElkUtil.translate(parent.getChildren().get(0), new KVector(120, 120), new KVector(100, 100));
            assertEquals(0 , innerBehavior.getX(), 1);
            assertEquals(0 , innerBehavior.getY(), 1);
        }
        // Test top left
        {
            ElkNode parent = createContentAlignmentTestGraph(ContentAlignment.topLeft());
            ElkNode node = parent.getChildren().get(0);
            ElkNode innerBehavior = node.getChildren().get(0);
            ElkUtil.translate(parent.getChildren().get(0), new KVector(120, 120), new KVector(100, 100));
            assertEquals(0 , innerBehavior.getX(), 1);
            assertEquals(0 , innerBehavior.getY(), 1);
        }
        // Test top center
        {
            ElkNode parent = createContentAlignmentTestGraph(ContentAlignment.topCenter());
            ElkNode node = parent.getChildren().get(0);
            ElkNode innerBehavior = node.getChildren().get(0);
            ElkUtil.translate(parent.getChildren().get(0), new KVector(120, 120), new KVector(100, 100));
            assertEquals(10 , innerBehavior.getX(), 1);
            assertEquals(0 , innerBehavior.getY(), 1);
        }
        // Test top right
        {
            ElkNode parent = createContentAlignmentTestGraph(EnumSet.of(ContentAlignment.V_TOP, ContentAlignment.H_RIGHT));
            ElkNode node = parent.getChildren().get(0);
            ElkNode innerBehavior = node.getChildren().get(0);
            ElkUtil.translate(parent.getChildren().get(0), new KVector(120, 120), new KVector(100, 100));
            assertEquals(20 , innerBehavior.getX(), 1);
            assertEquals(0 , innerBehavior.getY(), 1);
        }
        // Test center left
        {
            ElkNode parent = createContentAlignmentTestGraph(EnumSet.of(ContentAlignment.V_CENTER, ContentAlignment.H_LEFT));
            ElkNode node = parent.getChildren().get(0);
            ElkNode innerBehavior = node.getChildren().get(0);
            ElkUtil.translate(parent.getChildren().get(0), new KVector(120, 120), new KVector(100, 100));
            assertEquals(0 , innerBehavior.getX(), 1);
            assertEquals(10 , innerBehavior.getY(), 1);
        }
        // Test center center
        {
            ElkNode parent = createContentAlignmentTestGraph(ContentAlignment.centerCenter());
            ElkNode node = parent.getChildren().get(0);
            ElkNode innerBehavior = node.getChildren().get(0);
            ElkUtil.translate(parent.getChildren().get(0), new KVector(120, 120), new KVector(100, 100));
            assertEquals(10 , innerBehavior.getX(), 1);
            assertEquals(10 , innerBehavior.getY(), 1);
        }
        // Test center right
        {
            ElkNode parent = createContentAlignmentTestGraph(EnumSet.of(ContentAlignment.V_CENTER, ContentAlignment.H_RIGHT));
            ElkNode node = parent.getChildren().get(0);
            ElkNode innerBehavior = node.getChildren().get(0);
            ElkUtil.translate(parent.getChildren().get(0), new KVector(120, 120), new KVector(100, 100));
            assertEquals(20 , innerBehavior.getX(), 1);
            assertEquals(10 , innerBehavior.getY(), 1);
        }
        // Test bottom left
        {
            ElkNode parent = createContentAlignmentTestGraph(EnumSet.of(ContentAlignment.V_BOTTOM, ContentAlignment.H_LEFT));
            ElkNode node = parent.getChildren().get(0);
            ElkNode innerBehavior = node.getChildren().get(0);
            ElkUtil.translate(parent.getChildren().get(0), new KVector(120, 120), new KVector(100, 100));
            assertEquals(0 , innerBehavior.getX(), 1);
            assertEquals(20 , innerBehavior.getY(), 1);
        }
        // Test bottom center
        {
            ElkNode parent = createContentAlignmentTestGraph(EnumSet.of(ContentAlignment.V_BOTTOM, ContentAlignment.H_CENTER));
            ElkNode node = parent.getChildren().get(0);
            ElkNode innerBehavior = node.getChildren().get(0);
            ElkUtil.translate(parent.getChildren().get(0), new KVector(120, 120), new KVector(100, 100));
            assertEquals(10 , innerBehavior.getX(), 1);
            assertEquals(20 , innerBehavior.getY(), 1);
        }
        // Test bottom right
        {
            ElkNode parent = createContentAlignmentTestGraph(ContentAlignment.bottomRight());
            ElkNode node = parent.getChildren().get(0);
            ElkNode innerBehavior = node.getChildren().get(0);
            ElkUtil.translate(parent.getChildren().get(0), new KVector(120, 120), new KVector(100, 100));
            assertEquals(20 , innerBehavior.getX(), 1);
            assertEquals(20 , innerBehavior.getY(), 1);
        }
    }
    
    /**
     * Creates a graph with one node with inner behavior.
     * @param contentAlignment The content alignment of the graph.
     * @return The graph with one node with dimension (100, 100) and inner behavior of dimension (80, 80)
     */
    public ElkNode createContentAlignmentTestGraph(final EnumSet<ContentAlignment> contentAlignment) {
        ElkNode parent = ElkGraphUtil.createGraph();
        ElkNode node = ElkGraphUtil.createNode(parent);
        node.setProperty(CoreOptions.CONTENT_ALIGNMENT, contentAlignment);
        ElkNode innerBehavior = ElkGraphUtil.createNode(node);
        node.setDimensions(100, 100);
        innerBehavior.setDimensions(80, 80);
        return parent;
    }
}
