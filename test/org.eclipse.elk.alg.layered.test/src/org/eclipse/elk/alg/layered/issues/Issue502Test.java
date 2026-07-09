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
package org.eclipse.elk.alg.layered.issues;

import static org.junit.Assert.assertTrue;

import java.util.List;

import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.test.framework.LayoutTestRunner;
import org.eclipse.elk.alg.test.framework.annotations.Algorithm;
import org.eclipse.elk.alg.test.framework.annotations.DefaultConfiguration;
import org.eclipse.elk.alg.test.framework.annotations.GraphResourceProvider;
import org.eclipse.elk.alg.test.framework.io.AbstractResourcePath;
import org.eclipse.elk.alg.test.framework.io.ModelResourcePath;
import org.eclipse.elk.core.math.ElkRectangle;
import org.eclipse.elk.graph.ElkBendPoint;
import org.eclipse.elk.graph.ElkEdge;
import org.eclipse.elk.graph.ElkEdgeSection;
import org.eclipse.elk.graph.ElkNode;
import org.junit.Test;
import org.junit.runner.RunWith;

import com.google.common.collect.Lists;

/**
 * Test for issue 502.
 */
@RunWith(LayoutTestRunner.class)
@Algorithm(LayeredOptions.ALGORITHM_ID)
@DefaultConfiguration
public class Issue502Test {
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Sources

    @GraphResourceProvider
    public List<AbstractResourcePath> testGraphs() {
        return Lists.newArrayList(new ModelResourcePath("tickets/layered/502_collapsingCompoundNode.elkt"));
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Tests
    
    @Test
    public void testLargeEnoughNodes(final ElkNode graph) {
        graph.getChildren().stream()
            .forEach(node -> testLargeEnoughNode(node));
    }

    private void testLargeEnoughNode(final ElkNode node) {
        ElkRectangle requiredSize = new ElkRectangle();
        
        for (ElkNode child : node.getChildren()) {
            ElkRectangle childRect = new ElkRectangle(child.getX(), child.getY(), child.getWidth(), child.getHeight());
            requiredSize.union(childRect);
        }
        
        for (ElkEdge edge : node.getContainedEdges()) {
            for (ElkEdgeSection section : edge.getSections()) {
                for (ElkBendPoint bendPoint : section.getBendPoints()) {
                    requiredSize.width = Math.max(requiredSize.width, bendPoint.getX());
                    requiredSize.height = Math.max(requiredSize.height, bendPoint.getY());
                }
            }
        }
        
        assertTrue("Node not large enough for content. Size required: (" + requiredSize.width + ", "
                + requiredSize.height + "). " + "Actual size: (" + node.getWidth() + ", " + node.getHeight() + ")",
                node.getWidth() >= requiredSize.width && node.getHeight() >= requiredSize.height);
    }

}
