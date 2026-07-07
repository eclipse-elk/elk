/*******************************************************************************
 * Copyright (c) 2020 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.layered.issues;

import static org.junit.Assert.assertFalse;

import java.util.ArrayList;
import java.util.Deque;
import java.util.LinkedList;
import java.util.List;

import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.test.framework.LayoutTestRunner;
import org.eclipse.elk.alg.test.framework.annotations.Algorithm;
import org.eclipse.elk.alg.test.framework.annotations.Configurator;
import org.eclipse.elk.alg.test.framework.annotations.DefaultConfiguration;
import org.eclipse.elk.alg.test.framework.annotations.GraphResourceProvider;
import org.eclipse.elk.alg.test.framework.io.AbstractResourcePath;
import org.eclipse.elk.alg.test.framework.io.FileNameFilter;
import org.eclipse.elk.alg.test.framework.io.ModelResourcePath;
import org.eclipse.elk.core.math.ElkRectangle;
import org.eclipse.elk.core.options.NodeLabelPlacement;
import org.eclipse.elk.graph.ElkNode;
import org.junit.Test;
import org.junit.runner.RunWith;

import com.google.common.collect.Lists;

/**
 * Test for issue 603.
 */
@RunWith(LayoutTestRunner.class)
@Algorithm(LayeredOptions.ALGORITHM_ID)
@DefaultConfiguration
public class Issue603Test {

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Sources

    @GraphResourceProvider
    public List<AbstractResourcePath> testGraphs() {
        return Lists.newArrayList(
                new ModelResourcePath("tickets/layered/**").withFilter(new FileNameFilter("603.+\\.elkt")));
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Configuration

    @Configurator
    public void configureNodeLabelPlacement(final ElkNode graph) {
        graph.getChildren().stream()
            .forEach(node -> node.setProperty(
                    LayeredOptions.NODE_LABELS_PLACEMENT,
                    NodeLabelPlacement.insideTopCenter()));
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Tests

    @Test
    public void testNoLabelNodeOverlaps(final ElkNode graph) {
        Deque<ElkNode> nodeQueue = new LinkedList<>(graph.getChildren());
        while (!nodeQueue.isEmpty()) {
            ElkNode node = nodeQueue.poll();
            
            // We're looking for overlaps between the node's labels and its children
            if (!node.getLabels().isEmpty() && !node.getChildren().isEmpty()) {
                doTestNoLabelNodeOverlaps(node);
            }
            
            nodeQueue.addAll(node.getChildren());
        }
    }

    private void doTestNoLabelNodeOverlaps(ElkNode node) {
        List<ElkRectangle> nodeAndLabelRects = new ArrayList<>();
        
        // Convert labels and nodes to rectangles
        node.getLabels().stream()
            .map(l -> new ElkRectangle(l.getX(), l.getY(), l.getWidth(), l.getHeight()))
            .forEach(rect -> nodeAndLabelRects.add(rect));
        
        node.getChildren().stream()
            .map(n -> new ElkRectangle(n.getX(), n.getY(), n.getWidth(), n.getHeight()))
            .forEach(rect -> nodeAndLabelRects.add(rect));
        
        // Check for overlaps
        for (int first = 0; first < nodeAndLabelRects.size(); first++) {
            for (int second = first + 1; second < nodeAndLabelRects.size(); second++) {
                assertFalse(nodeAndLabelRects.get(first).intersects(nodeAndLabelRects.get(second)));
            }
        }
    }

}
