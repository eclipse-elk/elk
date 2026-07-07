/*******************************************************************************
 * Copyright (c) 2023 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.layered.issues;

import static org.junit.Assert.assertEquals;

import java.util.List;

import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.test.framework.LayoutTestRunner;
import org.eclipse.elk.alg.test.framework.annotations.Algorithm;
import org.eclipse.elk.alg.test.framework.annotations.DefaultConfiguration;
import org.eclipse.elk.alg.test.framework.annotations.GraphResourceProvider;
import org.eclipse.elk.alg.test.framework.io.AbstractResourcePath;
import org.eclipse.elk.alg.test.framework.io.ModelResourcePath;
import org.eclipse.elk.core.options.EdgeLabelPlacement;
import org.eclipse.elk.graph.ElkEdge;
import org.eclipse.elk.graph.ElkLabel;
import org.eclipse.elk.graph.ElkNode;
import org.junit.Test;
import org.junit.runner.RunWith;

import com.google.common.collect.Lists;

/**
 * Test for issue 905.
 */
@RunWith(LayoutTestRunner.class)
@Algorithm(LayeredOptions.ALGORITHM_ID)
@DefaultConfiguration
public class Issue905Test {

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Sources

    @GraphResourceProvider
    public List<AbstractResourcePath> testGraphs() {
        return Lists.newArrayList(
                // 2 identical graph except the location of the labels of the edge
                new ModelResourcePath("tickets/layered/905_beginLabelAboveEdge.elkt"),
                new ModelResourcePath("tickets/layered/905_beginLableBellowEdge.elkt"));
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Tests

    @Test
    public void testEdgeLabelAlignment(final ElkNode graph) {
        assertEquals("Wrong input graph: The graph must contained one node.", 1, graph.getChildren().size());
        assertEquals("Wrong input graph: The node of the graph must contained one edge.", 1, graph.getChildren().get(0).getContainedEdges().size());
        testBegingLabel(graph.getChildren().get(0).getContainedEdges().get(0));
    }
    
    /**
     * Check that the labels of the edge is at the same location after layout, independently of the location before the layout. Only tail label is problematic with issue 905.
     * @param edge The edge having labels to check
     */
    private void testBegingLabel(final ElkEdge edge) {
        double expectedXLocationForTailLabel = 10;
        double expectedYLocationForTailLabel = 92;
        
        double expectedXLocationForCenterLabel = 61;
        double expectedYLocationForCenterLabel = 91;
        
        double expectedXLocationForHeadLabel = 124;
        double expectedYLocationForHeadLabel = 92;

        
        for (ElkLabel label : edge.getLabels()) {
            EdgeLabelPlacement currentLabelPlacement = label.getProperty(LayeredOptions.EDGE_LABELS_PLACEMENT);
            if (EdgeLabelPlacement.TAIL.equals(currentLabelPlacement)) {
                assertEquals("The x location of the tail label is wrong after layout.", expectedXLocationForTailLabel, label.getX(), 0);
                assertEquals("The y location of the tail label is wrong after layout.", expectedYLocationForTailLabel, label.getY(), 0);
            } else if (EdgeLabelPlacement.CENTER.equals(currentLabelPlacement)) {
                assertEquals("The x location of the center label is wrong after layout.", expectedXLocationForCenterLabel, label.getX(), 0);
                assertEquals("The y location of the center label is wrong after layout.", expectedYLocationForCenterLabel, label.getY(), 0);
            } else if (EdgeLabelPlacement.HEAD.equals(currentLabelPlacement)) {
                assertEquals("The x location of the head label is wrong after layout.", expectedXLocationForHeadLabel, label.getX(), 0);
                assertEquals("The y location of the head label is wrong after layout.", expectedYLocationForHeadLabel, label.getY(), 0);
            }
        }
    }


}
