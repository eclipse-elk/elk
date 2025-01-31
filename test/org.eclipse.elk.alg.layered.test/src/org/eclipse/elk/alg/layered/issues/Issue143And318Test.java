/*******************************************************************************
 * Copyright (c) 2020 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.alg.layered.issues;

import static org.junit.Assert.assertNotEquals;

import java.awt.geom.Line2D;
import java.util.ArrayList;
import java.util.List;

import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.test.framework.LayoutTestRunner;
import org.eclipse.elk.alg.test.framework.annotations.Algorithm;
import org.eclipse.elk.alg.test.framework.annotations.DefaultConfiguration;
import org.eclipse.elk.alg.test.framework.annotations.GraphResourceProvider;
import org.eclipse.elk.alg.test.framework.io.AbstractResourcePath;
import org.eclipse.elk.alg.test.framework.io.ModelResourcePath;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.math.KVectorChain;
import org.eclipse.elk.core.util.ElkUtil;
import org.eclipse.elk.graph.ElkEdge;
import org.eclipse.elk.graph.ElkEdgeSection;
import org.eclipse.elk.graph.ElkNode;
import org.junit.Test;
import org.junit.runner.RunWith;

import com.google.common.collect.Lists;

/**
 * Test for issues 143 and 318. Both refer to similar problems, and thus 
 */
@RunWith(LayoutTestRunner.class)
@Algorithm(LayeredOptions.ALGORITHM_ID)
@DefaultConfiguration
public class Issue143And318Test {

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Sources

    @GraphResourceProvider
    public List<AbstractResourcePath> testGraphs() {
        return Lists.newArrayList(
                new ModelResourcePath("tickets/layered/143_smallVerticalEdgeSegmentSpacing.elkt"),
                new ModelResourcePath("tickets/layered/318_overlappingHorizontalEdgeSegments.elkt"));
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Tests
    
    /**
     * The closes that two horizontal edge segments may come to each other. This is hard-coded, but works for the two
     * test models. If this test is to be extended, this must be more complex since the thresholds are computed for
     * each layer separately.
     */
    private static final double PROXIMITY_THRESHOLD = 2.0;

    @Test
    public void testForCloseEdgeProximities(final ElkNode graph) {
        List<Line2D.Double> horizontalSegments = gatherHorizontalEdgeSegments(graph);
        
        // Pairwise comparison of all segments
        for (int segIndex1 = 0; segIndex1 < horizontalSegments.size(); segIndex1++) {
            Line2D.Double seg1 = horizontalSegments.get(segIndex1);
            
            for (int segIndex2 = segIndex1 + 1; segIndex2 < horizontalSegments.size(); segIndex2++) {
                Line2D.Double seg2 = horizontalSegments.get(segIndex2);
                
                if (intersectHorizontally(seg1, seg2)) {
                    assertNotEquals(seg1.y1, seg2.y1, PROXIMITY_THRESHOLD);
                }
            }
        }
    }

    private static final double TOLERANCE = 0.01;
    
    private List<Line2D.Double> gatherHorizontalEdgeSegments(final ElkNode graph) {
        List<Line2D.Double> horizontalSegments = new ArrayList<>();
        
        for (ElkEdge edge : graph.getContainedEdges()) {
            for (ElkEdgeSection edgeSection : edge.getSections()) {
                KVectorChain bendPoints = ElkUtil.createVectorChain(edgeSection);
                
                KVector prevPoint = null;
                KVector currPoint = null;
                
                for (KVector bendPoint : bendPoints) {
                    // Advance to new point
                    prevPoint = currPoint;
                    currPoint = bendPoint;
                    
                    if (prevPoint != null) {
                        // Check if the current segment is horizontal
                        if (Math.abs(prevPoint.y - currPoint.y) < TOLERANCE) {
                            // Be sure to always go from smaller x values to larger x values
                            if (prevPoint.x < currPoint.x) {
                                horizontalSegments.add(
                                        new Line2D.Double(prevPoint.x, prevPoint.y, currPoint.x, currPoint.y));
                                
                            } else {
                                horizontalSegments.add(
                                        new Line2D.Double(currPoint.x, currPoint.y, prevPoint.x, prevPoint.y));
                            }
                        }
                    }
                }
            }
        }
        
        return horizontalSegments;
    }
    
    private boolean intersectHorizontally(Line2D.Double seg1, Line2D.Double seg2) {
        assert seg1.x1 < seg1.x2;
        assert seg2.x1 < seg2.x2;
        
        return seg1.x1 < seg2.x2 && seg2.x1 < seg1.x2;
    }

}
