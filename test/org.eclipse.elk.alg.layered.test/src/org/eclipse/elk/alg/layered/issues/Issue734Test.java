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

import static org.junit.Assert.assertEquals;

import java.util.DoubleSummaryStatistics;
import java.util.List;

import org.eclipse.elk.alg.layered.intermediate.greedyswitch.TestGraphCreator;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.test.framework.LayoutTestRunner;
import org.eclipse.elk.alg.test.framework.annotations.Algorithm;
import org.eclipse.elk.alg.test.framework.annotations.GraphResourceProvider;
import org.eclipse.elk.alg.test.framework.io.AbstractResourcePath;
import org.eclipse.elk.alg.test.framework.io.ModelResourcePath;
import org.eclipse.elk.core.util.ElkUtil;
import org.eclipse.elk.graph.ElkEdge;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.emf.common.util.EList;
import org.junit.Test;
import org.junit.runner.RunWith;

import com.google.common.collect.Lists;

@RunWith(LayoutTestRunner.class)
@Algorithm(LayeredOptions.ALGORITHM_ID)
public class Issue734Test extends TestGraphCreator {

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Sources

    @GraphResourceProvider
    public List<AbstractResourcePath> testGraphs() {
        return Lists.newArrayList(
                new ModelResourcePath("tickets/layered/734_nodeEdgeSpacingFlat.elkt"),
                new ModelResourcePath("tickets/layered/734_nodeEdgeSpacingFlatMoreComplex.elkt"),
                new ModelResourcePath("tickets/layered/734_nodeEdgeSpacingHierarchical.elkt"));
    }
    

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Tests
    
    @Test
    public void testOrthogonalEdges(final ElkNode graph) {
        // Compute the area occupied by nodes on the top level
        CoordinateSpan nodeSpan = computeNodeCoordinateSpan(graph.getChildren());
        CoordinateSpan edgeSpan = computeEdgeCoordinateSpan(graph.getContainedEdges());
        
        // If the edge is routed properly, the edge span should have the same distance to the node span on both
        // sides. In other words, their respective centres should be aligned
        double nodeSpanCenter = (nodeSpan.min + nodeSpan.max) / 2;
        double edgeSpanCenter = (edgeSpan.min + edgeSpan.max) / 2;
        
        assertEquals(nodeSpanCenter, edgeSpanCenter, 0.5);
    }

    private CoordinateSpan computeNodeCoordinateSpan(EList<ElkNode> nodes) {
        CoordinateSpan span = new CoordinateSpan();
        
        span.min = nodes.stream()
            .mapToDouble(node -> node.getX())
            .min()
            .orElse(0);
        span.max = nodes.stream()
            .mapToDouble(node -> node.getX() + node.getWidth())
            .max()
            .orElse(0);
        
        return span;
    }

    private CoordinateSpan computeEdgeCoordinateSpan(EList<ElkEdge> containedEdges) {
        DoubleSummaryStatistics stats = containedEdges.stream()
            .flatMap(edge -> edge.getSections().stream())
            .map(section -> ElkUtil.createVectorChain(section))
            .flatMap(chain -> chain.stream())
            .mapToDouble(vec -> vec.x)
            .summaryStatistics();
        return new CoordinateSpan(stats.getMin(), stats.getMax());
    }
    
    private static class CoordinateSpan {
        private double min = Double.POSITIVE_INFINITY;
        private double max = Double.NEGATIVE_INFINITY;
        
        public CoordinateSpan() {
        }

        public CoordinateSpan(double min, double max) {
            this.min = min;
            this.max = max;
        }
    }
    
}
