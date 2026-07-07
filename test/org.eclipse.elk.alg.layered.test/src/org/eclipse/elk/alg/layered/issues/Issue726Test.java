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

import java.util.List;

import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.intermediate.greedyswitch.TestGraphCreator;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.layered.p5edges.BasicEdgeRouterTest;
import org.eclipse.elk.alg.layered.p5edges.OrthogonalEdgeRouter;
import org.eclipse.elk.alg.test.framework.LayoutTestRunner;
import org.eclipse.elk.alg.test.framework.annotations.Algorithm;
import org.eclipse.elk.alg.test.framework.annotations.DefaultConfiguration;
import org.eclipse.elk.alg.test.framework.annotations.GraphResourceProvider;
import org.eclipse.elk.alg.test.framework.annotations.TestAfterProcessor;
import org.eclipse.elk.alg.test.framework.io.AbstractResourcePath;
import org.eclipse.elk.alg.test.framework.io.ModelResourcePath;
import org.eclipse.elk.core.math.KVectorChain;
import org.junit.runner.RunWith;

import com.google.common.collect.Lists;
import com.google.common.collect.Streams;

@RunWith(LayoutTestRunner.class)
@Algorithm(LayeredOptions.ALGORITHM_ID)
@DefaultConfiguration()
public class Issue726Test extends TestGraphCreator {

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Sources

    @GraphResourceProvider
    public List<AbstractResourcePath> testGraphs() {
        return Lists.newArrayList(
                new ModelResourcePath("tickets/layered/726_edgesMissingBendPoints.elkt"));
    }
    

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Tests
    
    @TestAfterProcessor(OrthogonalEdgeRouter.class)
    public void testOrthogonalEdges(final Object graph) {
        // Inspect each edge
        ((LGraph) graph).getLayers().stream()
            .flatMap(layer -> layer.getNodes().stream())
            .flatMap(node -> Streams.stream(node.getOutgoingEdges()))
            .forEach(edge -> testOrthogonalEdge(edge));
    }
    
    private void testOrthogonalEdge(final LEdge edge) {
        // Assemble edge route
        KVectorChain route = new KVectorChain();
        route.add(edge.getSource().getAbsoluteAnchor());
        route.addAll(edge.getBendPoints());
        route.add(edge.getTarget().getAbsoluteAnchor());
        
        BasicEdgeRouterTest.checkEdgeRouteIsOrthogonal(edge, route);
    }
    
}
