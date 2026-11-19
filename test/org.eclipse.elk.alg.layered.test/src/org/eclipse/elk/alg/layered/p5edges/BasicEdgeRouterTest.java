/*******************************************************************************
 * Copyright (c) 2013, 2020 Kiel University and others.
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
package org.eclipse.elk.alg.layered.p5edges;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;

import java.util.Deque;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import org.eclipse.elk.alg.layered.intermediate.greedyswitch.TestGraphCreator;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.test.framework.LayoutTestRunner;
import org.eclipse.elk.alg.test.framework.annotations.Algorithm;
import org.eclipse.elk.alg.test.framework.annotations.ConfiguratorProvider;
import org.eclipse.elk.alg.test.framework.annotations.DefaultConfiguration;
import org.eclipse.elk.alg.test.framework.annotations.GraphResourceProvider;
import org.eclipse.elk.alg.test.framework.io.AbstractResourcePath;
import org.eclipse.elk.alg.test.framework.io.FileExtensionFilter;
import org.eclipse.elk.alg.test.framework.io.ModelResourcePath;
import org.eclipse.elk.core.LayoutConfigurator;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.math.KVectorChain;
import org.eclipse.elk.core.options.EdgeRouting;
import org.eclipse.elk.core.util.ElkUtil;
import org.eclipse.elk.graph.ElkEdge;
import org.eclipse.elk.graph.ElkEdgeSection;
import org.eclipse.elk.graph.ElkNode;
import org.junit.Test;
import org.junit.runner.RunWith;

import com.google.common.collect.Lists;

@RunWith(LayoutTestRunner.class)
@Algorithm(LayeredOptions.ALGORITHM_ID)
@DefaultConfiguration()
public class BasicEdgeRouterTest extends TestGraphCreator {

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Sources

    @GraphResourceProvider
    public List<AbstractResourcePath> testGraphs() {
        return Lists.newArrayList(
                new ModelResourcePath("realworld/ptolemy/**/").withFilter(new FileExtensionFilter("elkg")));
    }
    

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Configuration
    
    @ConfiguratorProvider
    public LayoutConfigurator orthogonalConfigurator() {
        return configuratorFor(EdgeRouting.ORTHOGONAL);
    }
    
//    @ConfiguratorProvider
//    public LayoutConfigurator polylineConfigurator() {
//        return configuratorFor(EdgeRouting.POLYLINE);
//    }
//    
//    @ConfiguratorProvider
//    public LayoutConfigurator splineConfigurator() {
//        return configuratorFor(EdgeRouting.SPLINES);
//    }
    
    private LayoutConfigurator configuratorFor(final EdgeRouting strategy) {
        LayoutConfigurator config = new LayoutConfigurator();
        config.configure(ElkNode.class).setProperty(LayeredOptions.EDGE_ROUTING, strategy);
        return config;
    }
    

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Tests
    
    @Test
    public void testOrthogonalEdges(final ElkNode graph) {
        // Iterate over all of the graph's edges, through the whole hierarchy
        Deque<ElkNode> nodeQueue = new LinkedList<>();
        nodeQueue.add(graph);
        
        while (!nodeQueue.isEmpty()) {
            ElkNode node = nodeQueue.poll();
            nodeQueue.addAll(node.getChildren());
            
            for (ElkEdge edge : node.getContainedEdges()) {
                for (ElkEdgeSection section : edge.getSections()) {
                    checkEdgeRouteIsOrthogonal(edge, ElkUtil.createVectorChain(section));
                }
            }
            
        }
    }
    
    private static final double TOLERANCE = 0.05;

    /**
     * Ensures that the given list of points describe an orthogonal route through the diagram. Utility method
     * that can be used by other tests as well.
     */
    public static void checkEdgeRouteIsOrthogonal(final Object edge, final KVectorChain route) {
        // The route should contain at least two points
        if (route.size() < 2) {
            fail(edge + ": Route with less than 2 points");
        }
        
        boolean horizontal = isHorizontal(route.get(0), route.get(1));
        Iterator<KVector> pointIter = route.iterator();
        
        KVector prev = pointIter.next();
        while (pointIter.hasNext()) {
            KVector curr = pointIter.next();
            
            String error = edge + ": horizontal = " + horizontal + " prev = " + prev + " curr = " + curr;
            
            if (horizontal) {
                assertEquals(error, prev.y, curr.y, TOLERANCE);
            } else {
                assertEquals(error, prev.x, curr.x, TOLERANCE);
            }
            
            horizontal = !horizontal;
            
            prev = curr;
        }
    }
    
    /**
     * Checks whether the line from p1 to p2 is horizontal.
     */
    private static boolean isHorizontal(final KVector p1, final KVector p2) {
        return Math.abs(p1.y - p2.y) < TOLERANCE;
    }
    
}
