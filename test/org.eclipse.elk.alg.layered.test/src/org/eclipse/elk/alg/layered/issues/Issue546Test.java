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

import java.awt.geom.Rectangle2D;
import java.util.LinkedList;
import java.util.List;

import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.test.framework.LayoutTestRunner;
import org.eclipse.elk.alg.test.framework.annotations.Algorithm;
import org.eclipse.elk.alg.test.framework.annotations.DefaultConfiguration;
import org.eclipse.elk.alg.test.framework.annotations.GraphResourceProvider;
import org.eclipse.elk.alg.test.framework.io.AbstractResourcePath;
import org.eclipse.elk.alg.test.framework.io.ModelResourcePath;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.util.ElkUtil;
import org.eclipse.elk.graph.ElkEdge;
import org.eclipse.elk.graph.ElkEdgeSection;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.ElkPort;
import org.eclipse.elk.graph.util.ElkGraphUtil;
import org.junit.Test;
import org.junit.runner.RunWith;

import com.google.common.collect.Lists;

/**
 * Test for issue 546.
 */
@RunWith(LayoutTestRunner.class)
@Algorithm(LayeredOptions.ALGORITHM_ID)
@DefaultConfiguration
public class Issue546Test {
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Sources

    @GraphResourceProvider
    public List<AbstractResourcePath> testGraphs() {
        return Lists.newArrayList(new ModelResourcePath("tickets/layered/546_borderGaps.elkt"));
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Tests
    
    private static final double TOLERANCE = 0.5;
    
    @Test
    public void testEdgeStartsAtPort(final ElkNode graph) {
        // We'll look at all edges
        LinkedList<ElkNode> nodeQueue = new LinkedList<>();
        nodeQueue.add(graph);
        
        while (!nodeQueue.isEmpty()) {
            ElkNode currNode = nodeQueue.poll();
            nodeQueue.addAll(currNode.getChildren());
            
            // The following code assumes that each edge has a single source and edge section, and that each
            // edge is connected to ports
            for (ElkEdge edge : currNode.getContainedEdges()) {
                // Retrieve absolute port rectangle, slightly enlarged
                ElkPort srcPort = ElkGraphUtil.connectableShapeToPort(edge.getSources().get(0));
                KVector srcPortPos = ElkUtil.absolutePosition(srcPort);
                Rectangle2D.Double srcPortRect = new Rectangle2D.Double(
                        srcPortPos.x - TOLERANCE,
                        srcPortPos.y - TOLERANCE,
                        srcPort.getWidth() + 2 * TOLERANCE,
                        srcPort.getHeight() + 2 * TOLERANCE);
                
                // Retrieve absolut edge source position
                ElkEdgeSection section = edge.getSections().get(0);
                KVector srcPoint = new KVector(section.getStartX(), section.getStartY());
                ElkUtil.toAbsolute(srcPoint, currNode);
                
                // Ensure that the source point is roughly on the port anchor
                assertTrue("Source point " + srcPoint + " is not close enough to port " + srcPortRect,
                        srcPortRect.contains(srcPoint.x, srcPoint.y));
            }
        }
    }

}
