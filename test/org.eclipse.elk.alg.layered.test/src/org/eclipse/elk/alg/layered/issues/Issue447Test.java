/*******************************************************************************
 * Copyright (c) 2019 Kiel University and others.
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
import org.eclipse.elk.alg.test.GraphTestUtils;
import org.eclipse.elk.alg.test.framework.LayoutTestRunner;
import org.eclipse.elk.alg.test.framework.annotations.Algorithm;
import org.eclipse.elk.alg.test.framework.annotations.Configurator;
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
import org.junit.Test;
import org.junit.runner.RunWith;

import com.google.common.collect.Lists;

/**
 * Test for issue 447.
 */
@RunWith(LayoutTestRunner.class)
@Algorithm(LayeredOptions.ALGORITHM_ID)
@DefaultConfiguration
public class Issue447Test {
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Sources

    @GraphResourceProvider
    public List<AbstractResourcePath> testGraphs() {
        return Lists.newArrayList(new ModelResourcePath("tickets/layered/447_hierarchicalPortAnchors.elkt"));
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Configuration
    
    private static final double PORT_SIZE = 10;
    private static final double PORT_ANCHOR = 5;

    @Configurator
    public void configurePortSizesAndAnchors(final ElkNode graph) {
        GraphTestUtils.allNodes(graph, false).stream()
            .flatMap(node -> node.getPorts().stream())
            .forEach(this::configurePort);
    }
    
    private void configurePort(final ElkPort port) {
        port.setWidth(PORT_SIZE);
        port.setHeight(PORT_SIZE);
        port.setProperty(LayeredOptions.PORT_ANCHOR, new KVector(PORT_ANCHOR, PORT_ANCHOR));
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Tests
    
    /** Allowed delta when comparing coordinates. */
    private final double COORDINATE_FUZZYNESS = 0.5;

    /**
     * Check that all edges connect to the port anchors, not to the edges.
     */
    @Test
    public void testEdgesConnectToAnchors(final ElkNode graph) {
        GraphTestUtils.allNodes(graph, true).stream()
            .flatMap(node -> node.getContainedEdges().stream())
            .forEach(this::testEdge);
    }
    
    private void testEdge(final ElkEdge edge) {
        for (ElkEdgeSection section : edge.getSections()) {
            // A section may be incident to a port
            if (section.getIncomingShape() instanceof ElkPort) {
                testEdgeEndPoint(
                        edge.getContainingNode(),
                        (ElkPort) section.getIncomingShape(),
                        section.getStartX(),
                        section.getStartY());
            }
            
            if (section.getOutgoingShape() instanceof ElkPort) {
                testEdgeEndPoint(
                        edge.getContainingNode(),
                        (ElkPort) section.getOutgoingShape(),
                        section.getEndX(),
                        section.getEndY());
            }
        }
    }
    
    private void testEdgeEndPoint(final ElkNode containingNode, final ElkPort port, final double x, final double y) {
        // Port anchor in absolute coordinates
        KVector portAnchor = new KVector(port.getX() + PORT_ANCHOR, port.getY() + PORT_ANCHOR);
        KVector absolutePortAnchor = ElkUtil.toAbsolute(portAnchor, port.getParent());
        
        // Edge end point in absolute coordinates
        KVector endPoint = new KVector(x, y);
        KVector absoluteEndPoint = ElkUtil.toAbsolute(endPoint, containingNode);
        
        // Check if they match
        assertEquals("X coordinate of edge does not match X coordinate of port anchor",
                absolutePortAnchor.x, absoluteEndPoint.x, COORDINATE_FUZZYNESS);
        assertEquals("Y coordinate of edge does not match Y coordinate of port anchor",
                absolutePortAnchor.y, absoluteEndPoint.y, COORDINATE_FUZZYNESS);
    }

}
