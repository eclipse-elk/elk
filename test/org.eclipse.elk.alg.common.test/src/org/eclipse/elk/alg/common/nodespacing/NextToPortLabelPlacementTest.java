/*******************************************************************************
 * Copyright (c) 2020 Kiel University and others.
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
package org.eclipse.elk.alg.common.nodespacing;

import static org.junit.Assert.assertTrue;

import java.util.Deque;
import java.util.LinkedList;
import java.util.List;

import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.test.framework.LayoutTestRunner;
import org.eclipse.elk.alg.test.framework.annotations.Algorithm;
import org.eclipse.elk.alg.test.framework.annotations.DefaultConfiguration;
import org.eclipse.elk.alg.test.framework.annotations.GraphResourceProvider;
import org.eclipse.elk.alg.test.framework.io.AbstractResourcePath;
import org.eclipse.elk.alg.test.framework.io.ModelResourcePath;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.options.PortLabelPlacement;
import org.eclipse.elk.graph.ElkEdge;
import org.eclipse.elk.graph.ElkLabel;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.ElkPort;
import org.eclipse.elk.graph.util.ElkGraphUtil;
import org.junit.Test;
import org.junit.runner.RunWith;

import com.google.common.collect.Lists;

/**
 * Tests whether nodes configured with {@link CoreOptions#PORT_LABELS_NEXT_TO_PORT_IF_POSSIBLE} actually have their port
 * labels placed next to their ports. If possible.
 */
@RunWith(LayoutTestRunner.class)
@Algorithm(LayeredOptions.ALGORITHM_ID)
@DefaultConfiguration(edges = false, nodes = false, ports = true)
public class NextToPortLabelPlacementTest {

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Sources

    @GraphResourceProvider
    public List<AbstractResourcePath> testGraphs() {
        return Lists.newArrayList(
                new ModelResourcePath("tests/core/label_placement/port_labels/next_to_port_if_possible_inside.elkt"),
                new ModelResourcePath("tests/core/label_placement/port_labels/next_to_port_if_possible_outside.elkt"));
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Tests
    
    @Test
    public void testNextToPortLabels(final ElkNode graph) {
        Deque<ElkNode> nodeQueue = new LinkedList<>(graph.getChildren());
        while (!nodeQueue.isEmpty()) {
            ElkNode currentNode = nodeQueue.poll();
            nodeQueue.addAll(currentNode.getChildren());
            
            boolean insideLabelPlacement =
                    currentNode.getProperty(CoreOptions.PORT_LABELS_PLACEMENT).contains(PortLabelPlacement.INSIDE);
            for (ElkPort port : currentNode.getPorts()) {
                testPortLabel(port, insideLabelPlacement);
            }
        }
    }

    private void testPortLabel(ElkPort port, final boolean insideLabelPlacement) {
        /* This code has a few assumptions:
         *   1. We're dealing with horizontal layouts.
         *   2. We're dealing with eastern and western ports.
         *   3. There's exactly one label per port.
         */
        ElkLabel label = port.getLabels().get(0);
        
        if (shouldLabelBePlacedNextToPort(port, insideLabelPlacement)) {
            assertTrue(label.getY() + label.getHeight() >= 0 && label.getY() <= port.getHeight());
        } else {
            // The label is to be placed above or below its port
            assertTrue(label.getY() + label.getHeight() <= 0 || label.getY() >= port.getHeight());
        }
    }
    
    private boolean shouldLabelBePlacedNextToPort(final ElkPort port, final boolean insideLabelPlacement) {
        // Shortcut: if there is no incident edge, we can place the label next to the port
        if (port.getIncomingEdges().isEmpty() && port.getOutgoingEdges().isEmpty()) {
            return true;
        }
        
        // Iterate over all edges and check whether there are (i) edges that connect to the parent node's insides, and
        // (ii) edges that connect to somewhere else
        boolean edgesToInsides = false;
        boolean edgesToSomewhereElse = false;
        
        for (ElkEdge outEdge : port.getOutgoingEdges()) {
            boolean insideEdge = ElkGraphUtil.isDescendant(
                    ElkGraphUtil.connectableShapeToNode(outEdge.getTargets().get(0)),
                    port.getParent());
            
            edgesToInsides |= insideEdge;
            edgesToSomewhereElse |= !insideEdge;
        }
        
        for (ElkEdge inEdge : port.getIncomingEdges()) {
            boolean insideEdge = ElkGraphUtil.isDescendant(
                    ElkGraphUtil.connectableShapeToNode(inEdge.getSources().get(0)),
                    port.getParent());
            
            edgesToInsides |= insideEdge;
            edgesToSomewhereElse |= !insideEdge;
        }
        
        return (insideLabelPlacement && !edgesToInsides) || (!insideLabelPlacement && !edgesToSomewhereElse);
    }

}
