/*******************************************************************************
 * Copyright (c) 2026 Kiel University and others.
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

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.List;

import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.core.RecursiveGraphLayoutEngine;
import org.eclipse.elk.core.math.ElkPadding;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.util.NullElkProgressMonitor;
import org.eclipse.elk.graph.ElkEdge;
import org.eclipse.elk.graph.ElkEdgeSection;
import org.eclipse.elk.graph.ElkLabel;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.util.ElkGraphUtil;
import org.junit.Test;

/**
 * Verifies that inside self-loop labels move with their node (issue #1232).
 */
public class Issue1232Test {

    private static final double EPSILON = 0.000001;

    @Test
    public void testLabelMovesWithNodeWhenRootPaddingChanges() {
        ElkNode baseline = layout(12);
        ElkNode moved = layout(60);
        ElkEdge baselineEdge = baseline.getPorts().get(0).getOutgoingEdges().get(0);
        ElkEdge movedEdge = moved.getPorts().get(0).getOutgoingEdges().get(0);

        assertTrue("Root padding must move the containing node",
                moved.getX() > baseline.getX() && moved.getY() > baseline.getY());
        assertEquals(baseline.getWidth(), moved.getWidth(), EPSILON);
        assertEquals(baseline.getHeight(), moved.getHeight(), EPSILON);
        assertEquals("Moving the node must preserve the local edge route",
                route(baselineEdge, baseline), route(movedEdge, moved));

        ElkLabel baselineLabel = baselineEdge.getLabels().get(0);
        ElkLabel movedLabel = movedEdge.getLabels().get(0);
        assertEquals("Moving the node must preserve the local label x coordinate",
                baselineLabel.getX() - baseline.getX(), movedLabel.getX() - moved.getX(), EPSILON);
        assertEquals("Moving the node must preserve the local label y coordinate",
                baselineLabel.getY() - baseline.getY(), movedLabel.getY() - moved.getY(), EPSILON);
    }

    private ElkNode layout(final double rootPadding) {
        ElkNode root = ElkGraphUtil.createGraph();
        root.setProperty(CoreOptions.ALGORITHM, LayeredOptions.ALGORITHM_ID);
        root.setProperty(CoreOptions.PADDING, new ElkPadding(rootPadding));

        ElkNode node = ElkGraphUtil.createNode(root);
        node.setProperty(CoreOptions.INSIDE_SELF_LOOPS_ACTIVATE, true);

        ElkEdge edge = ElkGraphUtil.createSimpleEdge(ElkGraphUtil.createPort(node), ElkGraphUtil.createPort(node));
        edge.setProperty(CoreOptions.INSIDE_SELF_LOOPS_YO, true);
        ElkLabel label = ElkGraphUtil.createLabel(edge);
        label.setText("1");
        label.setDimensions(7, 12);

        new RecursiveGraphLayoutEngine().layout(root, new NullElkProgressMonitor());
        return node;
    }

    private List<Double> route(final ElkEdge edge, final ElkNode node) {
        assertEquals("The self loop must have one routed section", 1, edge.getSections().size());
        ElkEdgeSection section = edge.getSections().get(0);
        List<Double> points = new ArrayList<>();
        points.add(section.getStartX() - node.getX());
        points.add(section.getStartY() - node.getY());
        section.getBendPoints().forEach(point -> {
            points.add(point.getX() - node.getX());
            points.add(point.getY() - node.getY());
        });
        points.add(section.getEndX() - node.getX());
        points.add(section.getEndY() - node.getY());
        return points;
    }
}
