/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
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

package org.eclipse.elk.alg.radial.p2routing;

import org.eclipse.elk.alg.radial.InternalProperties;
import org.eclipse.elk.alg.radial.RadialLayoutPhases;
import org.eclipse.elk.core.alg.ILayoutPhase;
import org.eclipse.elk.core.alg.LayoutProcessorConfiguration;
import org.eclipse.elk.core.math.ElkMath;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.util.IElkProgressMonitor;
import org.eclipse.elk.graph.ElkEdge;
import org.eclipse.elk.graph.ElkEdgeSection;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.ElkPort;
import org.eclipse.elk.graph.util.ElkGraphUtil;

/**
 * An edge router which draws a straight line between nodes. The line points from center to center.
 */
public class StraightLineEdgeRouter implements ILayoutPhase<RadialLayoutPhases, ElkNode> {

    @Override
    public void process(final ElkNode graph, final IElkProgressMonitor progressMonitor) {
        progressMonitor.begin("Straight Line Edge Routing", 1);
        progressMonitor.logGraph(graph, "Before");
        ElkNode root = graph.getProperty(InternalProperties.ROOT_NODE);
        routeEdges(root);
        progressMonitor.logGraph(graph, "After");
    }

    /**
     * Route edges from node center to node center. Then clip it, to not cross the node.
     */
    public void routeEdges(final ElkNode node) {
        for (ElkEdge edge : ElkGraphUtil.allOutgoingEdges(node)) {
            if (!(edge.getSources().get(0) instanceof ElkPort)) {
                ElkNode target = ElkGraphUtil.connectableShapeToNode(edge.getTargets().get(0));
                if (!edge.isHierarchical()) {
                    double sourceX = node.getX() + node.getWidth() / 2;
                    double sourceY = node.getY() + node.getHeight() / 2;

                    double targetX = target.getX() + target.getWidth() / 2;
                    double targetY = target.getY() + target.getHeight() / 2;

                    // Clipping
                    KVector vector = new KVector();
                    vector.x = targetX - sourceX;
                    vector.y = targetY - sourceY;
                    KVector sourceClip = new KVector(vector.x, vector.y);
                    ElkMath.clipVector(sourceClip, node.getWidth(), node.getHeight());
                    vector.x -= sourceClip.x;
                    vector.y -= sourceClip.y;

                    sourceX = targetX - vector.x;
                    sourceY = targetY - vector.y;

                    KVector targetClip = new KVector(vector.x, vector.y);
                    ElkMath.clipVector(targetClip, target.getWidth(), target.getHeight());
                    vector.x -= targetClip.x;
                    vector.y -= targetClip.y;

                    targetX = sourceX + vector.x;
                    targetY = sourceY + vector.y;

                    ElkEdgeSection section = ElkGraphUtil.firstEdgeSection(edge, true, true);
                    section.setStartLocation(sourceX, sourceY);
                    section.setEndLocation(targetX, targetY);
                    routeEdges(target);
                }
            }
        }
    }

    @Override
    public LayoutProcessorConfiguration<RadialLayoutPhases, ElkNode> getLayoutProcessorConfiguration(
            final ElkNode graph) {
        return null;
    }

}
