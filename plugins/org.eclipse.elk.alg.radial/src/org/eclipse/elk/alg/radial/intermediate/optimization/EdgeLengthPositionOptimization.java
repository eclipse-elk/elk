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
package org.eclipse.elk.alg.radial.intermediate.optimization;

import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.graph.ElkEdge;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.util.ElkGraphUtil;

/**
 * Calculate the edge length between the root node and the first radius nodes. The edge source is given by a position.
 * The algorithm expects the {@link CoreOptions.POSITION} option to be set. <em>Warning</em>: It makes the assumption
 * that the position points to a node <em>in</em> the tree-parent node!
 */
public class EdgeLengthPositionOptimization implements IEvaluation {
    /**
     * This evaluation method calculates the sum of the length of each outgoing edge from the root node. The edge
     * starting at the given position.
     */
    @Override
    public double evaluate(final ElkNode root) {
        double edgeLength = 0;
        for (ElkEdge edge : ElkGraphUtil.allOutgoingEdges(root)) {
            ElkNode target = ElkGraphUtil.connectableShapeToNode(edge.getTargets().get(0));
            double targetX = target.getX() + target.getWidth() / 2;
            double targetY = target.getY() + target.getHeight() / 2;

            KVector position = target.getProperty(CoreOptions.POSITION);
            double rootX = root.getX() + position.x + root.getWidth() / 2;
            double rootY = root.getY() + position.y + root.getHeight();

            double vectorX = targetX - rootX;
            double vectorY = targetY - rootY;
            edgeLength += Math.sqrt(vectorX * vectorX + vectorY * vectorY);
        }
        return edgeLength;
    }
}