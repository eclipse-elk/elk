/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.radial.intermediate.optimization;

import org.eclipse.elk.core.math.ElkMath;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.graph.ElkEdge;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.util.ElkGraphUtil;

/**
 * The edge length optimization provides an evaluation method which allows to compare different layouts by the length of all
 * edges.
 */
public class EdgeLengthOptimization implements IEvaluation {

    /**
     * This evaluation method calculates the sum of the length of each outgoing edge from the root node.
     */
    @Override
    public double evaluate(final ElkNode root) {
        double edgeLength = 0;
        for (ElkEdge edge : ElkGraphUtil.allOutgoingEdges(root)) {
            ElkNode target = ElkGraphUtil.connectableShapeToNode(edge.getTargets().get(0));

            double targetX = target.getX() + target.getWidth() / 2;
            double targetY = target.getY() + target.getHeight() / 2;

            double rootX = root.getX() + root.getWidth() / 2;
            double rootY = root.getY() + root.getHeight() / 2;

            // Clipping
            KVector vector = new KVector();
            vector.x = targetX - rootX;
            vector.y = targetY - rootY;
            KVector sourceClip = new KVector(vector.x, vector.y);
            ElkMath.clipVector(sourceClip, root.getWidth(), root.getHeight());
            vector.x -= sourceClip.x;
            vector.y -= sourceClip.y;

            rootX = targetX - vector.x;
            rootY = targetY - vector.y;

            KVector targetClip = new KVector(vector.x, vector.y);
            ElkMath.clipVector(targetClip, target.getWidth(), target.getHeight());
            vector.x -= targetClip.x;
            vector.y -= targetClip.y;

            targetX = rootX + vector.x;
            targetY = rootY + vector.y;

            double vectorX = targetX - rootX;
            double vectorY = targetY - rootY;
            edgeLength += Math.sqrt(vectorX * vectorX + vectorY * vectorY);
        }
        return edgeLength;
    }
}