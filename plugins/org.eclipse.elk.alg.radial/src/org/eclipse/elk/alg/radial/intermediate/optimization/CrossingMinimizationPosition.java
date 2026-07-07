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

import java.util.List;

import org.eclipse.elk.alg.radial.RadialUtil;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.graph.ElkNode;

/**
 * Count the number of crossings of the edges between root and the first radius. The algorithm expects the
 * {@link CoreOptions.POSITION} option to be set. <em>Warning</em>: It makes assumptions that the position points to a node <em>in</em>
 * the tree-parent node!
 */
public class CrossingMinimizationPosition implements IEvaluation {

    /** The root node of the graph. */
    private ElkNode root;

    @Override
    public double evaluate(final ElkNode rootNode) {
        this.root = rootNode;
        int crossings = 0;
        List<ElkNode> nodes = RadialUtil.getSuccessors(rootNode);

        int k = 0;
        for (ElkNode node1 : nodes) {
            k++;
            for (int i = k; i < nodes.size(); i++) {
                if (isCrossing(node1, nodes.get(i))) {
                    crossings += 1;
                }
            }
        }
        return crossings;
    }

    /**
     * Calculate if two edges are crossing by calculating the line between the nodes and their node position.
     * 
     * @param node1
     *            The first node for checking the crossing.
     * @param node2
     *            The second node for checking the crossing.
     * @return Returns if the lines of the nodes are crossing.
     */
    private boolean isCrossing(final ElkNode node1, final ElkNode node2) {
        // root
        double rootX = root.getX() + root.getWidth() / 2;
        double rootY = root.getX() + root.getWidth() / 2;

        // node1
        double xPos1 = node1.getX() + node1.getWidth() / 2;
        double yPos1 = node1.getY() + node1.getHeight() / 2;
        KVector node1Vector = new KVector(xPos1, yPos1);

        KVector position1 = node1.getProperty(CoreOptions.POSITION);
        position1.x = position1.x + rootX;
        position1.y = position1.y + rootY;

        double m1 = (node1Vector.y - position1.y) / (node1Vector.x - position1.x);
        double b1 = node1Vector.y - m1 * node1Vector.x;

        // node2
        double xPos2 = node2.getX() + node2.getWidth() / 2;
        double yPos2 = node2.getY() + node2.getHeight() / 2;
        KVector node2Vector = new KVector(xPos2, yPos2);

        KVector position2 = node2.getProperty(CoreOptions.POSITION);
        position2.x = position2.x + rootX;
        position2.y = position2.y + rootY;

        double m2 = (node2Vector.y - position2.y) / (node2Vector.x - position2.x);
        double b2 = node2Vector.y - m2 * node2Vector.x;

        double xCut = (b1 - b2) / (m2 - m1);
        // check whether the cut occurs on the relevant line segment
        if ((position1.x < xCut && node1Vector.x < xCut) || (xCut < position1.x && xCut < node1Vector.x)) {
            return false;
        } else if ((position2.x < xCut && node2Vector.x < xCut) || (xCut < position2.x && xCut < node2Vector.x)) {
            return false;
        }
        return true;
    }

}