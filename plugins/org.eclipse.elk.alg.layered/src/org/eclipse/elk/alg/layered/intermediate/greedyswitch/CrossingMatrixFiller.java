/*******************************************************************************
 * Copyright (c) 2015 Kiel University and others.
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
package org.eclipse.elk.alg.layered.intermediate.greedyswitch;

import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.intermediate.greedyswitch.SwitchDecider.CrossingCountSide;
import org.eclipse.elk.alg.layered.p3order.LayerSweepCrossingMinimizer.CrossMinType;

/**
 * This class manages the crossing matrix and fills it on demand. It needs to be reinitialized for
 * each free layer. For each layer the node.id fields MUST be set from 0 to layer.getSize() - 1!
 */
public final class CrossingMatrixFiller {
    private final boolean[][] isCrossingMatrixFilled;
    private final int[][] crossingMatrix;
    private final BetweenLayerEdgeTwoNodeCrossingsCounter inBetweenLayerCrossingCounter;
    private final CrossingCountSide direction;
    private final boolean oneSided;

    // SUPPRESS CHECKSTYLE NEXT 30 Javadoc
    /**
     * Constructs class which manages the crossing matrix.
     */
    public CrossingMatrixFiller(final CrossMinType greedySwitchType, final LNode[][] graph,
            final int freeLayerIndex, final CrossingCountSide direction) {
        
        this.direction = direction;
        oneSided = greedySwitchType == CrossMinType.ONE_SIDED_GREEDY_SWITCH;

        LNode[] freeLayer = graph[freeLayerIndex];
        isCrossingMatrixFilled = new boolean[freeLayer.length][freeLayer.length];
        crossingMatrix = new int[freeLayer.length][freeLayer.length];

        inBetweenLayerCrossingCounter = new BetweenLayerEdgeTwoNodeCrossingsCounter(graph, freeLayerIndex);
    }

    /**
     * Returns entry for crossings between edges incident to two nodes, where upperNode is above
     * lowerNode in the layer.
     */
    public int getCrossingMatrixEntry(final LNode upperNode, final LNode lowerNode) {
        if (!isCrossingMatrixFilled[upperNode.id][lowerNode.id]) {
            fillCrossingMatrix(upperNode, lowerNode);
            isCrossingMatrixFilled[upperNode.id][lowerNode.id] = true;
            isCrossingMatrixFilled[lowerNode.id][upperNode.id] = true;
        }
        return crossingMatrix[upperNode.id][lowerNode.id];
    }

    private void fillCrossingMatrix(final LNode upperNode, final LNode lowerNode) {
        if (oneSided) {
            switch (direction) {
            case EAST:
                inBetweenLayerCrossingCounter.countEasternEdgeCrossings(upperNode, lowerNode);
                break;
            case WEST:
                inBetweenLayerCrossingCounter.countWesternEdgeCrossings(upperNode, lowerNode);
            }
        } else {
            inBetweenLayerCrossingCounter.countBothSideCrossings(upperNode, lowerNode);
        }
        crossingMatrix[upperNode.id][lowerNode.id] =
                inBetweenLayerCrossingCounter.getUpperLowerCrossings();
        crossingMatrix[lowerNode.id][upperNode.id] =
                inBetweenLayerCrossingCounter.getLowerUpperCrossings();
    }
}
