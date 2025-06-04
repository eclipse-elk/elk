/*******************************************************************************
 * Copyright (c) 2022, 2025 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 * 
 * SPDX-License-Identifier: EPL-2.0 
 *******************************************************************************/
package org.eclipse.elk.alg.layered.p1cycles;

import java.util.List;

import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.options.GroupOrderStrategy;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.LayeredOptions;

/**
 * Greedy Cycle Breaker that behaves the same as {@link GreedyCycleBreaker} but does not randomly choose an edge to
 * reverse if multiple candidates exist but does so by model order.
 */
public final class GreedyModelOrderCycleBreaker extends GreedyCycleBreaker {
    
    /**
     * Choose the node with the minimum model order.
     */
    @Override
    protected LNode chooseNodeWithMaxOutflow(final List<LNode> nodes) {
        LNode returnNode = null;
        int minimumModelOrder = Integer.MAX_VALUE;
        int offset = Math.max(layeredGraph.getLayerlessNodes().size(), layeredGraph.getProperty(InternalProperties.MAX_MODEL_ORDER_NODES));
        int bigOffset = offset * layeredGraph.getProperty(InternalProperties.CB_NUM_MODEL_ORDER_GROUPS);
        GroupModelOrderCalculator moCalculator = new GroupModelOrderCalculator();
        boolean enforceGroupModelOrder = layeredGraph.getProperty(
                LayeredOptions.CONSIDER_MODEL_ORDER_GROUP_MODEL_ORDER_CB_GROUP_ORDER_STRATEGY) == GroupOrderStrategy.ENFORCED;
        for (LNode node : nodes) {
            // In this step nodes without a model order are disregarded.
            // One could of course think of a different strategy regarding this aspect.
            if (node.hasProperty(InternalProperties.MODEL_ORDER)) {
                int modelOrder = enforceGroupModelOrder
                ? moCalculator.computeConstraintGroupModelOrder(node, bigOffset, offset)
                : moCalculator.computeConstraintModelOrder(node, offset);
                if (minimumModelOrder > modelOrder) {                    
                    minimumModelOrder = modelOrder;
                    returnNode = node;
                }
            }
        }
        if (returnNode == null) {
            return super.chooseNodeWithMaxOutflow(nodes);
        }
        return returnNode;
    }

}
