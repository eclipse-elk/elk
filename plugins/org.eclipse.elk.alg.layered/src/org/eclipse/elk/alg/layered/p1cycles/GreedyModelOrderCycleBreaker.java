/*******************************************************************************
 * Copyright (c) 2022 Kiel University and others.
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
import org.eclipse.elk.alg.layered.options.InternalProperties;

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
        for (LNode node : nodes) {
            // In this step nodes without a model order are disregarded.
            // One could of course think of a different strategy regarding this aspect.
            // FUTURE WORK: If multiple model order groups exist, one has to chose based on the priority of the groups.
            // FIXME
            if (node.hasProperty(InternalProperties.MODEL_ORDER)
                    && node.getProperty(InternalProperties.MODEL_ORDER) < minimumModelOrder) {
                minimumModelOrder = node.getProperty(InternalProperties.MODEL_ORDER);
                returnNode = node;
            }
        }
        if (returnNode == null) {
            return super.chooseNodeWithMaxOutflow(nodes);
        }
        return returnNode;
    }

}
