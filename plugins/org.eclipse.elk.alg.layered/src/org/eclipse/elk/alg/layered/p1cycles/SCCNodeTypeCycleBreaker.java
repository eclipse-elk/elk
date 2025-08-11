/*******************************************************************************
 * Copyright (c) 2025 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 * 
 * SPDX-License-Identifier: EPL-2.0 
 *******************************************************************************/
package org.eclipse.elk.alg.layered.p1cycles;

import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.options.GroupOrderStrategy;
import org.eclipse.elk.alg.layered.options.LayeredOptions;

import com.google.common.collect.Iterables;

/**
 * This cycle breaking strategy extends the {@link org.eclipse.elk.alg.layered.p1cycles.SCCModelOrderCycleBreaker}.
 * The preferred node type for the minimum or maximum node can be defined.
 * 
 * <dl>
 *   <dt>Precondition:</dt>
 *      <dd>no self loops</dd>
 *   <dt>Postcondition:</dt>
 *      <dd>the graph has no cycles</dd>
 * </dl>
 *
 */
public class SCCNodeTypeCycleBreaker extends SCCModelOrderCycleBreaker {

    @Override
    public void findNodes(int offset, int bigOffset) {
        for (int i = 0; i < stronglyConnectedComponents.size(); i++) {
            // Nothing needs to be done if only one strongly connected component exists.
            if (stronglyConnectedComponents.get(i).size() <= 1) {
                continue;
            }
            LNode min = null;
            LNode max = null;
            int modelOrderMin = Integer.MAX_VALUE;
            int modelOrderMax = Integer.MIN_VALUE;
            // Check whether model order or the group model order is the primary criterion.
            // If it is group model order, we need to handle this differently.
            boolean enforceGroupModelOrder = this.graph.getProperty(
                    LayeredOptions.CONSIDER_MODEL_ORDER_GROUP_MODEL_ORDER_CB_GROUP_ORDER_STRATEGY) == GroupOrderStrategy.ENFORCED;
            GroupModelOrderCalculator calculator = new GroupModelOrderCalculator();
            // Iterate over all strongly connected components to find the maximum/minimum node to reverse edges.
            for (LNode n : stronglyConnectedComponents.get(i)) {
                // First calculate initial values
                if (min == null || max == null) {
                    min = n;
                    modelOrderMin = enforceGroupModelOrder
                            ? calculator.computeConstraintGroupModelOrder(n, bigOffset, offset)
                            : calculator.computeConstraintModelOrder(n, offset);
                    max = n;
                    modelOrderMax = modelOrderMin;
                } else {
                    // For all not first nodes find the group model order and model order with constraints and update the
                    // minimum and maximum node.
                    int modelOrderCurrent = enforceGroupModelOrder
                            ? calculator.computeConstraintGroupModelOrder(n, bigOffset, offset)
                            : calculator.computeConstraintModelOrder(n, offset);
                    if (modelOrderMin > modelOrderCurrent) {
                        min = n;
                        modelOrderMin = modelOrderCurrent;
                    }
                    if (modelOrderMax < modelOrderCurrent) {
                        max = n;
                        modelOrderMax = modelOrderCurrent;
                    }
                }
            }
            if (min.getProperty(LayeredOptions.CONSIDER_MODEL_ORDER_GROUP_MODEL_ORDER_CYCLE_BREAKING_ID) == graph
                    .getProperty(LayeredOptions.CONSIDER_MODEL_ORDER_GROUP_MODEL_ORDER_CB_PREFERRED_SOURCE_ID)) {
                // If the minimum model order node is a preferred source node, reverse all its incoming edges.
                for (LEdge edge : min.getIncomingEdges()) {
                    if (stronglyConnectedComponents.get(i).contains(edge.getSource().getNode())) {
                        revEdges.add(edge);
                    }
                }
            } else if (max.getProperty(LayeredOptions.CONSIDER_MODEL_ORDER_GROUP_MODEL_ORDER_CYCLE_BREAKING_ID) == graph
                    .getProperty(LayeredOptions.CONSIDER_MODEL_ORDER_GROUP_MODEL_ORDER_CB_PREFERRED_TARGET_ID)) {
                // If the maximum model order node is a preferred node reverse all its outgoing edges.
                for (LEdge edge : max.getOutgoingEdges()) {
                    if (stronglyConnectedComponents.get(i).contains(edge.getSource().getNode())) {
                        revEdges.add(edge);
                    }
                }
            } else {
                // Use connectivity to decide.
                if (Iterables.size(min.getIncomingEdges()) > Iterables.size(max.getOutgoingEdges())) {
                    for (LEdge edge : min.getIncomingEdges()) {
                        if (stronglyConnectedComponents.get(i).contains(edge.getSource().getNode())) {
                            revEdges.add(edge);
                        }
                    }
                } else {
                    for (LEdge edge : max.getOutgoingEdges()) {
                        if (stronglyConnectedComponents.get(i).contains(edge.getTarget().getNode())) {
                            revEdges.add(edge);
                        }
                    }
                }
            }
        }
    }
}