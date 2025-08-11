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
 * Based on the SCCModelOrderCycleBreaker. This finds the nodes with minimum and maximum model order and reverses the 
 * incoming nodes of the minimum, if its the in-degree is greater than the out-degree of the maximum node. Else it 
 * reverses the out-going edges of the maximum node.
 * If group model order should be enforced, this uses the group model order as a primary criterion.
 * 
 * <dl>
 *   <dt>Precondition:</dt>
 *      <dd>no self loops</dd>
 *   <dt>Postcondition:</dt>
 *      <dd>the graph has no cycles</dd>
 * </dl>
 *
 */
public class SCConnectivity extends SCCModelOrderCycleBreaker {

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
            // If the minimum node has more incoming edges than the maximum node has outgoing edges,
            // reverse all edges to the minimum node and remove it from the strongly connected component.
            // If it is the other way around, reverse all outgoing edges of the maximum node.
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