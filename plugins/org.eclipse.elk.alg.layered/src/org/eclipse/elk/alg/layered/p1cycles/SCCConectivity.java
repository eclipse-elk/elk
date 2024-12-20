/*******************************************************************************
 * Copyright (c) 2024 Sasuk and others.
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

import com.google.common.collect.Iterables;

/**
 * Based on the SCCModelOrderCycleBreaker. This finds the nodes with minimum and maximum model order and reverses the 
 * incoming nodes of the minimum, if its the in-degree is greater than the out-degree of the maximum node. Else it 
 * reverses the out-going edges of the maximum node.
 * @author Mwr
 *
 */
public class SCCConectivity extends SCCModelOrderCycleBreaker {

    @Override
    public void findNodes(int offset) {
      for (int i = 0; i < stronglyConnectedComponents.size(); i++) {
          if (stronglyConnectedComponents.get(i).size() <= 1) {
              continue;
          }
          LNode min = null;
          LNode max = null;
          int maxGroupModelOrder = 30; // This makes calculating the ordering values easier. FIXME
          int modelOrderMin = Integer.MAX_VALUE;
          int modelOrderMax = Integer.MIN_VALUE;
          int groupModelOrderMin = Integer.MAX_VALUE;
          int groupModelOrderMax = Integer.MIN_VALUE;
          // Iterate over all stringly connected components to find the maximum/minimum node to reverse edges.
          for (LNode n : stronglyConnectedComponents.get(i)) {
              // First calculate initial values FIXME why???
              if (min == null || max == null) {
                  min = n;
                  modelOrderMin = computeConstraintModelOrder(n,offset);
                  groupModelOrderMin = computeConstraintGroupModelOrder(n, maxGroupModelOrder);
                  max = n;
                  modelOrderMax = modelOrderMin;
                  groupModelOrderMax = groupModelOrderMin;
                  continue;
              }
              // For all not first nodes find the group model order and model order with constraints and update the
              // minimum and maximum node.
              int modelOrderCurrent = computeConstraintModelOrder(n, offset);
              int groupModelOrderCurrent = computeConstraintGroupModelOrder(n, maxGroupModelOrder);
              if (groupModelOrderMin < groupModelOrderCurrent) {
                  min = n;
                  groupModelOrderMin = groupModelOrderCurrent;
                  modelOrderMin = modelOrderCurrent;
              } else if (groupModelOrderMin == groupModelOrderCurrent) {                  
                  if (modelOrderMin > modelOrderCurrent) {
                      min = n;
                      modelOrderMin = modelOrderCurrent;
                  }
              } else if (groupModelOrderMax == groupModelOrderCurrent) {
                  if (modelOrderMax < modelOrderCurrent) {
                      max = n;
                      modelOrderMax = modelOrderCurrent;
                  }
              } else if (groupModelOrderMax > groupModelOrderCurrent) {
                  max = n;
                  groupModelOrderMax = groupModelOrderCurrent;
                  modelOrderMax = modelOrderCurrent;
              }
          }
          // If a minimum and maximum node are found, there must a node for which we reverse the edges.
          if (min != null && max != null) {
              // If the minimum node has more incoming edges than the maximum node has outgoing edges,
              // reverse all edges to the minimum node and remove it from the strongly connected component.
              // If it is the other wayaround, reverse all outgoing edges of the maximum node.
              if (Iterables.size(min.getIncomingEdges()) > 
                      Iterables.size(max.getOutgoingEdges())) {
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
