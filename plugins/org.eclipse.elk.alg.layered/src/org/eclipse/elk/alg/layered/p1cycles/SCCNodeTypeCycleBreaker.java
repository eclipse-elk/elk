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
import org.eclipse.elk.alg.layered.options.LayeredOptions;

/**
 * This Cycle Breaking Strategy extends the SCCModelOrderCycleBreaker. The preferred node type for the minimum or 
 * maximum node can be defined.
 *
 */
public class SCCNodeTypeCycleBreaker extends SCCModelOrderCycleBreaker {

    @Override
    public void findNodes(int offset) {
     // lowest model order only outgoing
      for (int i = 0; i < stronglyConnectedComponents.size(); i++) {
          if (stronglyConnectedComponents.get(i).size() <= 1) {
              continue;
          }
          LNode min = null;
          int maxGroupModelOrder = 100; // This is a hardcoded value.
          int modelOrderMin = Integer.MAX_VALUE;
          int modelOrderMax = Integer.MIN_VALUE;
          int groupModelOrderMin = Integer.MAX_VALUE;
          int groupModelOrderMax = Integer.MIN_VALUE;
          // FIXME explain what is happening here.
          LNode max = null;
          for (LNode n : stronglyConnectedComponents.get(i)) {
              if (min == null || max == null) {
                  min = n;
                  modelOrderMin = computeConstraintModelOrder(n,offset);
                  groupModelOrderMin = computeConstraintGroupModelOrder(n, maxGroupModelOrder);
                  max = n;
                  modelOrderMax = modelOrderMin;
                  groupModelOrderMax = groupModelOrderMin;
                  continue;
              }
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
          if (min != null) {
              if (min.getProperty(LayeredOptions.CONSIDER_MODEL_ORDER_GROUP_MODEL_ORDER_CYCLE_BREAKING_ID) == 1) { // FIXME
                  for (LEdge edge : min.getIncomingEdges()) {
                      if (stronglyConnectedComponents.get(i).contains(edge.getSource().getNode())) {
                          if(edge.getTarget().getNode().hasProperty(LayeredOptions.LAYERING_LAYER_CONSTRAINT)) {
                              continue;
                          }
                          revEdges.add(edge);
                      }
                  }
              }
              else {
                  for (LEdge edge : max.getIncomingEdges()) {
                      if (stronglyConnectedComponents.get(i).contains(edge.getSource().getNode())) {
                          if(edge.getTarget().getNode().hasProperty(LayeredOptions.LAYERING_LAYER_CONSTRAINT)) {
                              continue;
                          }
                          revEdges.add(edge);
                      }
                  }
              }
          }
      }
    }
}