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

import java.util.LinkedList;
import java.util.List;

import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.options.LayeredOptions;

/**
 * This Cycle Breaking Strategy extends the SCCModelOrderCycleBreaker, however instead of reversing all outgoing edges
 * of the node with the highest model order, it reverses all incoming edges of the node with the lowest model order.
 *
 */
public class SCCLowestModelOrderCycleBreaker extends SCCModelOrderCycleBreaker {

    @Override
    public void findNodes(int offset) {
     // lowest model order only outgoing
      for (int i = 0; i < stronglyConnectedComponents.size(); i++) {
          if (stronglyConnectedComponents.get(i).size() <= 1) {
              continue;
          }
          LNode min = null;
          int modelOrderMin = Integer.MAX_VALUE;
          for (LNode n : stronglyConnectedComponents.get(i)) {
              List<Integer> layermask = new LinkedList<Integer>();
              layermask = this.graph.getProperty(LayeredOptions.CONSIDER_MODEL_ORDER_GROUP_MODEL_ORDER_CB_ENFORCED_GROUP_ORDERS);
              int groupID = n.getProperty(LayeredOptions.CONSIDER_MODEL_ORDER_GROUP_MODEL_ORDER_CYCLE_BREAKING_ID);
              if (!layermask.contains(groupID)) {
                  continue;
              }
              if (min == null) {
                  min = n;
                  modelOrderMin = computeConstraintModelOrder(n,offset);
                  continue;
              }
              int modelOrderCurrent = computeConstraintModelOrder(n, offset);
              if (modelOrderMin > modelOrderCurrent) {
                  min = n;
                  modelOrderMin = modelOrderCurrent;
              }
          }
          System.out.println("MIN:" + min.toString() 
          + "MODEL ORDER: " + modelOrderMin + "Of SCC: " + i);
          if (min != null) {
              for (LEdge edge : min.getIncomingEdges()) {
                  if (stronglyConnectedComponents.get(i).contains(edge.getSource().getNode())) {
                      revEdges.add(edge);
                      System.out.println("Reversed: " + edge.toString());
                  }

             }
          }
      }
    }
}