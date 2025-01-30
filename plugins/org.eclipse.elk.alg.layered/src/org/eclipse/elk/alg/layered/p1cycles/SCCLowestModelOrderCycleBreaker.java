/*******************************************************************************
 * Copyright (c) 2023 mwr and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 * 
 * SPDX-License-Identifier: EPL-2.0 
 *******************************************************************************/
package org.eclipse.elk.alg.layered.p1cycles;

import java.util.Arrays;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;
import java.util.Set;
import java.util.Stack;

import org.eclipse.elk.alg.layered.LayeredPhases;
import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.intermediate.IntermediateProcessorStrategy;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.core.alg.ILayoutPhase;
import org.eclipse.elk.core.alg.LayoutProcessorConfiguration;
import org.eclipse.elk.core.util.IElkProgressMonitor;

import com.google.common.collect.Lists;

/**
 * This Cycle Breaking Strategy extends the SCCModelOrderCycleBreaker, however instead of reversing all outgoing edges
 * of the node with the highest model order, it reverses all incoming edges of the node with the lowest model order.
 * @author mwr
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
              layermask.add(1);
              layermask.add(4);
              int groupID = n.getProperty(LayeredOptions.CONSIDER_MODEL_ORDER_GROUP_I_D);
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