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

import java.util.Arrays;
import java.util.Comparator;
import java.util.List;

import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.LayeredOptions;

/**
 * Greedy Cycle Breaker that behaves the same as {@link GreedyCycleBreaker} but does not randomly choose an edge to
 * reverse if multiple candidates exist but does so by model order. Now with Group Model Order.
 */
public final class GreedyModelOrderCycleBreaker extends GreedyCycleBreaker {
    
    public static final int START_UP_GROUP_ID = 0;
    public static final int REACTOR_GROUP_ID = 1;
    public static final int ACTION_GROUP_ID = 2;
    public static final int TIMER_GROUP_ID = 3;
    public static final int REACTION_GROUP_ID = 4;
    public static final int DUMMY_GROUP_ID = 5;
    public static final int SHUT_DOWN_GROUP_ID = 6;
    
    /**
     * Choose the node with the minimum model order.
     */
    @Override
    protected LNode chooseNodeWithMaxOutflow(final List<LNode> nodes) {
        LNode returnNode = null;
        int minimumModelOrder = Integer.MAX_VALUE;
        for (LNode node : nodes) {
            if (node.hasProperty(InternalProperties.MODEL_ORDER)) {
                if (returnNode == null ||
                   (node.getProperty(InternalProperties.MODEL_ORDER) < minimumModelOrder &&
                     GROUP_ORDER.compare(node, returnNode) <=0)||
                   GROUP_ORDER.compare(node, returnNode) <0) {
                        minimumModelOrder = node.getProperty(InternalProperties.MODEL_ORDER);
                        returnNode = node;
                }       
            } else if (returnNode != null) {
                if (GROUP_ORDER.compare(node, returnNode) <0) {
                    minimumModelOrder = Integer.MAX_VALUE;
                    returnNode = node;
                }
            }
        }
        if (returnNode == null) {
            return super.chooseNodeWithMaxOutflow(nodes);
        }
        return returnNode;
    }

    /**
     * Define a static order between different Groups.
     */
    public static final Comparator<LNode> GROUP_ORDER =
            new Comparator<LNode>() {
          @Override
          public int compare(LNode node1, LNode node2) {
            assert(node1.hasProperty(LayeredOptions.CONSIDER_MODEL_ORDER_GROUP_I_D));
            assert(node2.hasProperty(LayeredOptions.CONSIDER_MODEL_ORDER_GROUP_I_D));
            int group1 = node1.getProperty(LayeredOptions.CONSIDER_MODEL_ORDER_GROUP_I_D);
            int group2 = node2.getProperty(LayeredOptions.CONSIDER_MODEL_ORDER_GROUP_I_D);
            if (group1 == group2) {
              return 0;
            }
            List<Integer> highestPriority = Arrays.asList(START_UP_GROUP_ID,TIMER_GROUP_ID);
            if(highestPriority.contains(group1)) return -1;
            if(highestPriority.contains(group2)) return 1;
            if(group1 == SHUT_DOWN_GROUP_ID) return 1;
            if(group2 == SHUT_DOWN_GROUP_ID) return -1;
            if(group1 == REACTION_GROUP_ID) return -1;
            if(group2 == REACTION_GROUP_ID) return 1;
            if(group1 == REACTOR_GROUP_ID) return -1;
            if(group2 == REACTOR_GROUP_ID) return 1;
            if(group1 == ACTION_GROUP_ID) return -1;
            if(group2 == ACTION_GROUP_ID) return 1;
            assert(false);
            return 0;
          }
        };
}
