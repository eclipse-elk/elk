/*******************************************************************************
 * Copyright (c) 2025 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 * 
 * SPDX-License-Identifier: EPL-2.0 
 *******************************************************************************/
package org.eclipse.elk.alg.layered.intermediate.preserveorder;

import java.util.List;

import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LGraphElement;
import org.eclipse.elk.alg.layered.options.GroupOrderStrategy;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.LayeredOptions;

/**
 * Helper class to calculate the model order and group model order for crossing minimization.
 */
public class CMGroupModelOrderCalculator {
    
    public static int calculateModelOrderOrGroupModelOrder(LGraphElement element, LGraphElement other, LGraph parent, int offset) {
        boolean enforceGroupModelOrder = parent.getProperty(
                LayeredOptions.CONSIDER_MODEL_ORDER_GROUP_MODEL_ORDER_CM_GROUP_ORDER_STRATEGY) == GroupOrderStrategy.ENFORCED;
        List<Integer> enforcedOrders = parent.getProperty(LayeredOptions.CONSIDER_MODEL_ORDER_GROUP_MODEL_ORDER_CM_ENFORCED_GROUP_ORDERS);
        if (!element.hasProperty(InternalProperties.MODEL_ORDER)) {
            return -1;
        } else if (enforceGroupModelOrder) {
            // Check whether the order between both is enforced. Otherwise use model order.
            if (enforcedOrders.contains(element.getProperty(LayeredOptions.CONSIDER_MODEL_ORDER_GROUP_MODEL_ORDER_CROSSING_MINIMIZATION_ID))
                && enforcedOrders.contains(other.getProperty(LayeredOptions.CONSIDER_MODEL_ORDER_GROUP_MODEL_ORDER_CROSSING_MINIMIZATION_ID))) {
                // This means that we need to calculate the model order by groupId * number of nodes/port/edges + real model order.
                return offset
                        * element.getProperty(LayeredOptions.CONSIDER_MODEL_ORDER_GROUP_MODEL_ORDER_CROSSING_MINIMIZATION_ID)
                        + element.getProperty(InternalProperties.MODEL_ORDER);
            }
            // Fallthrough case
        } else {
            // Case only model order.
            return element.getProperty(InternalProperties.MODEL_ORDER);
        }
        return element.getProperty(InternalProperties.MODEL_ORDER);
    }

}
