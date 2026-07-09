/*******************************************************************************
 * Copyright (c) 2016 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.

 * This Source Code may also be made available under the following Secondary
 * Licenses when the conditions for such availability set forth in the Eclipse
 * Public License v. 2.0 are satisfied: GPL-3.0 which is available at
 * https://www.gnu.org/licenses/gpl-3.0-standalone.html.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.layered.options;

import org.eclipse.elk.graph.properties.AdvancedPropertyValue;

/**
 * Definitions of strategies for the node promotion heuristic.
 * 
 * We recommend to use {@link #NO_BOUNDARY} if a small number of dummy nodes, i.e. short edges, is
 * desired. If the width (in a top-down layout) of a previous layering should not increase, we
 * recommend {@link #NIKOLOV_IMPROVED}.
 * 
 * Note that the term 'width' here refers to the maximum number of nodes in any layer.
 */
public enum NodePromotionStrategy {

    /**
     * Node promotion won't be applied to the graph.
     */
    NONE,
    /**
     * Strategy for promotion proposed by Nikolov et al. in their paper to keep the layering at a
     * width that's at worst as wide as the layering of the used layering algorithm.
     */
    NIKOLOV,
    /**
     * Strategy of Nikolov et al. (maintaining width while promoting nodes) but transferred to
     * approximated sizes of original and dummy nodes in pixels.
     */
    NIKOLOV_PIXEL,
    /**
     * This strategy tries a complete run of the node promotion heuristic without any boundaries and
     * checks afterwards if the maximal width has been exceeded. Only if it is exceeded, the result
     * will be dismissed and another run of the node promotion heuristic is started with the
     * strategy of Nikolov et al. (maintaining width while promoting nodes).
     */
    @AdvancedPropertyValue
    NIKOLOV_IMPROVED,
    /**
     * Operates like NIKOLOV_IMPROVED but with consideration of approximated sizes of the original
     * and dummy nodes.
     */
    @AdvancedPropertyValue
    NIKOLOV_IMPROVED_PIXEL,
    /**
     * Stops the promotion either if there are no more promotions to make or a certain percentage (given through
     * {@link LayeredOptions#LAYERING_NODE_PROMOTION_MAX_ITERATIONS}) of dummy nodes has been reduced while promoting.
     * If the percentage is reached there will be no more promotions made even if there are possible promotions left.
     */
    @AdvancedPropertyValue
    DUMMYNODE_PERCENTAGE,
    /**
     * Stops the promotion either if there are no more promotions to make or a certain number of
     * promotions were executed. The number is calculated by a percentage of
     * {@link LayeredOptions#LAYERING_NODE_PROMOTION_MAX_ITERATIONS} of the number of all nodes in the graph. 
     * If the number is reached there will be no more promotions made even if there are possible promotions left.
     */
    @AdvancedPropertyValue
    NODECOUNT_PERCENTAGE,
    /** 
     * The node promotion will run until there are no more promotions left to make. 
     */
    @AdvancedPropertyValue
    NO_BOUNDARY,
    /**
     * Promotes such that they are in a layer that conforms to the model order (if it is applicable).
     * 1 -> 2 -> 3 -> 5 is transformed to 1 -> 2 -> 3 -> 5
     *   -> 4 ------>                       ------> 4 ->
     * Since all nodes in layer 2 (node 2) are smaller and the next layer also has a smaller node,
     * 4 can be moved more to the right.
     * Model order does not care about maximum width but only about the model order.
     * Should be used together with LONGEST_PATH_TO_SOURCE layering.
     */
    MODEL_ORDER_LEFT_TO_RIGHT,
    /**
     * Promotes such that they are in a layer that conforms to the model order (if it is applicable).
     * 1 -> 3 -> 4 -> 5 is transformed to 1 -> 2 ------> 5
     *   ------> 2 ->                       -> 3 -> 4 ->
     * Since all nodes in layer 3 (node 3) are bigger and one in layer 2 is bigger, 2 can be moved more to the left.
     * Model order does not care about maximum width but only about the model order.
     * Should be used together with LONGEST_PATH layering.
     */
    MODEL_ORDER_RIGHT_TO_LEFT;
}
