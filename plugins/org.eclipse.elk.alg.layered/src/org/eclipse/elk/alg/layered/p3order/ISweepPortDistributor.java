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
package org.eclipse.elk.alg.layered.p3order;

import java.util.Random;

import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.p3order.LayerSweepCrossingMinimizer.CrossMinType;
import org.eclipse.elk.alg.layered.p3order.counting.IInitializable;

/**
 * PortDistributor to be used while sweeping in phase 3.
 * <p>
 * Must be initialized using {@link IInitializable#init(java.util.List)}!
 * </p>
 */
public interface ISweepPortDistributor extends IInitializable {

    /**
     * Distribute ports in one layer. To be used in the context of layer sweep.
     *
     * @param order
     *            the current order of the nodes
     * @param freeLayerIndex
     *            the index of the layer the node is in
     * @param isForwardSweep
     *            whether we are sweeping forward or not.
     */
    boolean distributePortsWhileSweeping(LNode[][] order, int freeLayerIndex,
            boolean isForwardSweep);
    
    /**
     * Make a port distributor.
     * 
     * @param cmt
     *            the crossing minimization type
     * @param r
     *            random number generator for this graph
     * @param currentOrder
     *            the current order of the nodes.
     * @return the port distributor
     */
    static ISweepPortDistributor create(final CrossMinType cmt, final Random r, final LNode[][] currentOrder) {
        if (cmt == CrossMinType.TWO_SIDED_GREEDY_SWITCH) {
            return new GreedyPortDistributor();
        } else if (r.nextBoolean()) {
            // Since both methods lead to different results, but neither is clearly better, we choose randomly.
            return new NodeRelativePortDistributor(currentOrder.length);
        } else {
            return new LayerTotalPortDistributor(currentOrder.length);
        }
    }
}

