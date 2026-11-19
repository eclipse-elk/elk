/*******************************************************************************
 * Copyright (c) 2022 - 2023 Kiel University and others.
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
package org.eclipse.elk.alg.topdownpacking;

import org.eclipse.elk.core.alg.ILayoutPhase;
import org.eclipse.elk.core.alg.ILayoutPhaseFactory;

/**
 * Node arrangement strategy to use during topdown packing algorithm.
 *
 */
public enum NodeArrangementStrategy implements ILayoutPhaseFactory<TopdownPackingPhases, GridElkNode> {
    /**
     * Places nodes from left to right, top to bottom in a grid.
     */
    LEFT_RIGHT_TOP_DOWN_NODE_PLACER;

    /**
     * {@inheritDoc}
     */
    @Override
    public ILayoutPhase<TopdownPackingPhases, GridElkNode> create() {
        switch (this) {
        case LEFT_RIGHT_TOP_DOWN_NODE_PLACER:
            return new LeftRightTopDownNodePlacer();
        default:
            throw new IllegalArgumentException(
                    "No implementation is available for the node placer " + this.toString());
        }
    }

}
