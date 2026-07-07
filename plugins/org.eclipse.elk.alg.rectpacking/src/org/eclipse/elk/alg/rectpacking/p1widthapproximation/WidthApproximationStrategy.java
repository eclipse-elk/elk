/*******************************************************************************
 * Copyright (c) 2022 Kiel University and others.
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
package org.eclipse.elk.alg.rectpacking.p1widthapproximation;

import org.eclipse.elk.alg.rectpacking.RectPackingLayoutPhases;
import org.eclipse.elk.core.alg.ILayoutPhase;
import org.eclipse.elk.core.alg.ILayoutPhaseFactory;
import org.eclipse.elk.graph.ElkNode;

/**
 * Strategy factory to determine a target width to transform the rectangle placement problem into a strip-packing
 * problem. Each strategy needs to set the target width property to a sensible value.
 */
public enum WidthApproximationStrategy implements ILayoutPhaseFactory<RectPackingLayoutPhases, ElkNode> {
    GREEDY,
    TARGET_WIDTH;

    /* (non-Javadoc)
     * @see org.eclipse.elk.core.alg.ILayoutPhaseFactory#create()
     */
    @Override
    public ILayoutPhase<RectPackingLayoutPhases, ElkNode> create() {
        switch (this) {
        case GREEDY:
            return new GreedyWidthApproximator();
            
        case TARGET_WIDTH:
            return new TargetWidthWidthApproximator();
            
        default:
            throw new IllegalArgumentException(
                    "No implementation is available for the width approximator " + this.toString());
        }
    }

}
