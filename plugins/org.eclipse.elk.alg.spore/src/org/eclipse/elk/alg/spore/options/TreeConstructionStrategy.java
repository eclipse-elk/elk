/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.spore.options;

import org.eclipse.elk.alg.spore.SPOrEPhases;
import org.eclipse.elk.alg.spore.graph.Graph;
import org.eclipse.elk.alg.spore.p2processingorder.MaxSTPhase;
import org.eclipse.elk.alg.spore.p2processingorder.MinSTPhase;
import org.eclipse.elk.core.alg.ILayoutPhase;
import org.eclipse.elk.core.alg.ILayoutPhaseFactory;

/**
 * Definition of the spanning tree construction strategy for the processing order phase.
 */
public enum TreeConstructionStrategy implements ILayoutPhaseFactory<SPOrEPhases, Graph> {
    /** Minimum spanning tree construction. */
    MINIMUM_SPANNING_TREE,
    /** Maximum spanning tree construction. */
    MAXIMUM_SPANNING_TREE;

    @Override
    public ILayoutPhase<SPOrEPhases, Graph> create() {
        switch (this) {
        case MINIMUM_SPANNING_TREE:
            return new MinSTPhase();

        case MAXIMUM_SPANNING_TREE:
            return new MaxSTPhase();

        default:
            throw new IllegalArgumentException(
                    "No implementation available for " + this.toString());
        }
    }

}
