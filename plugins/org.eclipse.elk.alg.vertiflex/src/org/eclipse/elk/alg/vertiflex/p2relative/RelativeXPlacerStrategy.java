/*******************************************************************************
 * Copyright (c) 2023 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 * 
 * SPDX-License-Identifier: EPL-2.0 
 *******************************************************************************/
package org.eclipse.elk.alg.vertiflex.p2relative;

import org.eclipse.elk.alg.vertiflex.VertiFlexLayoutPhases;
import org.eclipse.elk.core.alg.ILayoutPhase;
import org.eclipse.elk.core.alg.ILayoutPhaseFactory;
import org.eclipse.elk.graph.ElkNode;

/**
 * Horizontal node placement strategies.
 *
 */
public enum RelativeXPlacerStrategy implements ILayoutPhaseFactory<VertiFlexLayoutPhases, ElkNode> {

    /**
     * Simple strategy for setting the horizontal positions of nodes. These positions are relative to their parents
     * and chosen such that overlaps are avoided.
     */
    SIMPLE_X_PLACING;

    @Override
    public ILayoutPhase<VertiFlexLayoutPhases, ElkNode> create() {
        switch (this) {
        case SIMPLE_X_PLACING:
            return new RelativeXPlacer();

        default:
            throw new IllegalArgumentException(
                    "No implementation is available for the node placer " + this.toString());
        }
    }

}