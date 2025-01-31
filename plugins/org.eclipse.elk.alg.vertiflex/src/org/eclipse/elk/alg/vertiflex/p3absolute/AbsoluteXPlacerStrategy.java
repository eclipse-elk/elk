/*******************************************************************************
 * Copyright (c) 2023 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 * 
 * SPDX-License-Identifier: EPL-2.0 
 *******************************************************************************/
package org.eclipse.elk.alg.vertiflex.p3absolute;

import org.eclipse.elk.alg.vertiflex.VertiFlexLayoutPhases;
import org.eclipse.elk.core.alg.ILayoutPhase;
import org.eclipse.elk.core.alg.ILayoutPhaseFactory;
import org.eclipse.elk.graph.ElkNode;

/**
 * Strategies for absolute node placement.
 *
 */
public enum AbsoluteXPlacerStrategy implements ILayoutPhaseFactory<VertiFlexLayoutPhases, ElkNode> {

    /**
     * Compute absolute x-coordinates based on relative coordinates
     */
    ABSOLUTE_XPLACING;

    @Override
    public ILayoutPhase<VertiFlexLayoutPhases, ElkNode> create() {
        switch (this) {
        case ABSOLUTE_XPLACING:
            return new AbsoluteXPlacer();

        default:
            throw new IllegalArgumentException(
                    "No implementation is available for the node placer " + this.toString());
        }
    }

}