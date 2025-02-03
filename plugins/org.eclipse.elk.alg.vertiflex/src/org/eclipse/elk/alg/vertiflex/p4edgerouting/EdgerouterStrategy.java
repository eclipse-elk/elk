/*******************************************************************************
 * Copyright (c) 2023 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 * 
 * SPDX-License-Identifier: EPL-2.0 
 *******************************************************************************/
package org.eclipse.elk.alg.vertiflex.p4edgerouting;

import org.eclipse.elk.alg.vertiflex.VertiFlexLayoutPhases;
import org.eclipse.elk.core.alg.ILayoutPhase;
import org.eclipse.elk.core.alg.ILayoutPhaseFactory;
import org.eclipse.elk.graph.ElkNode;

/**
 * Edge routing strategies.
 *
 */
public enum EdgerouterStrategy implements ILayoutPhaseFactory<VertiFlexLayoutPhases, ElkNode> {

    /**
     * Straight edge routing.
     */
    DIRECT_ROUTING,
    /**
     * Routing with bendpoints if nodes are in the way.
     */
    BEND_ROUTING;

    @Override
    public ILayoutPhase<VertiFlexLayoutPhases, ElkNode> create() {
        switch (this) {
        case DIRECT_ROUTING:
            return new StraightEdgeRouter();
        case BEND_ROUTING:
            return new BendEdgeRouter();
        default:
            throw new IllegalArgumentException(
                    "No implementation is available for the edge router " + this.toString());
        }
    }

}