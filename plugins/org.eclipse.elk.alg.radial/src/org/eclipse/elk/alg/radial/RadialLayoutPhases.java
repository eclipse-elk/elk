/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.radial;

import org.eclipse.elk.alg.radial.p1position.EadesRadial;
import org.eclipse.elk.alg.radial.p2routing.StraightLineEdgeRouter;
import org.eclipse.elk.core.alg.ILayoutPhase;
import org.eclipse.elk.core.alg.ILayoutPhaseFactory;
import org.eclipse.elk.graph.ElkNode;

/**
 * The phases of the radial layout provider.
 *
 */
public enum RadialLayoutPhases implements ILayoutPhaseFactory<RadialLayoutPhases, ElkNode> {

    /** Place the nodes. */
    P1_NODE_PLACEMENT,
    /** Route the edges of the graph. */
    P2_EDGE_ROUTING;

   @Override
    public ILayoutPhase<RadialLayoutPhases, ElkNode> create() {
        switch (this) {
        case P1_NODE_PLACEMENT:
            return new EadesRadial();
        case P2_EDGE_ROUTING:
            return new StraightLineEdgeRouter();
        default:
            throw new IllegalArgumentException(
                    "No implementation is available for the layout processor " + this.toString());
        }
    }

}