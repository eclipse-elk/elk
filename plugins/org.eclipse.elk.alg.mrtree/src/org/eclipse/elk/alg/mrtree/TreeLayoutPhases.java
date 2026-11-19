/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
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
package org.eclipse.elk.alg.mrtree;

import org.eclipse.elk.alg.mrtree.graph.TGraph;
import org.eclipse.elk.alg.mrtree.p1treeify.DFSTreeifyer;
import org.eclipse.elk.alg.mrtree.p2order.NodeOrderer;
import org.eclipse.elk.alg.mrtree.p3place.NodePlacer;
import org.eclipse.elk.alg.mrtree.p4route.EdgeRouter;
import org.eclipse.elk.core.alg.ILayoutPhase;
import org.eclipse.elk.core.alg.ILayoutPhaseFactory;

/**
 * This algorithm's layout phases. Since none of the phases have alternative implementations, this enumeration also
 * serves as a layout phase factory (this is not quite true for the node ordering phase, but there was never a way
 * to switch between the two implementations; if there is one some day, each phase will have to get its own factory).
 */
public enum TreeLayoutPhases implements ILayoutPhaseFactory<TreeLayoutPhases, TGraph> {
    
    /** Phase 1. */
    P1_TREEIFICATION,
    /** Phase 2. */
    P2_NODE_ORDERING,
    /** Phase 3. */
    P3_NODE_PLACEMENT,
    /** Phase 4. */
    P4_EDGE_ROUTING;

    /* (non-Javadoc)
     * @see org.eclipse.elk.core.alg.ILayoutPhaseFactory#create()
     */
    @Override
    public ILayoutPhase<TreeLayoutPhases, TGraph> create() {
        switch (this) {
        case P1_TREEIFICATION:
            return new DFSTreeifyer();
            
        case P2_NODE_ORDERING:
            return new NodeOrderer();
            
        case P3_NODE_PLACEMENT:
            return new NodePlacer();
            
        case P4_EDGE_ROUTING:
            return new EdgeRouter();
            
        default:
            throw new IllegalArgumentException(
                    "No implementation is available for the layout phase " + toString());
        }
    }
    
}
