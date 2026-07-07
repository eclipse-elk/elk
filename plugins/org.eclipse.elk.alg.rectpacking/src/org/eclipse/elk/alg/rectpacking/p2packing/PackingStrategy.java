/*******************************************************************************
 * Copyright (c) 2022 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 * 
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later 
 *******************************************************************************/
package org.eclipse.elk.alg.rectpacking.p2packing;

import org.eclipse.elk.alg.rectpacking.RectPackingLayoutPhases;
import org.eclipse.elk.core.alg.ILayoutPhase;
import org.eclipse.elk.core.alg.ILayoutPhaseFactory;
import org.eclipse.elk.graph.ElkNode;

/**
 * Packing strategy and factory for strip packing algorithms.
 */
public enum PackingStrategy implements ILayoutPhaseFactory<RectPackingLayoutPhases, ElkNode> {
    COMPACTION,
    SIMPLE,
    NONE;
    
    /* (non-Javadoc)
     * @see org.eclipse.elk.core.alg.ILayoutPhaseFactory#create()
     */
    @Override
    public ILayoutPhase<RectPackingLayoutPhases, ElkNode> create() {
        switch (this) {
        case COMPACTION:
            return new Compactor();
        case SIMPLE:
            return new SimplePlacement();
        case NONE:
            return new NoPlacement();
        default:
            return null;
        }
    }

}
