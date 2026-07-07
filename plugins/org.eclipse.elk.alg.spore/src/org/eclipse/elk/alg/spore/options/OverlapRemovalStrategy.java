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
import org.eclipse.elk.alg.spore.p3execution.GrowTreePhase;
import org.eclipse.elk.core.alg.ILayoutPhase;
import org.eclipse.elk.core.alg.ILayoutPhaseFactory;

/**
 * Definition of the strategy for the third phase of overlap removal by growing a tree.
 */
public enum OverlapRemovalStrategy implements ILayoutPhaseFactory<SPOrEPhases, Graph> {
    /** Overlap removal by growing a tree. */
    GROW_TREE;
    
    @Override
    public ILayoutPhase<SPOrEPhases, Graph> create() {
        return new GrowTreePhase();
    }
}
