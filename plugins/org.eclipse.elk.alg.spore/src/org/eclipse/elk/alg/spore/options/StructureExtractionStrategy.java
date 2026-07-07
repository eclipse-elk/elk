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
import org.eclipse.elk.alg.spore.p1structure.DelaunayTriangulationPhase;
import org.eclipse.elk.core.alg.ILayoutPhase;
import org.eclipse.elk.core.alg.ILayoutPhaseFactory;

/**
 * Definition of the structure extraction strategy for SPOrE.
 */
public enum StructureExtractionStrategy implements ILayoutPhaseFactory<SPOrEPhases, Graph> {
    /** A Delaunay triangulation. */
    DELAUNAY_TRIANGULATION;

    @Override
    public ILayoutPhase<SPOrEPhases, Graph> create() {
        switch (this) {
        case DELAUNAY_TRIANGULATION:
            return new DelaunayTriangulationPhase();

        default:
            throw new IllegalArgumentException(
                    "No implementation available for " + this.toString());
        }
    }

}
