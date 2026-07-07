/*******************************************************************************
 * Copyright (c) 2022 - 2023 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 * 
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later  
 *******************************************************************************/
package org.eclipse.elk.alg.topdownpacking;

import org.eclipse.elk.core.alg.ILayoutPhase;
import org.eclipse.elk.core.alg.ILayoutPhaseFactory;

/**
 * During the white space elimination phase nodes are expanded beyond their initial size
 * to completely fill the space provided by their parent.
 *
 */
public enum WhitespaceEliminationStrategy
        implements ILayoutPhaseFactory<TopdownPackingPhases, GridElkNode> {
    /**
     * Equally expands nodes on the bottom row to fill up remaining space.
     */
    BOTTOM_ROW_EQUAL_WHITESPACE_ELIMINATOR;

    /**
     * {@inheritDoc}
     */
    @Override
    public ILayoutPhase<TopdownPackingPhases, GridElkNode> create() {
        switch (this) {
        case BOTTOM_ROW_EQUAL_WHITESPACE_ELIMINATOR:
            return new BottomRowEqualWhitespaceEliminator();
        default:
            throw new IllegalArgumentException(
                    "No implementation is available for the whitespace eliminator" + this.toString());
        }
    }

}
