/*******************************************************************************
 * Copyright (c) 2022 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 * 
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later 
 *******************************************************************************/
package org.eclipse.elk.alg.rectpacking.p3whitespaceelimination;

import org.eclipse.elk.alg.rectpacking.RectPackingLayoutPhases;
import org.eclipse.elk.core.alg.ILayoutPhase;
import org.eclipse.elk.core.alg.ILayoutPhaseFactory;
import org.eclipse.elk.graph.ElkNode;

/**
 * Strategy factory for whitespace elimination.
 *
 */
public enum WhiteSpaceEliminationStrategy implements ILayoutPhaseFactory<RectPackingLayoutPhases, ElkNode> {
    /**
     * The packing that is divided into rows, stacks, blocks, and subrows equally divides available space between this
     * structures by moving the structures and the rectangles before expanding the rectangles such that all gaps are
     * filled.
     * The size of the parent remains unchanged.
     */
    EQUAL_BETWEEN_STRUCTURES,
    /**
     * Same as the equal between structures but the strategy expands to the desired aspect ratio.
     * Hence, this increases the size of the parent.
     */
    TO_ASPECT_RATIO,
    /**
     * Explicitly set that no whitespace elimination is done.
     * This can be done explicitly, since some language do not support the use of null.
     */
    NONE;

    /* (non-Javadoc)
     * @see org.eclipse.elk.core.alg.ILayoutPhaseFactory#create()
     */
    @Override
    public ILayoutPhase<RectPackingLayoutPhases, ElkNode> create() {
        switch (this) {
        case EQUAL_BETWEEN_STRUCTURES:
            return new EqualWhitespaceEliminator();
        case TO_ASPECT_RATIO:
            return new ToAspectratioNodeExpander();
        case NONE:
        default:
            return null;
        }
    }

}
