/*******************************************************************************
 * Copyright (c) 2016 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.layered.options;

/**
 * Definition of strategies for horizontal compaction.
 */
public enum GraphCompactionStrategy {
    /** doesn't apply compaction. */
    NONE,
    /** compacts to the left. */
    LEFT,
    /** compacts to the right. */
    RIGHT,
    /** compacts left, locks {@link CNode}s that are not constrained, compacts right. */
    LEFT_RIGHT_CONSTRAINT_LOCKING,
    /** compacts left, locks {@link CNode}s based on their {@link CompactionLock}, compacts right.
     *  Yields better results for average edge length because {@link CNode}s are locked in the
     *  direction of fewer connections. */
    LEFT_RIGHT_CONNECTION_LOCKING,
    /** compacts as much as possible, however instead of minimizing the width it minimizes edge
     *  length.
     */
    EDGE_LENGTH;
}
