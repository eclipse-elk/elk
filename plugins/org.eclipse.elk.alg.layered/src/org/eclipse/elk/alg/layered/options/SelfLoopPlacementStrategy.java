/*******************************************************************************
 * Copyright (c) 2014, 2015 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.layered.options;

/**
 * Enumeration to distinguish the possible ways to distribute the selfLoops around a node.
 * 
 * @author tit
 */
public enum SelfLoopPlacementStrategy {
    /** Distributes the loops equally around the node. */
    EQUALLY_DISTRIBUTED,
    /** Stacks all loops to the north side of the node. */
    NORTH_STACKED,
    /** Loops are placed sequentially (next to each other) to the north side of the node. */
    NORTH_SEQUENCE;
}
