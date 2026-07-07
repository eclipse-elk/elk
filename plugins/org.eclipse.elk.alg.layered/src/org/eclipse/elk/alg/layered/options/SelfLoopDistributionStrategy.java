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
 * The possible ways to distribute self loops around a node.
 */
public enum SelfLoopDistributionStrategy {
    
    /** Distributes the loops equally around the node. */
    EQUALLY,
    /** Puts all loops to the north side of the node. */
    NORTH,
    /** Loops are distributed over the north and the south side of the node. */
    NORTH_SOUTH;
    
}
