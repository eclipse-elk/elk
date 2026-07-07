/*******************************************************************************
 * Copyright (c) 2020 Kiel University and others.
 *
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.layered.p5edges.orthogonal.direction;

import org.eclipse.elk.core.options.Direction;

/**
 * Enumeration of available routing directions. This is different from a graph's {@link Direction} because different
 * routing directions may arise in a single graph. Also, we don't have an east-to-west routing direction here.
 */
public enum RoutingDirection {
    
    /** west to east routing direction. */
    WEST_TO_EAST,
    /** north to south routing direction. */
    NORTH_TO_SOUTH,
    /** south to north routing direction. */
    SOUTH_TO_NORTH;
    
}