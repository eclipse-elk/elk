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
package org.eclipse.elk.alg.spore.options;

import org.eclipse.elk.alg.common.TEdge;

/**
 * An enumeration of possible cost function that may be used for the construction of a spanning tree.
 */
public enum SpanningTreeCostFunction {
    /** The euclidean distance between the center points of the two nodes connected by a {@link TEdge}. */
    CENTER_DISTANCE,
    /** Distance between the circumcircles of the two graph elements pertaining to a {@link TEdge}. */
    CIRCLE_UNDERLAP,
    /** Underlap between the rectangles of the two graph elements pertaining to a {@link TEdge}. 
     * Restrict the movement of the rectangles to the line connecting their center points, then the underlap is the 
     * distance they could be moved towards each other without colliding.*/
    RECTANGLE_UNDERLAP,
    /** Cost function for overlap removal by Nachmanson et al. "Node Overlap Removal by Growing a Tree".
     *  If the nodes do overlap, -overlap is returned otherwise the shortest distance. */
    INVERTED_OVERLAP,
    /** The minimum of the euclidean distances between the nodes of an edge and the root of the spanning tree. */
    MINIMUM_ROOT_DISTANCE; 
}
