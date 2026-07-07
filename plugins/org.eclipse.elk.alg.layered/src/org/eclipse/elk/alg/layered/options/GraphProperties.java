/*******************************************************************************
 * Copyright (c) 2011, 2015 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.layered.options;

/**
 * An enumeration of properties a graph may have. These can be used as part of an {@code EnumSet} to
 * base decisions on graph properties. For example, self-loop processing may be skipped if the graph
 * doesn't contain self-loops in the first place.
 * 
 * <p>An {@code EnumSet} for this enumeration can be attached to a graph via the
 * {@link LayeredOptions#GRAPH_PROPERTIES} property.</p>
 * 
 * <p>{@link org.eclipse.elk.graph.properties.GraphFeature GraphFeature} serves a similar purpose in ELK.
 * However, since this enumeration is more specific to what ELK Layered supports, the two should stay
 * separate.
 * 
 * @author cds
 */
public enum GraphProperties {
    
    /** The graph contains comment boxes. */
    COMMENTS,
    /** The graph contains dummy nodes representing external ports. */
    EXTERNAL_PORTS,
    /** The graph contains hyperedges. */
    HYPEREDGES,
    /** The graph contains hypernodes (nodes that are marked as such). */
    HYPERNODES,
    /** The graph contains ports that are not free for positioning. */
    NON_FREE_PORTS,
    /** The graph contains ports on the northern or southern side. */
    NORTH_SOUTH_PORTS,
    /** The graph contains self-loops. */
    SELF_LOOPS,
    /** The graph contains node labels. */
    CENTER_LABELS,
    /** The graph contains head or tail edge labels. */
    END_LABELS,
    /** The graph's nodes are partitioned. */
    PARTITIONS;
    
}
