/*******************************************************************************
 * Copyright (c) 2024 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.core.options;

/**
 * Edge coordinate systems for JSON output. To be accessed using {@link CoreOptions#JSON_EDGE_COORDS}.
 * Applies to edges, and to labels of edges.
 */
public enum EdgeCoords {

    /**
     * Inherit the parent node's coordinate system. The root node has no parent node; here, this setting defaults to
     * {@link #CONTAINER}.
     */
    INHERIT,
    /** relative to the edge's proper container node. */
    CONTAINER,
    /** relative to the edge's JSON parent node. */
    PARENT,
    /** relative to the root node, a.k.a. global coordinates. */
    ROOT;
    
}
