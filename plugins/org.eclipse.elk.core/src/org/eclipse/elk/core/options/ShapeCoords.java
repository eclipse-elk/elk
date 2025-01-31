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
 * Shape coordinate systems for JSON output. To be accessed using {@link CoreOptions#JSON_SHAPE_COORDS}.
 * Applies to nodes, ports, and labels of nodes and ports.
 */
public enum ShapeCoords {

    /**
     * Inherit the parent shape's coordinate system. The root node has no parent shape; here, this setting defaults to
     * {@link #PARENT}.
     */
    INHERIT,
    /** relative to the shape's JSON parent shape. */
    PARENT,
    /** relative to the root node, a.k.a. global coordinates. */
    ROOT;
    
}
