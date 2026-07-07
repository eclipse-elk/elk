/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.spore.options;

/**
 * An enumeration of methods for selecting a root node when constructing a spanning tree.
 */
public enum RootSelection {
    /** A Node is manually selected by label as root. */
    FIXED,
    /** The node most central to the drawing is selected as root. */
    CENTER_NODE;
}
