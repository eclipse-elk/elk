/*******************************************************************************
 * Copyright (c) 2016, 2020 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.core.options;

/**
 * Options for setting how children of nodes should be handled in the current layout run.
 * 
 * <p>
 * The basic idea is this: If you want nodes to be laid out together across hierarchy levels, set hierarchy handling to
 * {@link #INCLUDE_CHILDREN} on all nodes that should be laid out in one go. As soon as a node's content should be laid
 * out in a separate layout run, set the node's hierarchy handling to {@link #SEPARATE_CHILDREN}.
 * </p>
 * 
 * <p>
 * If the layout algorithm doesn't support hierarchical layout, this property is ignored and the layout is calculated
 * separately for each child hierarchy.
 * </p>
 * <p>
 * <i>Note:</i> Layout algorithms only need to differentiate between {@link #INCLUDE_CHILDREN} and
 * {@link #SEPARATE_CHILDREN} as {@link #INHERIT} is evaluated and set to the appropriate more specific value by ELK.
 * </p>
 */
public enum HierarchyHandling {

    /**
     * Inherit the parent node's hierarchy handling. The root node has no parent node; here, this setting defaults to
     * {@link #SEPARATE_CHILDREN}.
     */
    INHERIT,

    /**
     * Allows the node's children to be included in the current layout run. Which children are included in the layout
     * run is determined by their hierarchy handling setting. For a child to actually be included, its hierarchy
     * handling must be set to either {@link #INHERIT} or {@link #INCLUDE_CHILDREN}.
     */
    INCLUDE_CHILDREN,

    /**
     * Lays out the node with a new layout run. Even if its parent node is set to {@link #INCLUDE_CHILDREN}, this node
     * will trigger a separate layout run and will thus not be included in the parent node's layout run.
     */
    SEPARATE_CHILDREN

}
