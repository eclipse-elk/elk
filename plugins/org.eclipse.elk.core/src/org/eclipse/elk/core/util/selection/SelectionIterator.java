/*******************************************************************************
 * Copyright (c) 2016 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.core.util.selection;

import java.util.Set;

import org.eclipse.elk.core.util.ElkUtil;
import org.eclipse.elk.graph.ElkEdge;
import org.eclipse.elk.graph.ElkGraphElement;
import org.eclipse.elk.graph.ElkPort;
import org.eclipse.emf.common.util.AbstractTreeIterator;

/**
 * This class can be used to implement Iterators for selections. The {@link SelectionIterator} is used in
 * {@link ElkUtil} to get the {@link ElkGraphElement graph elements} connected to a selected edge.
 * 
 * @see DefaultSelectionIterator
 */
public abstract class SelectionIterator extends AbstractTreeIterator<ElkGraphElement> {

    private static final long serialVersionUID = 478793714459586388L;

    /**
     * The list of already visited ports. Used to break infinite loops.
     */
    protected Set<ElkPort> visited; // SUPPRESS CHECKSTYLE VisibilityModifier

    /**
     * Creates a {@link SelectionIterator} which needs to be configured afterwards by adding the set
     * of visited nodes. The iterator won't include the starting object.
     * 
     * @param edge
     *            The object to iterate from
     */
    public SelectionIterator(final ElkEdge edge) {
        super(edge, false);
    }

    /**
     * Attach a set of {@link ElkPortt ports} to the iterator to be used as a set of visited ports.
     * Can be used to share a set of nodes across multiple iterators (target and source iterator).
     * 
     * @param visitedSet
     *            The set of nodes to be used as visited set.
     */
    public void attachVisitedSet(final Set<ElkPort> visitedSet) {
        this.visited = visitedSet;
    }
}
