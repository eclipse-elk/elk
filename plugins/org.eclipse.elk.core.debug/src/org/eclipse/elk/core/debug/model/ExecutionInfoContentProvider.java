/*******************************************************************************
 * Copyright (c) 2019 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.core.debug.model;

import java.util.List;
import java.util.function.Predicate;

import org.eclipse.jface.viewers.ITreeContentProvider;

/**
 * Content provider for tree viewers that wish to display execution infos. The input is always assumed to be a list of
 * {@link ExecutionInfo} objects. An optional filter can be installed which determines for a given execution info
 * whether it should be included in the tree viewer or not.
 */
public class ExecutionInfoContentProvider implements ITreeContentProvider {

    /** A filter that always returns {@code true}. */
    private static final Predicate<ExecutionInfo> ALL_PASS_FILTER = (info) -> true;

    /** Our currently active filter. If we have none, it always returns true. */
    private Predicate<ExecutionInfo> filter = ALL_PASS_FILTER;

    /**
     * Determines whether this content provider filters out elements that are not on a path towards elements that have
     * log entries.
     */
    public boolean isFiltering() {
        return filter == ALL_PASS_FILTER;
    }

    /**
     * The filter to be used to filter execution infos.
     * 
     * @param filter
     *            the filter to be used or {@code null} to deactivate filtering.
     */
    public void setFilter(Predicate<ExecutionInfo> filter) {
        if (filter == null) {
            this.filter = ALL_PASS_FILTER;
        } else {
            this.filter = filter;
        }
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Content Provider

    @Override
    public Object[] getElements(Object inputElement) {
        if (inputElement instanceof List<?>) {
            return ((List<?>) inputElement).stream()
                    .map(o -> (ExecutionInfo) o)
                    .filter(filter)
                    .toArray();
        } else {
            return new Object[0];
        }
    }

    @Override
    public Object[] getChildren(Object parentElement) {
        if (parentElement instanceof ExecutionInfo) {
            // If our filter is active, filter out those that are not on a path towards log messages
            return ((ExecutionInfo) parentElement).getChildren().stream()
                    .filter(filter)
                    .toArray();

        } else {
            return new ExecutionInfo[0];
        }
    }

    @Override
    public Object getParent(Object element) {
        if (element instanceof ExecutionInfo) {
            return ((ExecutionInfo) element).getParent();
        }
        return null;
    }

    @Override
    public boolean hasChildren(Object element) {
        if (element instanceof ExecutionInfo) {
            // If our filter is turned on, this thing only has children if at least one of them has log items or has
            // ancestors with log items
            return ((ExecutionInfo) element).getChildren().stream().anyMatch(filter);

        } else {
            return false;
        }
    }

}
