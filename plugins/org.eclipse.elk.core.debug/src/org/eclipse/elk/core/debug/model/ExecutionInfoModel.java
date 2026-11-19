/*******************************************************************************
 * Copyright (c) 2019 Kiel University and others.
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
package org.eclipse.elk.core.debug.model;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * Central repository of execution infos that store information about past layout runs to be displayed in debug views.
 */
public class ExecutionInfoModel {

    /** List of execution infos. */
    private final List<ExecutionInfo> executions = new ArrayList<>();
    /** List of event listeners. */
    private final Set<IExecutionInfoModelListener> listeners = new HashSet<>();

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Model Access

    /**
     * Returns an unmodifiable view on our list of execution infos.
     */
    public List<ExecutionInfo> getExecutionInfos() {
        return Collections.unmodifiableList(executions);
    }

    /**
     * Adds the information contained in the given progress monitor to our list of executions and notifies any
     * listeners.
     */
    public void addExecution(final ExecutionInfo info) {
        executions.add(info);
        notifyListeners();
    }

    /**
     * Replaces an existing execution by a new one. Used to reload graphs loaded from file.
     * 
     * @param oldInfo
     *            old execution to replace.
     * @param newInfo
     *            new execution.
     */
    public void replaceExecution(final ExecutionInfo oldInfo, final ExecutionInfo newInfo) {
        int index = executions.indexOf(oldInfo);
        if (index >= 0) {
            executions.set(index, newInfo);
            notifyListeners();
        }
    }

    /**
     * Removes the execution infos as well as their subtrees and notifies any listeners.
     * 
     * @param infos
     *            the execution infos to remove.
     */
    public void removeExecutionInfos(final ExecutionInfo... infos) {
        for (ExecutionInfo info : infos) {
            if (info.getParent() == null) {
                // This thing is in our top-level list
                executions.remove(info);
            } else {
                // This is not a top-level element, so simply remove it from the parent
                info.getParent().getChildren().remove(info);
            }
        }

        notifyListeners();
    }

    /**
     * Removes all execution infos and notifies any listeners.
     */
    public void removeAllExecutionInfos() {
        executions.clear();
        notifyListeners();
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Events

    /**
     * Adds the given listener to our set of event listeners.
     */
    public void addExecutionInfoModelListener(final IExecutionInfoModelListener listener) {
        listeners.add(listener);
    }

    /**
     * Removes the given listener from our set of event listeners.
     */
    public void removeExecutionInfoModelListener(final IExecutionInfoModelListener listener) {
        listeners.remove(listener);
    }

    /**
     * Notifies listeners of a change in our execution info tree.
     */
    private void notifyListeners() {
        for (IExecutionInfoModelListener listener : listeners) {
            listener.executionInfoChanged();
        }
    }

}
