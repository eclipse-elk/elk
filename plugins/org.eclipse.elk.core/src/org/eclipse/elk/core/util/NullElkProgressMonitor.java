/*******************************************************************************
 * Copyright (c) 2019 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.core.util;

import java.nio.file.Path;
import java.util.List;

import org.eclipse.elk.core.util.LoggedGraph.Type;
import org.eclipse.elk.graph.ElkNode;

/**
 * An implementation of {@link IElkProgressMonitor} which does not do anything. The primary purpose of this monitor is
 * for it to be used with unit tests. The only method which returns something sensible is {@link #subTask(float)}, which
 * simply returns the instance it was called on (since the progress monitor doesn't do anything anyway, we don't bother
 * creating a new instance).
 */
public class NullElkProgressMonitor implements IElkProgressMonitor {

    @Override
    public boolean isCanceled() {
        return false;
    }

    @Override
    public boolean begin(final String name, final float totalWork) {
        return false;
    }

    @Override
    public void worked(final float work) {
    }

    @Override
    public void done() {
    }

    @Override
    public boolean isRunning() {
        return false;
    }

    @Override
    public String getTaskName() {
        return null;
    }

    @Override
    public IElkProgressMonitor subTask(final float work) {
        return this;
    }

    @Override
    public List<IElkProgressMonitor> getSubMonitors() {
        return null;
    }

    @Override
    public IElkProgressMonitor getParentMonitor() {
        return null;
    }

    @Override
    public boolean isLoggingEnabled() {
        return false;
    }

    @Override
    public boolean isLogPersistenceEnabled() {
        return false;
    }

    @Override
    public void log(final Object object) {
    }

    @Override
    public List<String> getLogs() {
        return null;
    }

    @Override
    public void logGraph(final ElkNode graph, final String tag) {
    }

    @Override
    public void logGraph(final Object graph, final String tag, final Type graphType) {
    }

    @Override
    public List<LoggedGraph> getLoggedGraphs() {
        return null;
    }

    // elkjs-exclude-start
    @Override
    public Path getDebugFolder() {
        return null;
    }
    // elkjs-exclude-end

    @Override
    public boolean isExecutionTimeMeasured() {
        return false;
    }

    @Override
    public double getExecutionTime() {
        return 0;
    }

}
