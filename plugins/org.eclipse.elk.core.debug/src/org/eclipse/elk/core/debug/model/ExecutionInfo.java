/*******************************************************************************
 * Copyright (c) 2016, 2019 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.core.debug.model;

import java.nio.file.Path;
import java.time.Instant;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.eclipse.elk.core.util.IElkProgressMonitor;
import org.eclipse.elk.core.util.LoggedGraph;

import com.google.common.collect.Lists;

/**
 * Keeps information about a layout execution. The information are derived from {@link IElkProgressMonitor} instances.
 * Since progress monitors can form a hierarchy, so can executions. Once constructed based on a progress monitor, an
 * execution is immutable.
 */
public final class ExecutionInfo {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Fields

    /** The name of this execution. */
    private String name;
    /** When this info object was created. */
    private Instant creationTime;
    /** Our parent execution, if any. */
    private ExecutionInfo parent;
    /** Unmodifiable list of child executions. */
    private List<ExecutionInfo> children;

    /** If this execution represents a graph loaded from a file, this is the file name. */
    private String fileName;
    /** If layout was performed on the graph after it was loaded. */
    private boolean laidOutAfterLoad;

    /** Whether execution times were measured. */
    private boolean executionTimeMeasured;
    /** Amount of time spent in this execution and its child executions. */
    private double executionTimeIncludingChildren;
    /** Amount of time spent in this execution locally, without its child executions. */
    private double executionTimeLocal;
    
    /** Folder the monitor would save debugging information to. */
    private Path debugFolder;

    /** List of log messages. */
    private List<String> logMessages;
    /** Whether there are any descendants that have a non-empty list of log messages. */
    private boolean hasDescendantsWithLogMessages;

    /** List of logged graphs. */
    private List<LoggedGraph> logGraphs;
    /** Whether there are any descendants that have a non-empty list of logged graphs. */
    private boolean hasDescendantsWithLogGraphs;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Construction

    /**
     * Private constructor.
     */
    private ExecutionInfo() {
        creationTime = Instant.now();
    }

    /**
     * Create an Execution from a given {@link IElkProgressMonitor}.
     * 
     * @param monitor
     *            that holds information about execution time.
     * @return the execution with its children.
     */
    public static ExecutionInfo fromProgressMonitor(final IElkProgressMonitor monitor) {
        return fromProgressMonitor(monitor, null);
    }

    /**
     * Create an execution from a given {@link IElkProgressMonitor} and set the file name attribute. Used by the layout
     * graph view to indicate executions loaded through the load graph action.
     * 
     * @param monitor
     *            monitor that holds information about execution time.
     * @param fileName
     *            absolute path to the file this execution is based on.
     * @return the execution with its children.
     */
    public static ExecutionInfo fromProgressMonitorAndFile(final IElkProgressMonitor monitor, String fileName,
            boolean laidOutAfterLoad) {
        
        ExecutionInfo info = fromProgressMonitor(monitor);
        info.fileName = fileName;
        info.laidOutAfterLoad = laidOutAfterLoad;
        return info;
    }

    /**
     * Create a child execution for the given parent.
     * 
     * @param monitor
     *            that holds information about execution time.
     * @param parent
     *            for the child created from the given monitor.
     * @return a child execution with its children.
     */
    private static ExecutionInfo fromProgressMonitor(final IElkProgressMonitor monitor, final ExecutionInfo parent) {
        ExecutionInfo execution = new ExecutionInfo();

        // Basic properties
        execution.name = monitor.getTaskName() != null ? monitor.getTaskName() : "Unnamed";
        execution.executionTimeMeasured = monitor.isExecutionTimeMeasured();
        execution.executionTimeIncludingChildren = monitor.getExecutionTime();
        execution.debugFolder = monitor.getDebugFolder();
        execution.parent = parent;

        // Log messages
        List<String> log = monitor.getLogs();
        if (log == null) {
            log = Collections.emptyList();
        }
        execution.logMessages = Collections.unmodifiableList(Lists.newArrayList(log));

        // Log graphs
        List<LoggedGraph> graphs = monitor.getLoggedGraphs();
        if (graphs == null) {
            graphs = Collections.emptyList();
        }
        execution.logGraphs = Collections.unmodifiableList(Lists.newArrayList(graphs));

        // We compute the local execution time by starting with the total execution time and subtracting child
        // execution times as we iterate over them
        execution.executionTimeLocal = execution.executionTimeIncludingChildren;

        // Add executions for child monitors
        List<ExecutionInfo> childExecutions = new ArrayList<>(monitor.getSubMonitors().size());
        for (IElkProgressMonitor child : monitor.getSubMonitors()) {
            ExecutionInfo childExecution = ExecutionInfo.fromProgressMonitor(child, execution);
            childExecutions.add(childExecution);

            execution.executionTimeLocal -= child.getExecutionTime();
            execution.hasDescendantsWithLogMessages |=
                    childExecution.hasLogMessages() || childExecution.hasDescendantsWithLogMessages();
            execution.hasDescendantsWithLogGraphs |=
                    childExecution.hasLoggedGraphs() || childExecution.hasDescendantsWithLoggedGraphs();
        }
        execution.children = childExecutions;

        // Ensure the local execution time does not drop below zero due to double arithmetic
        execution.executionTimeLocal = Math.max(execution.executionTimeLocal, 0.0);

        return execution;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Accessors

    /**
     * Returns the name of this execution that should be displayed.
     */
    public String getName() {
        return name;
    }

    /**
     * Returns the instant when this object was created.
     */
    public Instant getCreationTime() {
        return creationTime;
    }

    /**
     * Returns the parent execution.
     * 
     * @return the parent, or {@code null} if there is none.
     */
    public ExecutionInfo getParent() {
        return parent;
    }

    /**
     * Returns the children of this execution. Can be modified.
     */
    public List<ExecutionInfo> getChildren() {
        return children;
    }
    
    /**
     * Returns {@code true} if this execution was based on a graph loaded from a file. If this method returns
     * {@code true}, {@link #getFileName()} will return a non-{@code null} value.
     */
    public boolean isLoadedFromFile() {
        return fileName != null;
    }
    
    /**
     * Returns the path to the file this execution was based on or {@code null}.
     */
    public String getFileName() {
        return fileName;
    }
    
    /**
     * Returns whether, if this execution was based on a file, the graph loaded from that file was laid out after having
     * been loaded or not.
     */
    public boolean isLaidOutAfterLoad() {
        return laidOutAfterLoad;
    }

    /**
     * Returns whether execution times were measured. That is, if calls to {@link #getExecutionTimeLocal()} and
     * {@link #getExecutionTimeIncludingChildren()} return anything meaningful.
     */
    public boolean isExecutionTimeMeasured() {
        return executionTimeMeasured;
    }

    /**
     * Returns how much time this execution took to execute, including its child executions.
     */
    public double getExecutionTimeIncludingChildren() {
        return executionTimeIncludingChildren;
    }

    /**
     * Returns how much time this execution took to execute locally, without child executions.
     */
    public double getExecutionTimeLocal() {
        return executionTimeLocal;
    }
    
    /**
     * Returns the path the execution's persisted debug files would be found in.
     * 
     * @return debug folder path.
     */
    public Path getDebugFolder() {
        return debugFolder;
    }

    /**
     * Returns the list of log messages.
     * 
     * @return the list of log messages, as an unmodifiable list.
     */
    public List<String> getLogMessages() {
        return logMessages;
    }

    /**
     * Returns whether or not this execution info contains at least one log message.
     */
    public boolean hasLogMessages() {
        return !logMessages.isEmpty();
    }

    /**
     * Returns whether or not any descendants of this execution would return {@code true} upon a call to
     * {@link #hasLogMessages()}.
     */
    public boolean hasDescendantsWithLogMessages() {
        return hasDescendantsWithLogMessages;
    }

    /**
     * Returns the list of logged graphs.
     * 
     * @return the list of logged graphs, as an unmodifiable list.
     */
    public List<LoggedGraph> getLoggedGraphs() {
        return logGraphs;
    }

    /**
     * Returns whether or not this execution info contains at least one logged graph.
     */
    public boolean hasLoggedGraphs() {
        return !logGraphs.isEmpty();
    }

    /**
     * Returns whether or not any descendants of this execution would return {@code true} upon a call to
     * {@link #hasLoggedGraphs()}.
     */
    public boolean hasDescendantsWithLoggedGraphs() {
        return hasDescendantsWithLogGraphs;
    }

}
