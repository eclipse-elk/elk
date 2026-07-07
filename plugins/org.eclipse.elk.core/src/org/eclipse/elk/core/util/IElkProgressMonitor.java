/*******************************************************************************
 * Copyright (c) 2009, 2015 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.core.util;

// elkjs-exclude-start
import java.nio.file.Path;
// elkjs-exclude-end
import java.util.List;

import org.eclipse.elk.graph.ElkNode;

/**
 * Interface for monitors of progress of a job. Besides monitoring the progress of operations (and allowing clients to
 * cancel them), progress monitors also have debugging capabilities by providing execution time measurement as well as
 * ways to log messages and graphs.
 */
public interface IElkProgressMonitor extends IElkCancelIndicator {

    /** constant indicating an unknown amount of work. */
    float UNKNOWN_WORK = -1;
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Work

    /**
     * Notifies that the task will begin after this method has been called.
     * 
     * @param name readable name of the new task
     * @param totalWork total amount of work units, or <code>UNKNOWN_WORK</code>
     *            if this is not specified
     * @return true if the task has not begun before, false otherwise
     */
    boolean begin(String name, float totalWork);

    /**
     * Notifies that the given number of work units has been completed. This
     * method will have no effect if the monitor is closed.
     * 
     * @param work number of work units
     */
    void worked(float work);

    /**
     * Notifies that the current task is done and closes the monitor. This method may be called
     * multiple times after the task has begun, without any effect after the first time.
     */
    void done();
    
    /**
     * Returns true if the task has already begun and is not done yet.
     * 
     * @return true if the task is running
     */
    boolean isRunning();

    /**
     * Returns the name of the task associated with this progress monitor.
     * 
     * @return task name
     */
    String getTaskName();
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Sub-Tasks

    /**
     * Creates a new sub-task that will complete the given amount of work units
     * when it is done. The sub-task begins when {@link #begin(String, int)} is
     * called for the new progress monitor instance, and it ends when
     * {@link #done()} is called for that instance.
     * 
     * @param work number of work units that are completed in the current
     *            monitor instance when the sub-task is done
     * @return a progress monitor for the new sub-task, or null if the monitor
     *         is closed
     */
    IElkProgressMonitor subTask(float work);

    /**
     * Returns a list of all monitors associated with direct sub-tasks.
     * 
     * @return list of sub-task monitors
     */
    List<IElkProgressMonitor> getSubMonitors();

    /**
     * Returns the parent monitor. The parent monitor is the one for which a
     * call to {@link #subTask(int)} resulted in the current monitor instance.
     * 
     * @return the parent monitor, or null if there is none
     */
    IElkProgressMonitor getParentMonitor();
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Debugging
    
    /**
     * Returns whether the logging methods will log anything. If not, clients should not bother producing debug output
     * in the first place.
     * 
     * @return {@code true} if logging methods will actually log stuff.
     */
    boolean isLoggingEnabled();
    
    /**
     * Whether the progress monitor will save logged data on the file system. If so, things are saved only if logging
     * itself is enabled as well (see {@link #isLoggingEnabled()}).
     * 
     * @return {@code true} if logged objects are saved to the file system.
     */
    boolean isLogPersistenceEnabled();
    
    /**
     * Logs the given object.
     * 
     * @param object to be logged.
     */
    void log(Object object);

    /**
     * Returns the collected logs for the task associated with this monitor.
     * 
     * @return list of logs
     */
    List<String> getLogs();
    
    /**
     * Saves the given ElkGraph together with its tag.
     * 
     * @param graph the graph to be logged.
     * @param tag for identifying the graph.
     */
    void logGraph(ElkNode graph, String tag);
    
    /**
     * Saves the given graph together with its tag and type.
     * 
     * @param graph the graph to be logged.
     * @param tag for identifying the graph.
     * @param graphType of the given graph.
     * @throws ClassCastException if {@code graph} does not conform to the expected type.
     */
    void logGraph(Object graph, String tag, LoggedGraph.Type graphType);
    
    /**
     * Returns the collected intermediate graphs for the task associated with this monitor.
     * 
     * @return list of graphs
     */
    List<LoggedGraph> getLoggedGraphs();
    
    /**
     * Returns the path where this monitor would put debug files. Clients can put custom debug files there as well. The
     * path is {@code null} if {@link #isLoggingEnabled()} or {@link #isLogPersistenceEnabled()} return {@code false}.
     * If the path is non-{@code null}, it is not required to exist yet.
     * 
     * @return debug path.
     */
    // elkjs-exclude-start
    Path getDebugFolder();
    // elkjs-exclude-end
    
    /**
     * Returns whether this monitor measures execution time. If it is, the value returned by {@link #getExecutionTime()}
     * must be meaningful.
     * 
     * @return {@code true} if execution time is measured.
     */
    boolean isExecutionTimeMeasured();

    /**
     * Returns the measured execution time for the task associated with this monitor.
     * This is optional: implementations may just return 0 in order to avoid the additional
     * overhead of measuring the execution time.
     * 
     * @return number of seconds used for execution
     */
    double getExecutionTime();

}
