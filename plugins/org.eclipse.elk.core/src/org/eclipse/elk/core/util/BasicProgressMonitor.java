/*******************************************************************************
 * Copyright (c) 2009, 2019 Kiel University and others.
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
package org.eclipse.elk.core.util;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;

import org.eclipse.elk.graph.ElkNode;
import org.eclipse.emf.ecore.util.EcoreUtil;

/**
 * Base class for implementations of progress monitors. This class performs execution time measurement, logs debug
 * output, keeps track of the amount of completed work, and handles sub-tasks properly.
 * 
 * <p>Progress monitors can be configured to persist log entries. What this means depends on whether the log entries
 * were logged via a call to {@link #log(Object)} or to the {@code logGraph(...)} methods. All regular log entries end
 * up in a file called {@code log.txt} in the monitor's folder. All graphs are persisted into their own files, with the
 * file names based on their tags.</p>
 * 
 * <p>Each monitor has a separate folder to save logs to. The folders are nested to reflect the monitor structure. All
 * will end up in a folder called {@code logs} in ELK's {@link ElkUtil#debugFolderPath(String...) root debug folder}.
 * The folder name for a monitor will then depend on whether the monitor has a parent or not.</p>
 * 
 * <ul>
 *   <li>If the monitor has a parent, the name will be built according to this template:
 *       {@code <timestamp><two-random-letters>-<name>}.</li>
 *   <li>If it does not have a parent, the name will be built according to this template:
 *       {@code <index>-<name>}.</li>
 * </ul>
 * 
 * <p>In both cases, the name is the monitor's task name with space characters replaced by underscores and everything
 * filtered out which isn't a regular character or a number.</p>
 */
public class BasicProgressMonitor implements IElkProgressMonitor {
    
    /** name of the parent folder for all debug data. */
    public static final String ROOT_DEBUG_FOLDER_NAME = "logs";
    /** indicates an infinite number of hierarchy levels for progress reporting. */
    private static final int INFINITE_HIERARCHY_LEVELS = -1;
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Properties

    /** the parent monitor. */
    private BasicProgressMonitor parentMonitor;
    /** list of child monitors. */
    private final List<IElkProgressMonitor> children = new LinkedList<IElkProgressMonitor>();
    /** the maximal number of hierarchy levels for which progress is reported. */
    private int maxLevels = INFINITE_HIERARCHY_LEVELS;
    
    /** the name of the associated task. */
    private String taskName;
    /** indicates whether the monitor has been closed. */
    private boolean closed = false;
    /** the number of work units that can be completed in total. */
    private float totalWork;
    /** the amount of work that is completed. */
    private float completedWork = 0.0f;
    /** the number of work units that will be consumed after completion of the currently active child task. */
    private float currentChildWork = -1;
    
    /** whether logging should be enabled (if not, logging objects does nothing). */
    private boolean recordLogs = false;
    /** whether logged data are saved on the file system (provided that logging itself is enabled). */
    private boolean persistLogs = false;
    /** list of log messages (only non-null if logs are enabled). */
    private List<String> logMessages = null;
    /** list of logged graphs (only non-null if logs are enabled). */
    private List<LoggedGraph> logGraphs = null;
    // elkjs-exclude-start
    /** path to save debug output in. */
    private Path debugFolder = null;
    /** file to save log messages to. */
    private Path logFile = null;
    // elkjs-exclude-end
    /** whether the execution time shall be measured when the task is done. */
    private boolean recordExecutionTime = false;
    /** the start time of the associated task, in nanoseconds. */
    private long startTime;
    /** the total time of the associated task, in seconds. */
    private double totalTime;

    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Creation

    /**
     * Creates a progress monitor with infinite number of hierarchy levels.
     */
    public BasicProgressMonitor() {
        
    }

    /**
     * Creates a progress monitor with the given maximal number of hierarchy levels. If the
     * number is negative, the hierarchy levels are infinite. Otherwise progress is
     * reported to parent monitors only up to the specified number of levels.
     *
     * @deprecated Use {@link #withMaxHierarchyLevels(int)} instead.
     * @param themaxLevels the maximal number of hierarchy levels for which progress is
     *     reported
     */
    public BasicProgressMonitor(final int themaxLevels) {
        withMaxHierarchyLevels(themaxLevels);
    }

    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Configuration
    
    /**
     * Configures the monitor for the given maximum number of hierarchy levels. If the number is negative, the hierarchy
     * levels are infinite. Otherwise, progress is reported to parent monitors only up to the specified number of
     * levels.This method should be called right after the monitor has been created and allows for method chaining.
     * 
     * @param levels
     *            the maximum number of hierarchy levels for which progress is reported.
     * @return this progress monitor.
     */
    public BasicProgressMonitor withMaxHierarchyLevels(final int levels) {
        if (levels < 0) {
            this.maxLevels = INFINITE_HIERARCHY_LEVELS;
        } else {
            this.maxLevels = levels;
        }
        
        return this;
    }
    
    /**
     * Enable or disable logging on this monitor. This method should be called right after the monitor has been created
     * and allows for method chaining.
     * 
     * @param enabled {@code true} if log messages and graphs should be saved.
     * @return this progress monitor.
     */
    public BasicProgressMonitor withLogging(final boolean enabled) {
        recordLogs = enabled;
        
        if (recordLogs) {
            logMessages = new ArrayList<>();
            logGraphs = new ArrayList<>();
        } else {
            logMessages = null;
            logGraphs = null;
        }
        
        return this;
    }
    
    /**
     * Enable or disable log persistence on this monitor. This method should be called right after the monitor has been
     * created and allows for method chaining.
     * 
     * @param enabled {@code true} if log messages and graphs should be saved to the file system.
     * @return this progress monitor.
     */
    public BasicProgressMonitor withLogPersistence(final boolean enabled) {
        persistLogs = enabled;
        return this;
    }

    /**
     * Enable or disable execution time measurement on this monitor. This method should be called right after the
     * monitor has been created and allows for method chaining.
     * 
     * @param enabled {@code true} if execution times should be measured.
     * @return this progress monitor.
     */
    public BasicProgressMonitor withExecutionTimeMeasurement(final boolean enabled) {
        recordExecutionTime = enabled;
        return this;
    }
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Work

    @Override
    public final boolean begin(final String name, final float thetotalWork) {
        if (closed) {
            throw new IllegalStateException("The task is already done.");
        } else if (this.taskName != null) {
            return false;
        } else {
            if (name == null) {
                throw new NullPointerException();
            }
            this.taskName = name;
            this.totalWork = thetotalWork;
            doBegin(name, thetotalWork, parentMonitor == null, maxLevels);
            if (recordExecutionTime) {
                // elkjs-exclude-start
                startTime = System.nanoTime();
                // CHECKSTYLEOFF
                // In GWT 'System.nanoTime()' is not available, so we resort to milliseconds here.
                if (false)
                // elkjs-exclude-end
                    startTime = System.currentTimeMillis() * 1000l;
                // CHECKSTYLEON
            }
            return true;
        }
    }

    /**
     * Invoked when a task begins, to be overridden by subclasses. This implementation does nothing.
     *
     * @param name task name
     * @param newTotalWork total amount of work for the new task
     * @param topInstance if true, this progress monitor is the top instance
     * @param maxHierarchyLevels the maximal number of reported hierarchy levels, or -1 for infinite levels
     */
    protected void doBegin(final String name, final float newTotalWork,
            final boolean topInstance, final int maxHierarchyLevels) {
    }
    
    @Override
    public final void worked(final float work) {
        if (work > 0 && !closed) {
            internalWorked(work);
        }
    }

    /**
     * Sets the current work counters of this monitor and all parent monitors.
     *
     * @param work amount of work that has been completed
     */
    private void internalWorked(final float work) {
        if (totalWork > 0 && completedWork < totalWork) {
            completedWork += work;
            doWorked(completedWork, totalWork, parentMonitor == null);
            if (parentMonitor != null && parentMonitor.currentChildWork > 0 && maxLevels != 0) {
                parentMonitor.internalWorked(work / totalWork * parentMonitor.currentChildWork);
            }
        }
    }

    /**
     * Invoked when work is done for this progress monitor, to be overridden by
     * subclasses. This implementation does nothing.
     *
     * @param thecompletedWork total number of work that is done for this task
     * @param thetotalWork total number of work that is targeted for completion
     * @param topInstance if true, this progress monitor is the top instance
     */
    protected void doWorked(final float thecompletedWork, final float thetotalWork,
            final boolean topInstance) {
    }
    
    /** factor for nanoseconds. */
    private static final double NANO_FACT = 1e-9;

    @Override
    public final void done() {
        if (taskName == null) {
            throw new IllegalStateException("The task has not begun yet.");
        }
        if (!closed) {
            if (recordExecutionTime) {
                long endTime;
                // elkjs-exclude-start
                endTime = System.nanoTime();
                // CHECKSTYLEOFF
                // In GWT 'System.nanoTime()' is not available, so we resort to milliseconds here.
                if (false)
                // elkjs-exclude-end
                    endTime = System.currentTimeMillis() * 1000l;
                // CHECKSTYLEON

                totalTime = (endTime - startTime) * NANO_FACT;
            }
            if (completedWork < totalWork) {
                internalWorked(totalWork - completedWork);
            }
            doDone(parentMonitor == null, maxLevels);
            closed = true;
        }
    }

    /**
     * Invoked when a task ends, to be overridden by subclasses. This
     * implementation does nothing.
     *
     * @param topInstance if true, this progress monitor is the top instance
     * @param maxHierarchyLevels the maximal number of reported hierarchy levels, or -1
     *     for infinite levels
     */
    protected void doDone(final boolean topInstance, final int maxHierarchyLevels) {
    }
    
    @Override
    public boolean isRunning() {
        return taskName != null && !closed;
    }

    /**
     * This implementation always returns {@code false}.
     *
     * @return {@code false}
     */
    public boolean isCanceled() {
        return false;
    }

    @Override
    public String getTaskName() {
        return taskName;
    }
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Sub-Tasks

    @Override
    public final IElkProgressMonitor subTask(final float work) {
        if (!closed) {
            BasicProgressMonitor subMonitor = doSubTask(work, maxLevels);
            children.add(subMonitor);
            subMonitor.parentMonitor = this;
            currentChildWork = work;
            return subMonitor;
        } else {
            return null;
        }
    }

    /**
     * Invoked when a sub-task is created, to be overridden by subclasses. Regarding debug logging, the new instance
     * should be configured just like the instance the method is called on. This implementation creates a new
     * {@code BasicProgressMonitor} instance.
     *
     * @param work amount of work that is completed in the current monitor
     *         instance when the sub-task ends
     * @param maxHierarchyLevels the maximal number of reported hierarchy levels for the parent
     *         progress monitor, or -1 for infinite levels
     * @return a new progress monitor instance
     */
    protected BasicProgressMonitor doSubTask(final float work, final int maxHierarchyLevels) {
        int newMaxHierarchyLevels = maxHierarchyLevels > 0 ? maxHierarchyLevels - 1 : maxHierarchyLevels;
        return new BasicProgressMonitor()
                .withMaxHierarchyLevels(newMaxHierarchyLevels)
                .withLogging(recordLogs)
                .withLogPersistence(persistLogs)
                .withExecutionTimeMeasurement(recordExecutionTime);
    }

    @Override
    public final List<IElkProgressMonitor> getSubMonitors() {
        assert children != null;
        return Collections.unmodifiableList(children);
    }

    @Override
    public final IElkProgressMonitor getParentMonitor() {
        return parentMonitor;
    }
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Debugging
    
    @Override
    public boolean isLoggingEnabled() {
        return recordLogs;
    }
    
    @Override
    public boolean isLogPersistenceEnabled() {
        return persistLogs;
    }

    @Override
    public void log(final Object object) {
        if (recordLogs && object != null) {
            String logMessage = object.toString();
            logMessages.add(logMessage);

            // elkjs-exclude-start
            if (persistLogs) {
                // Write to the log file
                Path outputFile = retrieveLogFilePath();
                if (outputFile != null) {
                    try {
                        Files.write(
                                outputFile,
                                new ArrayList<>(java.util.Arrays.asList(logMessage)),
                                StandardOpenOption.APPEND,
                                StandardOpenOption.CREATE);
                    } catch (IOException e) {
                        // We ignore problems writing to the log file
                    }
                }
            }
            // elkjs-exclude-end
        }
    }

    @Override
    public List<String> getLogs() {
        return logMessages == null
                ? null
                : Collections.unmodifiableList(logMessages);
    }
    
    @Override
    public void logGraph(final ElkNode graph, final String tag) {
        if (recordLogs && graph != null) {
            logGraph(EcoreUtil.copy(graph), tag, LoggedGraph.Type.ELK);
        }
    }
    
    @Override
    public void logGraph(final Object object, final String tag, final LoggedGraph.Type graphType) {
        if (recordLogs && object != null && graphType != null) {
            LoggedGraph loggedGraph = new LoggedGraph(object, tag, graphType);
            logGraphs.add(loggedGraph);

            // elkjs-exclude-start
            if (persistLogs) {
                // Find out which file to write to
                String actualTag = (tag == null || tag.isEmpty()) ? "Unnamed" : tag;
                Path filePath = retrieveFilePath(actualTag, graphType.getFileExtension());
                
                // Write to the file
                if (filePath != null) {
                    try {
                        Files.write(
                                filePath,
                                new ArrayList<>(java.util.Arrays.asList(loggedGraph.serialize())),
                                StandardOpenOption.WRITE,
                                StandardOpenOption.CREATE);
                    } catch (IOException e) {
                        // We ignore problems writing to the log file
                    }
                }
            }
            // elkjs-exclude-end
        }
    }

    @Override
    public List<LoggedGraph> getLoggedGraphs() {
        return logGraphs == null
                ? null
                : Collections.unmodifiableList(logGraphs);
    }

    // elkjs-exclude-start
    @Override
    public Path getDebugFolder() {
        initDebugFolder(false);
        return debugFolder;
    }
    // elkjs-exclude-end

    @Override
    public boolean isExecutionTimeMeasured() {
        return recordExecutionTime;
    }

    @Override
    public final double getExecutionTime() {
        return totalTime;
    }
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Debugging Persistence
    
    /**
     * Initializes {@link #debugFolder}.
     * 
     * @param ensureExistence if {@code true}, tries to ensure that the folder exists.
     */
    private void initDebugFolder(final boolean ensureExistence) {
        // elkjs-exclude-start
        if (recordLogs && persistLogs && debugFolder == null) {
            // Our debug path hasn't been computed yet. How we do so depends on whether we have a parent or not
            if (getParentMonitor() == null) {
                initRootMonitorDebugFolder();
            } else {
                initChildMonitorDebugFolder();
            }
            
            if (ensureExistence) {
                // Check various error conditions
                if (Files.isRegularFile(debugFolder)) {
                    // The folder name refers to a file -- panic!
                    recordLogs = false;
                    log("Debug folder '" + debugFolder.toString() + "' refers to a file! Not persisting logs.");
                    debugFolder = null;
                    
                } else if (!Files.exists(debugFolder)) {
                    debugFolder.toFile().mkdirs();
                }
                
                // Only use the debug folder if it exists now
                if (!Files.isDirectory(debugFolder)) {
                    recordLogs = false;
                    log("Unable to create debug folder '" + debugFolder.toString() + "'! Not persisting logs.");
                    debugFolder = null;
                }
            }
        }
        // elkjs-exclude-end
    }
    
    /**
     * Sets {@link #debugFolder} for a root monitor. Does not create the folder.
     */
    private void initRootMonitorDebugFolder() {
        // Retrieve the components of our name
        long timestamp = System.currentTimeMillis();
        
        int validCharacterRange = 'z' - 'a' + 1;
        char randChar1 = (char) ('a' + (int) (Math.random() * validCharacterRange));
        char randChar2 = (char) ('a' + (int) (Math.random() * validCharacterRange));
        
        String name = getTaskName();
        if (name == null || name.isEmpty()) {
            name = "Unnamed";
        }

        // elkjs-exclude-start
        name = ElkUtil.toSafePathName(name);
        
        // Assemble our folder name
        String monitorFolder = String.format("%d%c%c-%s",
                timestamp,
                randChar1,
                randChar2,
                name);
        debugFolder = Paths.get(ElkUtil.debugFolderPath(ROOT_DEBUG_FOLDER_NAME, monitorFolder));
        // elkjs-exclude-end
    }
    
    /**
     * Sets {@link #debugFolder} for a non-root monitor. Does not create the folder.
     */
    private void initChildMonitorDebugFolder() {
        // Retrieve the components of our name
        int index = getParentMonitor().getSubMonitors().indexOf(this);
        
        String name = getTaskName();
        if (name == null || name.isEmpty()) {
            name = "Unnamed";
        }
        // elkjs-exclude-start
        name = ElkUtil.toSafePathName(name);
        
        // Assemble our folder name
        String monitorFolder = String.format("%02d-%s",
                index,
                name);
        debugFolder = getParentMonitor().getDebugFolder().resolve(monitorFolder);
	    // elkjs-exclude-end
    }
    
    /**
     * Returns the path to this monitor's log text file.
     */
    // elkjs-exclude-start
    private Path retrieveLogFilePath() {
        // Determine the log file's path if we haven't already
        if (logFile == null) {
            logFile = retrieveFilePath("log", "txt");
        }
        
        return logFile;
    }
    // elkjs-exclude-end

    /**
     * Returns the path to a file in our monitor's debug output folder. The file name will be based on the given base
     * name and the given extension. If a file with that name already exists, numbers are appended to the name until we
     * find one that doesn't exist.
     * 
     * @param name the base file name (the part before the extension).
     * @param extension the file's extension, without leading period.
     * @return path to a file that doesn't exist.
     */
    // elkjs-exclude-start
    private Path retrieveFilePath(final String name, final String extension) {
        initDebugFolder(true);
        if (debugFolder == null) {
            return null;
        }
        
        // We try different file names until we find one that doesn't already exist
        Path filePath = debugFolder.resolve(name + "." + extension);
        
        int number = 0;
        while (Files.exists(filePath)) {
            // Try the next possible file name
            number++;
            filePath = debugFolder.resolve(name + "-" + number + "." + extension);
        }
        
        // We now have a path to a file that doesn't exist
        return filePath;
    }
    // elkjs-exclude-end

}
