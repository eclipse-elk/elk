/*******************************************************************************
 * Copyright (c) 2009, 2015 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.core.service.util;

import org.eclipse.core.runtime.IProgressMonitor;
import org.eclipse.elk.core.util.BasicProgressMonitor;

/**
 * Wrapper class for Eclipse progress monitors.
 * 
 * @author msp
 */
public class ProgressMonitorAdapter extends BasicProgressMonitor {

    /** the Eclipse progress monitor used by this monitor. */
    private IProgressMonitor progressMonitor;
    /** the number of work units that were already submitted. */
    private int submittedWork = 0;

    /**
     * Creates a progress monitor wrapper for a given Eclipse progress monitor. The number
     * of hierarchy levels for which progress is reported to the Eclipse monitor is infinite.
     * 
     * @param theprogressMonitor the progress monitor
     */
    public ProgressMonitorAdapter(final IProgressMonitor theprogressMonitor) {
        super();
        this.progressMonitor = theprogressMonitor;
    }
    
    /**
     * Creates a progress monitor wrapper with given maximal number of hierarchy levels. Progress
     * is reported to parent monitors only up to the specified hierarchy level.
     * 
     * @deprecated Use {@link #withMaxHierarchyLevels(int)} instead.
     * @param theprogressMonitor the progress monitor
     * @param maxLevels maximal number of hierarchy levels for which progress is reported
     */
    public ProgressMonitorAdapter(final IProgressMonitor theprogressMonitor, final int maxLevels) {
        super(maxLevels);
        this.progressMonitor = theprogressMonitor;
    }
    
    /**
     * Reports to the integrated Eclipse progress monitor that the current task begins.
     * 
     * @param name task name
     * @param totalWork total amount of work for the new task
     * @param topInstance if true, this progress monitor is the top instance
     * @param maxHierarchyLevels maximal number of reported hierarchy levels
     */
    @Override
    protected void doBegin(final String name, final float totalWork, final boolean topInstance,
            final int maxHierarchyLevels) {
        
        if (topInstance) {
            progressMonitor.beginTask(name, (int) (totalWork <= 0
                    ? IProgressMonitor.UNKNOWN : totalWork));
        } else if (maxHierarchyLevels != 0) {
            progressMonitor.subTask(name);
        }
    }

    /**
     * Reports to the integrated Eclipse progress monitor that the main task is
     * done, if this is the top instance.
     * 
     * @param topInstance if true, this progress monitor is the top instance
     * @param maxHierarchyLevels maximal number of reported hierarchy levels
     */
    @Override
    protected void doDone(final boolean topInstance, final int maxHierarchyLevels) {
        if (topInstance) {
            progressMonitor.done();
        }
    }

    /**
     * Returns true if the integrated Eclipse progress monitor reports
     * cancellation.
     * 
     * @return true if the user has requested to cancel the operation
     */
    public boolean isCanceled() {
        return progressMonitor.isCanceled();
    }

    /**
     * Creates a new instance of {@code ProgressMonitorAdapter}.
     * 
     * @param work amount of work that is completed in the current monitor
     *            instance when the sub-task ends
     * @param maxHierarchyLevels the maximal number of hierarchy levels for the parent
     *         progress monitor
     * @return a new progress monitor instance
     */
    @Override
    protected BasicProgressMonitor doSubTask(final float work, final int maxHierarchyLevels) {
        int newMaxHierarchyLevels = maxHierarchyLevels > 0 ? maxHierarchyLevels - 1 : maxHierarchyLevels;
        return new ProgressMonitorAdapter(progressMonitor)
                .withMaxHierarchyLevels(newMaxHierarchyLevels)
                .withLogging(isLoggingEnabled())
                .withLogPersistence(isLogPersistenceEnabled())
                .withExecutionTimeMeasurement(isExecutionTimeMeasured());
    }

    /**
     * Reports to the integrated Eclipse progress monitor that some work was
     * done, if this is the top instance.
     * 
     * @param completedWork total number of work that is done for this task
     * @param totalWork total number of work that is targeted for completion
     * @param topInstance if true, this progress monitor is the top instance
     */
    @Override
    protected void doWorked(final float completedWork, final float totalWork, final boolean topInstance) {
        if (topInstance) {
            int newWork = (int) completedWork;
            if (newWork > submittedWork) {
                progressMonitor.worked(newWork - submittedWork);
                submittedWork = newWork;
            }
        }
    }

}
