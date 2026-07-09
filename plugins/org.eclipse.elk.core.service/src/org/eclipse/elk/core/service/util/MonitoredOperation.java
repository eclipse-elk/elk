/*******************************************************************************
 * Copyright (c) 2010, 2015 Kiel University and others.
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
package org.eclipse.elk.core.service.util;

import java.util.LinkedList;
import java.util.concurrent.ExecutorService;

import org.eclipse.core.runtime.IProgressMonitor;
import org.eclipse.core.runtime.IStatus;
import org.eclipse.core.runtime.NullProgressMonitor;
import org.eclipse.core.runtime.Status;
import org.eclipse.elk.core.service.DiagramLayoutEngine;
import org.eclipse.elk.core.service.ElkServicePlugin;
import org.eclipse.elk.core.util.BasicProgressMonitor;
import org.eclipse.elk.core.util.IElkCancelIndicator;
import org.eclipse.elk.core.util.IElkProgressMonitor;
import org.eclipse.elk.core.util.Maybe;
import org.eclipse.jface.operation.IRunnableWithProgress;
import org.eclipse.jface.preference.IPreferenceStore;
import org.eclipse.swt.widgets.Display;
import org.eclipse.ui.PlatformUI;
import org.eclipse.ui.statushandlers.StatusManager;

/**
 * An operation that can be tracked with a progress bar. This operation can be called from any
 * thread, UI or non-UI, and adds only very little overhead for progress tracking compared to the
 * usual synchronous tracking of the workbench progress service. Just implement
 * {@link #execute(IProgressMonitor)} and call {@link #runMonitored()}.
 * 
 * <p>
 * In case no workbench is set up and, thus, the operation shall run off-screen call
 * {@link #runUnmonitored()}.
 * </p>
 * 
 * @author msp
 */
public abstract class MonitoredOperation {

    /** maximal number of recursion levels for which progress is displayed. */
    private static final int MAX_PROGRESS_LEVELS = 4;
    
    /** The executor service used to perform operations. */
    private final ExecutorService executorService;
    /** Timestamp of the operation's starting time. */
    private long timestamp;
    /** Whether the operation has been canceled. */
    private boolean isCanceled;
    /** Additional cancel indicator to query. */
    private final IElkCancelIndicator cancelIndicator;
    
    /**
     * Create a monitored operation with the given executor service.
     * 
     * @param service an executor service for performing operations
     */
    public MonitoredOperation(final ExecutorService service) {
        this(service, null);
    }
    
    /**
     * Create a monitored operation with the given executor service and cancel indicator.
     * 
     * @param service an executor service for performing operations
     * @param cancelIndicator an indicator queried regularly for cancelation
     */
    public MonitoredOperation(final ExecutorService service, final IElkCancelIndicator cancelIndicator) {
        if (service == null) {
            throw new NullPointerException();
        }
        this.executorService = service;
        this.cancelIndicator = cancelIndicator;
    }
    
    /**
     * Determine whether the operation has been canceled.
     */
    protected boolean isCanceled() {
        return isCanceled || (cancelIndicator != null && cancelIndicator.isCanceled());
    }
    
    /**
     * Execute the monitored operation.
     * 
     * @param monitor the progress monitor for the operation
     * @return a status indicating success or failure
     */
    protected abstract IStatus execute(IElkProgressMonitor monitor);
    
    /**
     * Executed in the UI thread before the operation starts. The default
     * implementation does nothing.
     */
    protected void preUIexec() {
    }
    
    /**
     * Executed in the UI thread after the operation has ended. The default
     * implementation does nothing.
     */
    protected void postUIexec() {
    }

    /**
     * Factory method providing {@link IElkProgressMonitor IElkProgressMonitors}.
     * May be overridden in order to contributed specialized implementations.
     *
     * @return a new {@link IElkProgressMonitor} instance each time it is called
     */
    protected IElkProgressMonitor createMonitor() {
        return new CancelableProgressMonitor();
    }
    
    /**
     * Run the operation. If the current thread is the UI thread, the actual operation
     * is executed in a new thread that runs in parallel. Otherwise the operation is
     * executed directly and a handler for the progress monitor is executed in the
     * parallel UI thread. In either case the method returns only after execution of
     * the operation is done.
     */
    public final void runMonitored() {
        timestamp = System.currentTimeMillis();
        Display display = Display.getCurrent();
        if (display == null) {
            display = PlatformUI.getWorkbench().getDisplay();
            runMonitored(display, false);
        } else {
            runMonitored(display, true);
        }
    }
    
    /**
     * Run the operation without any progress monitor. If the current thread is the UI thread,
     * the actual operation is executed in a new thread that runs in parallel. Otherwise the
     * operation is executed directly and a handler for the progress monitor is executed in the
     * parallel UI thread. In either case the method returns only after execution of
     * the operation is done.
     * 
     */
    public final void runUnmonitored() {
        timestamp = System.currentTimeMillis();

        if (PlatformUI.isWorkbenchRunning()) {
            final Display display = Display.getCurrent();

            if (display == null) {
                runUnmonitored(PlatformUI.getWorkbench().getDisplay(), false);
            } else {
                runUnmonitored(display, true);
            }
        } else {
            runOffscreen();
        }
    }
    
    /**
     * Return the timestamp of the operation's start, or 0 if it has not started yet.
     * 
     * @return the timestamp
     */
    public long getTimestamp() {
        return timestamp;
    }
    
    /**
     * Cancel the operation. This flag can be queried with the
     * {@link IElkProgressMonitor#isCanceled()} method of the progress monitor passed in
     * {@link #execute(IElkProgressMonitor)}.
     */
    public void cancel() {
        isCanceled = true;
    }


    /**
     * Run the operation without switching to any other thread.
     * 
     * @author chsch
     */
    private void runOffscreen() {
        final Maybe<IStatus> status = new Maybe<IStatus>();
        
        try {
            // execute the preparation prior to the actual operation 
            preUIexec();
        } catch (Throwable throwable) {
            status.set(new Status(IStatus.ERROR, ElkServicePlugin.PLUGIN_ID,
                    "Error in monitored operation running offscreen during prepration", throwable));
        }

        if (status.get() == null && !isCanceled()) {
            // execute the actual operation without progress monitor
            status.set(execute(createMonitor()));
        }

        if (status.get() != null && status.get().getSeverity() == IStatus.OK && !isCanceled()) {
            // execute the post processing code after the actual operation
            try {
                postUIexec();
            } catch (Throwable throwable) {
                status.set(new Status(IStatus.ERROR, ElkServicePlugin.PLUGIN_ID,
                        "Error in monitored operation running offscreen during post processing",
                        throwable));
            }
        }

        handleStatus(status);
    }
    
    /**
     * Run the operation from the current thread without a progress monitor. 
     * 
     * @param display the display that runs the UI thread
     * @param isUiThread if true, the current thread is the UI thread
     */
    private void runUnmonitored(final Display display, final boolean isUiThread) {
        final Maybe<IStatus> status = new Maybe<IStatus>();
        
        if (isUiThread) {
            try {
                // execute UI code prior to the actual operation
                preUIexec();

                // execute the actual operation without progress monitor
                synchronized (executorService) {
                    executorService.execute(new Runnable() {
                        public void run() {
                            status.set(execute(createMonitor()));
                            assert status.get() != null;
                            display.wake();
                        }
                    });
                }
                
                while (status.get() == null && !isCanceled()) {
                    boolean hasMoreToDispatch;
                    do {
                        hasMoreToDispatch = display.readAndDispatch();
                    } while (hasMoreToDispatch && status.get() == null && !isCanceled());
                    if (status.get() == null && !isCanceled()) {
                        display.sleep();
                    }
                }
                
                if (status.get() != null && status.get().getSeverity() == IStatus.OK && !isCanceled()) {
                    // execute UI code after the actual operation
                    postUIexec();
                }
            } catch (Throwable throwable) {
                status.set(new Status(IStatus.ERROR, ElkServicePlugin.PLUGIN_ID,
                        "Error in monitored operation", throwable));
            }
        } else {
            // execute UI code prior to the actual operation
            display.syncExec(new Runnable() {
                public void run() {
                    try {
                        preUIexec();
                    } catch (Throwable throwable) {
                        status.set(new Status(IStatus.ERROR, ElkServicePlugin.PLUGIN_ID,
                                "Error in monitored operation", throwable));
                    }
                }
            });
            
            if (status.get() == null && !isCanceled()) {
                // execute the actual operation without progress monitor
                status.set(execute(createMonitor()));
                
                if (status.get().getSeverity() == IStatus.OK && !isCanceled()) {
                    // execute UI code after the actual operation
                    display.syncExec(new Runnable() {
                        public void run() {
                            try {
                                postUIexec();
                            } catch (Throwable throwable) {
                                status.set(new Status(IStatus.ERROR,
                                        ElkServicePlugin.PLUGIN_ID,
                                        "Error in monitored operation", throwable));
                            }
                        }
                    });
                }
            }
        }
        
        handleStatus(status);
    }
    
    /**
     * Run the operation from the current thread. 
     * 
     * @param display the display that runs the UI thread
     * @param isUiThread if true, the current thread is the UI thread
     */
    private void runMonitored(final Display display, final boolean isUiThread) {
        final Maybe<IProgressMonitor> monitor = new Maybe<IProgressMonitor>();
        final Maybe<IStatus> status = new Maybe<IStatus>();
        
        if (isUiThread) {
            // execute the operation in a new thread and the UI code in the current thread
            synchronized (executorService) {
                executorService.execute(new Runnable() {
                    public void run() {
                        runOperation(display, monitor, status);
                    }
                });
            }
            runUiHandler(display, monitor, status);
        } else {
            // execute the operation in the current thread and the UI code in the UI thread
            display.asyncExec(new Runnable() {
                public void run() {
                    runUiHandler(display, monitor, status);
                }
            });
            runOperation(display, monitor, status);
        }
        
        handleStatus(status);
    }
    
    /**
     * Handle a status object.
     * 
     * @param status a status
     */
    protected void handleStatus(final Maybe<IStatus> status) {
        if (status.get() != null) {
            int handlingStyle = StatusManager.NONE;
            switch (status.get().getSeverity()) {
            case IStatus.ERROR:
                handlingStyle = StatusManager.SHOW | StatusManager.LOG;
                break;
            case IStatus.WARNING:
            case IStatus.INFO:
                handlingStyle = StatusManager.LOG;
                break;
            }
            StatusManager.getManager().handle(status.get(), handlingStyle);
        }
    }
    
    /**
     * Runs the operation after synchronization with the UI handler.
     * 
     * @param display the current display
     * @param monitor the progress monitor wrapper
     * @param status the returned status
     */
    private void runOperation(final Display display, final Maybe<IProgressMonitor> monitor,
            final Maybe<IStatus> status) {
        try {
            synchronized (monitor) {
                while (monitor.get() == null && !isCanceled()) {
                    try {
                        monitor.wait();
                    } catch (InterruptedException exception) {
                        // ignore exception
                    }
                }
            }
            if (status.get() == null && !isCanceled()) {
                IPreferenceStore prefStore = ElkServicePlugin.getInstance().getPreferenceStore();
                IElkProgressMonitor monAdapt = new ProgressMonitorAdapter(monitor.get(), MAX_PROGRESS_LEVELS)
                        .withLogging(prefStore.getBoolean(DiagramLayoutEngine.PREF_DEBUG_LOGGING))
                        .withExecutionTimeMeasurement(prefStore.getBoolean(DiagramLayoutEngine.PREF_DEBUG_EXEC_TIME));
                status.set(execute(monAdapt));
                assert status.get() != null;
            }
        } finally {
            synchronized (status) {
                if (status.get() == null) {
                    status.set(Status.OK_STATUS);
                }
                display.wake();
            }
        }
    }
    
    /**
     * Runs the UI handler.
     * 
     * @param display the current display
     * @param monitor the progress monitor wrapper
     * @param status the returned status
     */
    private void runUiHandler(final Display display, final Maybe<IProgressMonitor> monitor,
            final Maybe<IStatus> status) {
        try {
            
            // execute UI code prior to the actual operation
            preUIexec();
            
            // execute UI handler until the operation has finished
            PlatformUI.getWorkbench().getProgressService().run(
                    false, true, new IRunnableWithProgress() {
                public void run(final IProgressMonitor uiMonitor) {
                    ProgressMonitorWrapper monitorWrapper = new ProgressMonitorWrapper(display);
                    synchronized (monitor) {
                        monitor.set(monitorWrapper);
                        monitor.notify();
                    }
                    while (status.get() == null && !isCanceled()) {
                        boolean hasMoreToDispatch = false;
                        do {
                            hasMoreToDispatch = display.readAndDispatch();
                            isCanceled = uiMonitor.isCanceled();
                        } while (hasMoreToDispatch && status.get() == null && !isCanceled());
                        if (status.get() == null && monitorWrapper.commands.isEmpty() && !isCanceled()) {
                            display.sleep();
                        }
                        while (!monitorWrapper.commands.isEmpty() && status.get() == null
                                && !isCanceled()) {
                            WrapperCommand command;
                            synchronized (monitorWrapper.commands) {
                                command = monitorWrapper.commands.removeFirst();
                            }
                            switch (command.type) {
                            case BEGIN_TASK:
                                uiMonitor.beginTask(command.name, command.work);
                                break;
                            case SET_TASK_NAME:
                                uiMonitor.setTaskName(command.name);
                                break;
                            case SUB_TASK:
                                uiMonitor.subTask(command.name);
                                break;
                            case WORKED:
                                uiMonitor.worked(command.work);
                                break;
                            case INTERNAL_WORKED:
                                uiMonitor.internalWorked(command.work);
                            case DONE:
                                uiMonitor.done();
                                return;
                            }
                        }
                    }
                }
            });
            while (status.get() == null && !isCanceled()) {
                boolean hasMoreToDispatch;
                do {
                    hasMoreToDispatch = display.readAndDispatch();
                } while (hasMoreToDispatch && status.get() == null && !isCanceled());
                if (status.get() == null && !isCanceled()) {
                    display.sleep();
                }
            }
            
            // execute UI code after the actual operation
            if (status.get() != null && status.get().getSeverity() == IStatus.OK && !isCanceled()) {
                postUIexec();
            }
            
        } catch (Throwable throwable) {
            synchronized (monitor) {
                if (monitor.get() == null) {
                    monitor.set(new NullProgressMonitor());
                    monitor.notify();
                }
            }
            synchronized (status) {
                if (status.get() == null || status.get().getSeverity() == IStatus.OK) {
                    status.set(new Status(IStatus.ERROR, ElkServicePlugin.PLUGIN_ID,
                            "Error in monitored operation", throwable));
                    handleStatus(status);
                }
            }
        }
    }
    
    
    /**
     * Data type for progress monitor wrapper commands.
     */
    private static final class WrapperCommand {
        
        /**
         * Enumeration of progress monitor wrapper command types.
         */
        enum Type {
            BEGIN_TASK, SET_TASK_NAME, SUB_TASK, WORKED, INTERNAL_WORKED, DONE
        }
        
        /** the command type. */
        private Type type;
        /** the task name for new subtasks. */
        private String name;
        /** the amount of work for progress reports. */
        private int work;
        
        /**
         * Create a progress monitor wrapper command.
         * 
         * @param thetype the command type
         * @param thename the task name for new subtasks
         * @param thework the amount of work for progress reports
         */
        private WrapperCommand(final Type thetype, final String thename, final int thework) {
            this.type = thetype;
            this.name = thename;
            this.work = thework;
        }
        
    }

    /**
     * Wrapper class for Eclipse progress monitors.
     */
    private class ProgressMonitorWrapper implements IProgressMonitor {
        
        /** queue of commands that have to be delegated to the wrapped Eclipse progress monitor. */
        private final LinkedList<WrapperCommand> commands = new LinkedList<WrapperCommand>();
        /** the display that is woken after each incoming command. */
        private final Display display;
        
        /**
         * Create a progress monitor wrapper.
         * 
         * @param thedisplay the display that is woken after each incoming monitor command
         */
        ProgressMonitorWrapper(final Display thedisplay) {
            this.display = thedisplay;
        }
        
        @Override
        public void beginTask(final String name, final int totalWork) {
            synchronized (commands) {
                commands.addLast(new WrapperCommand(WrapperCommand.Type.BEGIN_TASK, name, totalWork));
            }
            display.wake();
        }

        @Override
        public void done() {
            synchronized (commands) {
                commands.addLast(new WrapperCommand(WrapperCommand.Type.DONE, null, 0));
            }
            display.wake();
        }

        @Override
        public void internalWorked(final double work) {
            synchronized (commands) {
                commands.addLast(new WrapperCommand(WrapperCommand.Type.INTERNAL_WORKED, null,
                        (int) work));
            }
            display.wake();
        }

        @Override
        public void setTaskName(final String name) {
            synchronized (commands) {
                commands.addLast(new WrapperCommand(WrapperCommand.Type.SET_TASK_NAME, name, 0));
            }
            display.wake();
        }

        @Override
        public void subTask(final String name) {
            synchronized (commands) {
                commands.addLast(new WrapperCommand(WrapperCommand.Type.SUB_TASK, name, 0));
            }
            display.wake();
        }

        @Override
        public void worked(final int work) {
            synchronized (commands) {
                commands.addLast(new WrapperCommand(WrapperCommand.Type.WORKED, null, work));
            }
            display.wake();
        }

        @Override
        public boolean isCanceled() {
            return isCanceled;
        }

        @Override
        public void setCanceled(final boolean value) {
            isCanceled = value;
        }
        
    }
    
    
    /**
     * A progress monitor that can be canceled with the operation's {@link #cancel()} method.
     */
    protected class CancelableProgressMonitor extends BasicProgressMonitor {
        
        /**
         * Create a cancelable progress monitor.
         */
        public CancelableProgressMonitor() {
            super(0);

            IPreferenceStore prefStore = ElkServicePlugin.getInstance().getPreferenceStore();
            withLogging(prefStore.getBoolean(DiagramLayoutEngine.PREF_DEBUG_LOGGING));
            withLogPersistence(prefStore.getBoolean(DiagramLayoutEngine.PREF_DEBUG_STORE));
            withExecutionTimeMeasurement(prefStore.getBoolean(DiagramLayoutEngine.PREF_DEBUG_EXEC_TIME));
        }
        
        @Override
        public boolean isCanceled() {
            return MonitoredOperation.this.isCanceled();
        }
        
    }
    
}
