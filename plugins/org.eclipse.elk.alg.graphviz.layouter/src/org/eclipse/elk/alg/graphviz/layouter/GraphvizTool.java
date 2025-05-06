/*******************************************************************************
 * Copyright (c) 2008, 2015 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.alg.graphviz.layouter;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileFilter;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.List;

import org.eclipse.elk.alg.graphviz.dot.transform.Command;
import org.eclipse.elk.alg.graphviz.layouter.preferences.GraphvizLayouterPreferenceStoreAccess;
import org.eclipse.elk.alg.graphviz.layouter.preferences.GraphvizPreferencePage;
import org.eclipse.elk.core.util.Maybe;
import org.eclipse.elk.core.util.WrappedException;
import org.eclipse.jface.preference.PreferenceDialog;
import org.eclipse.swt.widgets.Display;
import org.eclipse.ui.PlatformUI;
import org.eclipse.ui.dialogs.PreferencesUtil;

import com.google.common.collect.Lists;

/**
 * Handler for accessing Graphviz via a separate process.
 * 
 * @author msp
 */
public class GraphvizTool {
    
    /**
     * Available cleanup modes.
     */
    public static enum Cleanup {
        /** normal cleanup. */
        NORMAL,
        /** read error output and stop the Graphviz process. */
        ERROR,
        /** stop the Graphviz process and the watcher thread. */
        STOP;
    }

    /** preference constant for Graphviz executable. */
    public static final String PREF_GRAPHVIZ_EXECUTABLE = "graphviz.executable";
    /** preference constant for timeout. */
    public static final String PREF_TIMEOUT = "graphviz.timeout";
    /** default timeout for waiting for Graphviz to give some output. */
    public static final AtomicInteger PROCESS_DEF_TIMEOUT = new AtomicInteger(20000);
    /** minimal timeout for waiting for Graphviz to give some output. */
    public static final int PROCESS_MIN_TIMEOUT = 200;
    
    /** different names for the Windows program files folder. */
    private static final String[] PROGRAM_FILES_FOLDERS = {
        "Program Files", "Program Files (x86)", "Programme", "Programme (x86)"
    };
    /**
     * Default locations of the dot executable. Each entry ends with the path separator, so that the
     * dot executable's file name can be directly appended.
     */
    private static final List<String> DEFAULT_LOCS = new ArrayList<String>();

    /** argument used to specify the command. */
    private static final String ARG_COMMAND = "-K";
    /** argument to suppress warnings of the executable. */
    private static final String ARG_NOWARNINGS = "-q";
    /** argument to invert the Y axis to conform with SWT. */
    private static final String ARG_INVERTYAXIS = "-y";

    /** the process instance that is used for multiple layout runs. */
    private Process process;
    /** the command that is used to create the process. */
    private Command command;
    /** the watcher thread used to cancel a blocked read operation. */
    private Watchdog watchdog;
    /** the input stream given by the Graphviz process. */
    private InputStream graphvizStream;
    
    
    static {
        // Add all paths from the system PATH variable to the list of paths we will look for dot in
        // to our list of default locations
        String envPath = System.getenv("PATH");
        if (envPath != null) {
            String[] envPaths = envPath.split(File.pathSeparator);
            
            for (int i = 0; i < envPaths.length; i++) {
                if (envPaths[i].trim().length() > 0) {
                    String path; 
                    if (envPaths[i].startsWith("\"") && envPaths[i].endsWith("\"")) {
                        path = envPaths[i].substring(1, envPaths[i].length() - 1);
                    } else {
                        path = envPaths[i];
                    }
                    if (path.endsWith(File.separator)) {
                        DEFAULT_LOCS.add(path);
                    } else {
                        DEFAULT_LOCS.add(path + File.separator);
                    }
                }
            }
        }
        
        if (File.separator.equals("/")) {
            // Fallback list of default locations for Unix-like environments
            DEFAULT_LOCS.add("/opt/local/bin/");
            DEFAULT_LOCS.add("/usr/local/bin/");
            DEFAULT_LOCS.add("/usr/bin/");
        } else if (File.separator.equals("\\")) {
            // If we're on Windows, we try to find the default Graphviz installation directory in the
            // program files folder
            for (String programFilesName : PROGRAM_FILES_FOLDERS) {
                File programFilesFolder = new File("C:\\" + programFilesName);
                if (programFilesFolder.exists()
                        && programFilesFolder.isDirectory()
                        && programFilesFolder.canRead()) {
                    
                    // Find Graphviz directories
                    File[] graphvizDirs = programFilesFolder.listFiles(new FileFilter() {
                        public boolean accept(final File pathname) {
                            return pathname.isDirectory()
                                    && pathname.canRead()
                                    && pathname.getName().toLowerCase().startsWith("graphviz");
                        }
                    });
                    
                    // Add each Graphviz directory
                    if (graphvizDirs != null) {
                        for (File graphvizDir : graphvizDirs) {
                            DEFAULT_LOCS.add(
                                    graphvizDir.toString() + File.separator + "bin" + File.separator);
                        }
                    }
                }
            }
        }
    }
    
    
    /**
     * Create a Graphviz tool instance for the given command.
     * 
     * @param thecommand a Graphviz command
     */
    public GraphvizTool(final Command thecommand) {
        if (thecommand == Command.INVALID) {
            throw new IllegalArgumentException("Invalid Graphviz command.");
        }
        this.command = thecommand;
    }

    /**
     * Initialize the Graphviz tool instance by starting the dot process and the
     * watcher thread as necessary.
     */
    public void initialize() {
        initialize(null);
    }
    
    /**
     * Initialize the Graphviz tool instance by starting the dot process and the watcher
     * thread as necessary. The given command line arguments are appended to the default
     * arguments.
     * 
     * @param arguments command line arguments to be added to the default list of arguments.
     *                  May be {@code null} or empty.
     */
    public synchronized void initialize(final List<String> arguments) {
        if (watchdog == null) {
            // start the watcher thread for timeout checking
            watchdog = new Watchdog();
            watchdog.setName("Graphviz Watchdog");
            watchdog.start();
        }

        if (process == null) {
            String dotExecutable = getDotExecutable();
            
            // assemble the final list of command-line arguments
            List<String> args = Lists.newArrayList(
                    dotExecutable,
                    ARG_NOWARNINGS,
                    ARG_INVERTYAXIS,
                    ARG_COMMAND + command);
            
            if (arguments != null) {
                args.addAll(arguments);
            }
            
            // create the process
            try {
                process = Runtime.getRuntime().exec(args.toArray(new String[args.size()]));
            } catch (IOException exception) {
                throw new WrappedException("Failed to start Graphviz process."
                        + " Please check your Graphviz installation.", exception);
            }
        }
    }

    /**
     * Returns the dot executable path. If it is not found, the user is asked to provide it. Calling
     * this method is equivalent to calling {@code GraphvizTool.getDotExecutable(true)}.
     * 
     * @return path to the dot executable.
     */
    public static String getDotExecutable() {
        String executable = getDotExecutable(true);
        if (executable == null) {
            throw new RuntimeException("The Dot executable was not found in default paths."
                    + " Please check your Graphviz installation.");
        }
        return executable;
    }

    /**
     * Returns the dot executable path.
     * 
     * @param promptUser if the dot executable is not found and this parameter is {@code true}, the user
     *                   is asked to provide the path to the executable. If it is not found and this
     *                   parameter is {@code false}, this method returns {@code null}.
     * @return path to the dot executable, or {@code null} if the executable was not found
     */
    public static String getDotExecutable(final boolean promptUser) {
        String dotExecutable = null;
        File dotFile;
        
        // Load the graphviz path from the preferences, if any. However, do this only if we're really
        // running in an Eclipse context
        if (EclipseRuntimeDetector.isEclipseRunning()) {
            dotExecutable =
                    GraphvizLayouterPreferenceStoreAccess.getUISaveString(PREF_GRAPHVIZ_EXECUTABLE);
            if (dotExecutable != null) {
                dotFile = new File(dotExecutable);
                if (dotFile.exists() && dotFile.canExecute()) {
                    return dotExecutable;
                }
            }
        }
        
        // look in a selection of default locations where it might be installed
        for (String location : DEFAULT_LOCS) {
            // Linux
            dotExecutable = location + "dot";
            dotFile = new File(dotExecutable);
            if (dotFile.exists() && dotFile.canExecute()) {
                return dotExecutable;
            }
            
            // Windows
            dotExecutable = location + "dot.exe";
            dotFile = new File(dotExecutable);
            if (dotFile.exists() && dotFile.canExecute()) {
                return dotExecutable;
            }
        }
        
        // If we haven't found an executable yet, ask the user if so requested and if Eclipse is running
        if (promptUser && EclipseRuntimeDetector.isEclipseRunning()) {
            if (handleExecPath()) {
                // fetch the executable string again after the user has entered a new path
                dotExecutable =
                        GraphvizLayouterPreferenceStoreAccess.getUISaveString(PREF_GRAPHVIZ_EXECUTABLE);
                if (dotExecutable != null) {
                    dotFile = new File(dotExecutable);
                    if (dotFile.exists() && dotFile.canExecute()) {
                        return dotExecutable;
                    }
                }
            }
        }
        
        return null;
    }

    /**
     * Handle missing path to the dot executable. The Graphviz preference page
     * is opened so the user can enter the correct path. The method returns
     * after the preference page has been closed.
     * 
     * @return true if the user has selected "Ok" in the shown dialog, false otherwise
     */
    private static boolean handleExecPath() {
        try {
            final Display display = PlatformUI.getWorkbench().getDisplay();
            final Maybe<Integer> dialogResult = Maybe.create();
            display.syncExec(new Runnable() {
                public void run() {
                    PreferenceDialog preferenceDialog =
                            PreferencesUtil.createPreferenceDialogOn(display.getActiveShell(),
                                    GraphvizPreferencePage.ID, new String[] {}, null);
                    dialogResult.set(preferenceDialog.open());
                }
            });
            return dialogResult.get() == PreferenceDialog.OK;
        } catch (NoClassDefFoundError e) {
            // silent
        }
        return false;
    }

    /**
     * Return the stream that is used to give input to Graphviz.
     * 
     * @return an output stream for writing to the tool
     */
    public OutputStream input() {
        if (process != null) {
            return new BufferedOutputStream(process.getOutputStream());
        }
        throw new IllegalStateException("Graphviz tool has not been initialized.");
    }
    
    /**
     * Return the stream for reading the output of the Graphviz process.
     * 
     * @return an input stream for reading from the tool
     */
    public InputStream output() {
        if (process != null) {
            synchronized (nextJob) {
                // create an input stream and make it visible to the watcher thread
                graphvizStream = new GraphvizStream(process.getInputStream());
                // wake the watcher, which will then sleep until a timeout occurs
                nextJob.notify();
            }
            return graphvizStream;
        }
        throw new IllegalStateException("Graphviz tool has not been initialized.");
    }
    
    /** maximal number of characters to read from error stream. */
    private static final int MAX_ERROR_OUTPUT = 512;
    /** time to wait before checking process errors. */
    private static final int PROC_ERROR_TIME = 500;
    
    /**
     * Clean up, optionally preparing the tool for the next use.
     * 
     * @param c the cleanup option
     */
    public synchronized void cleanup(final Cleanup c) {
        StringBuilder error = null;
        if (process != null) {
            InputStream errorStream = process.getErrorStream();
            try {
                if (c == Cleanup.ERROR && graphvizStream != null) {
                    // wait a bit so the process can either terminate or generate error
                    Thread.sleep(PROC_ERROR_TIME);
                    // read the error stream to display a meaningful error message
                    error = new StringBuilder();
                    int ch;
                    do {
                        ch = errorStream.read();
                        if (ch >= 0) {
                            error.append((char) ch);
                        }
                    } while (error.length() < MAX_ERROR_OUTPUT && ch >= 0);
                    if (error.length() == 0) {
                        // no error message -- check for exit value
                        int exitValue = process.exitValue();
                        if (exitValue != 0) {
                            exitValueError(exitValue, error);
                        }
                    }
                }
                // if error stream is not empty, the process may not terminate
                while (errorStream.available() > 0) {
                    errorStream.read();
                }
            } catch (Exception ex) {
                // ignore exception
            }
            // terminate the Graphviz process if requested
            if (c == Cleanup.ERROR || c == Cleanup.STOP) {
                try {
                    process.getOutputStream().close();
                    process.getInputStream().close();
                } catch (IOException exception) {
                    // ignore exception
                }
                process.destroy();
                process = null;
            }
        }
        
        if (error == null) {
            synchronized (nextJob) {
                // reset the stream to indicate that the job is done
                graphvizStream = null;
                if (watchdog != null) {
                    // wake the watcher if it is still waiting for timeout
                    watchdog.interrupt();
                    // if requested, reset the watcher to indicate that it should terminate
                    if (c == Cleanup.STOP) {
                        watchdog = null;
                    }
                }
            }
        } else if (error.length() > 0) {
            // an error output could be read from Graphviz, so display that to the user
            throw new GraphvizException("Graphviz error: " + error.toString());
        }
    }
    
    /**
     * Generate an error message for the given exit value.
     * 
     * @param exitValue an exit value
     * @param error a string builder for error messages
     */
    private void exitValueError(final int exitValue, final StringBuilder error) {
        error.append("Process terminated with exit value ").append(exitValue);
        // CHECKSTYLEOFF MagicNumber
        if (exitValue > 128) {
            switch (exitValue - 128) {
            case 2: // SIGINT
                error.append(" (interrupted)");
                break;
            case 3: // SIGQUIT
                error.append(" (quit)");
                break;
            case 4: // SIGILL
                error.append(" (illegal instruction)");
                break;
            case 6: // SIGABRT
                error.append(" (aborted)");
                break;
            case 8: // SIGFPE
                error.append(" (floating point error)");
                break;
            case 9: // SIGKILL
                error.append(" (killed)");
                break;
            case 11: // SIGSEGV
                error.append(" (segmentation fault)");
                break;
            case 13: // SIGPIPE
                error.append(" (broken pipe)");
                break;
            case 15: // SIGTERM
                error.append(" (terminated)");
                break;
            }
        }
        // CHECKSTYLEON MagicNumber
        error.append('.');
    }

    /**
     * A specialized input stream for reading data from the Graphviz process.
     */
    private static class GraphvizStream extends InputStream {
        
        /** the stream of process data output. */
        private InputStream stream;
        /** how many opening curly braces have been read that haven't closed yet. */
        private int depth = 0;
        /** how many graphs have been completely processed. */
        private int finished = 0;
        /** buffered character to return on the next read. */
        private int buf = -1;
        
        /**
         * Create a Graphviz input stream.
         * 
         * @param thestream the process stream to read from
         */
        GraphvizStream(final InputStream thestream) {
            this.stream = thestream;
        }

        @Override
        public int read() throws IOException {
            if (buf >= 0) {
                int c = buf;
                buf = -1;
                return c;
            }
            // don't block if we already finished reading a graph
            if (finished > 0 && stream.available() == 0) {
                return -1;
            }
            
            // track the opening and closing braces while reading the stream
            int c = stream.read();
            if (c == '{') {
                depth++;
            } else if (c == '}') {
                depth--;
                if (depth == 0) {
                    finished++;
                }
            } else if (c == '\\') {
                // discard any line breaks that have been escaped
                buf = stream.read();
                if (buf == '\n' || buf == '\r') {
                    c = stream.read();
                    if (buf == '\r' && c == '\n') {
                        c = stream.read();
                    }
                    buf = -1;
                }
            }
            return c;
        }
        
        @Override
        public int available() throws IOException {
            return stream.available();
        }
            
    }
    
    /** synchronization object between the main thread and the watcher thread. */
    private Object nextJob = new Object();
    
    /**
     * A watcher thread that takes action when a timeout occurs.
     */
    private class Watchdog extends Thread {
        
        @Override
        public void run() {
            do {
                synchronized (nextJob) {
                    // the watcher starts working as soon as a stream is made visible
                    while (graphvizStream == null) {
                        try {
                            // wait for notification by the main thread
                            nextJob.wait();
                        } catch (InterruptedException ex) {
                            // an interrupt can happen when a shutdown is requested
                            if (watchdog == null) {
                                return;
                            }
                        }
                    }
                }
                
                // retrieve the current timeout value
                int timeout = PROCESS_DEF_TIMEOUT.get();
                if (EclipseRuntimeDetector.isEclipseRunning()) {
                    int timeoutPreference =
                            GraphvizLayouterPreferenceStoreAccess.getUISaveInt(PREF_TIMEOUT);
                    if (timeoutPreference >= PROCESS_MIN_TIMEOUT) {
                        timeout = timeoutPreference;
                    }
                }

                boolean interrupted = false;
                try {
                    Thread.sleep(timeout);
                }  catch (InterruptedException ex) {
                    // this means the main thread has done a cleanup before the timeout occurred
                    interrupted = true;
                }
                
                if (!interrupted) {
                    synchronized (nextJob) {
                        // timeout has occurred! kill the process so the main thread will wake
                        Process myProcess = process;
                        if (myProcess != null) {
                            graphvizStream = null;
                            myProcess.destroy();
                        }
                    }
                }
                
            } while (watchdog != null);
        }
    }

}
