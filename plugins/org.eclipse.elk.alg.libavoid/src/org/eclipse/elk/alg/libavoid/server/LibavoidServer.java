/*******************************************************************************
 * Copyright (c) 2013, 2022 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.alg.libavoid.server;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.net.URL;
import java.util.HashMap;
import java.util.Map;

import org.eclipse.core.runtime.FileLocator;
import org.eclipse.core.runtime.Path;
import org.eclipse.elk.alg.libavoid.LibavoidPlugin;
import org.osgi.framework.Bundle;

/**
 * Wraps the execution of the libavoid-server binary. Also employs an watchdog in case of errors.
 * 
 * @author uru
 * @author msp
 */
public class LibavoidServer {

    /**
     * Available cleanup modes.
     */
    public static enum Cleanup {
        /** normal cleanup. */
        NORMAL,
        /** read error output and stop the Libavoid process and the watcher thread. */
        ERROR,
        /** stop the Libavoid process and the watcher thread. */
        STOP;
    }

    /**
     * A helper enumeration for identifying the operating system.
     */
    private enum OS {
        LINUX32, LINUX64, WIN32, WIN64, OSX32, OSX64, SOLARIS, UNKNOWN
    }

    /**
     * Detect the operating system from system properties.
     * 
     * @return the operating system
     */
    private static OS detectOS() {
        String os = System.getProperty("os.name").toLowerCase();
        String arch = System.getProperty("os.arch").toLowerCase();
        if (os.contains("linux")) {
            if (arch.contains("64")) {
                return OS.LINUX64;
            } else if (arch.contains("86")) {
                return OS.LINUX32;
            }
        } else if (os.contains("win")) {
            if (arch.contains("64")) {
                return OS.WIN64;
            } else if (arch.contains("86")) {
                return OS.WIN32;
            }
        } else if (os.contains("mac")) {
            if (arch.contains("64")) {
                return OS.OSX64;
            } else if (arch.contains("86")) {
                return OS.OSX32;
            }
        } else if (os.contains("solaris")) {
            return OS.SOLARIS;
        }
        return OS.UNKNOWN;
    }

    /**
     * Constructor only has package visibility. Use {@link LibavoidServerPool} to create instances.
     */
    LibavoidServer() {
    }

    /** The ogdf server executable. */
    private String executable;
    /** The ogdf server process. */
    private Process process;
    /** The watcher thread used to cancel a blocked read operation. */
    private Watchdog watchdog;
    /** The input stream given by the Libavoid process. */
    private InputStream libavoidStream;
    /** A temporary file that should be removed after closing the process. */
    private File tempFile;
    /** Timeout waiting for the Libavoid process */
    private int processTimeout = PROCESS_DEF_TIMEOUT;

    /** the relative path for the linux64 executable. */
    public static final String EXECUTABLE_PATH_LINUX64 = "/libavoid-server/libavoid-server-linux";
    /** the relative path for the win64 executable. */
    public static final String EXECUTABLE_PATH_WIN64 = "/libavoid-server/libavoid-server-win.exe";
    /** the relative path for the osx64 executable. */
    public static final String EXECUTABLE_PATH_OSX64 = "/libavoid-server/libavoid-server-macos";

    /** the size for file transfer buffers. */
    public static final int BUFFER_SIZE = 512;

    /**
     * Resolve the Libavoid server executable.
     * 
     * @param an
     *            executable file
     * @throws IOException
     *             when the executable could not be located
     */
    @SuppressWarnings("incomplete-switch")
    private File resolveExecutable() throws IOException {
        String path = null;
        OS os = detectOS();
        switch (os) {
        case LINUX64:
            path = EXECUTABLE_PATH_LINUX64;
            break;
        case WIN64:
            path = EXECUTABLE_PATH_WIN64;
            break;
        case OSX64:
            path = EXECUTABLE_PATH_OSX64;
            break;
        default:
            throw new LibavoidServerException("Unsupported operating system.");
        }
        URL url = null;
        if (LibavoidPlugin.getDefault() != null) {
	        Bundle bundle = LibavoidPlugin.getDefault().getBundle();
	        url = FileLocator.find(bundle, new Path(path), null);
        }
        if (url == null) {
        	url = getClass().getResource(path);
        }
        if (url == null) {
            throw new LibavoidServerException("Libavoid binary could not be located.");
        }
        File execFile = new File(FileLocator.resolve(url).getFile());

        // if the plug-in is in a jar archive, create a temporary file to execute
        if (!execFile.exists()) {
            execFile = File.createTempFile("libavoid-server", ".exe");
            OutputStream dest = new FileOutputStream(execFile);
            InputStream source = url.openStream();
            byte[] buffer = new byte[BUFFER_SIZE];
            int count;
            do {
                count = source.read(buffer);
                if (count > 0) {
                    dest.write(buffer, 0, count);
                }
            } while (count > 0);
            dest.close();
            tempFile = execFile;
        }

        // set the file permissions if necessary

        switch (os) {
        case LINUX32:
        case LINUX64:
        case OSX32:
        case OSX64:
        case SOLARIS:
            if (!execFile.canExecute()) {
                boolean success = execFile.setExecutable(true);
                if (!success) {
                    throw new LibavoidServerException("Failed to set executable permission for "
                            + execFile.getPath());
                }
            }
            break;
        }
        return execFile;
    }

    /**
     * Initialize the Libavoid server instance by starting the Libavoid process and the watcher
     * thread as necessary.
     */
    public synchronized void initialize() {
        if (watchdog == null) {
            // start the watcher thread for timeout checking
            watchdog = new Watchdog();
            watchdog.setName("Libavoid Watchdog");
            watchdog.start();
        }

        if (process == null) {
            try {
                if (executable == null) {
                    executable = resolveExecutable().getPath();
                }
                process = Runtime.getRuntime().exec(new String[] { executable });
            } catch (IOException exception) {
                throw new LibavoidServerException("Failed to start libavoid server process.", exception);
            } finally {
                if (process == null) {
                    cleanup(Cleanup.STOP);
                }
            }
        }
    }
    
    /**
     * Sets the timeout to an externally defined value.
     * @param timeout Time to wait for the server to return any output.
     */
    public void setProcessTimeout(int timeout) {
        processTimeout = timeout;
    }

    /**
     * Return the stream that is used to give input to Libavoid.
     * 
     * @return an output stream for writing to the tool
     */
    public OutputStream input() {
        if (process != null) {
            return new BufferedOutputStream(process.getOutputStream());
        }
        throw new IllegalStateException("Libavoid server has not been initialized.");
    }

    /**
     * Return the stream for reading the output of the Libavoid process.
     * 
     * @return an input stream for reading from the tool
     */
    private InputStream output() {
        if (process != null) {
            synchronized (nextJob) {
                // create an input stream and make it visible to the watcher thread
                libavoidStream = process.getInputStream();
                // wake the watcher, which will then sleep until a timeout occurs
                nextJob.notify();
            }
            return libavoidStream;
        }
        throw new IllegalStateException("Libavoid server has not been initialized.");
    }

    /**
     * An enumeration for keeping track of the current parser state.
     */
    private enum ParseState {
        TYPE, DATA, ERROR
    }

    /**
     * Read output data from the Libavoid server process.
     * 
     * @return key-value map of output data, or {@code null} if the process output was not complete
     */
    public Map<String, String> readOutputData() {
        Map<String, String> data = new HashMap<String, String>();
        BufferedReader reader = new BufferedReader(new InputStreamReader(output()));
        ParseState state = ParseState.TYPE;
        boolean parseMore = true;
        StringBuilder error = null;

        while (parseMore) {
            String line = null;
            try {
                line = reader.readLine();
                // System.out.println(line);
            } catch (IOException exception) {
                // most probably the stream was closed due to a timeout of the watchdog thread
            }
            if (line == null) {
                // the stream is empty although more input is expected
                return null;
            }
            
            // capture debug output
            if (line.startsWith("DEBUG")){
                System.out.println(line);
                continue;
            } 

            switch (state) {
            case TYPE:
                if (line.equals("LAYOUT")) {
                    state = ParseState.DATA;
                } else if (line.equals("ERROR")) {
                    state = ParseState.ERROR;
                    error = new StringBuilder();
                }
                break;

            case DATA:
                if (line.equals("DONE")) {
                    parseMore = false;
                } else {
                    String[] tokens = line.split("=");
                    if (tokens.length == 2 && tokens[0].length() > 0) {
                        data.put(tokens[0], tokens[1]);
                    }
                }
                break;

            case ERROR:
                if (line.equals("DONE")) {
                    cleanup(Cleanup.STOP);
                    throw new LibavoidServerException(error.toString());
                } else {
                    if (error.length() > 0) {
                        error.append('\n');
                    }
                    error.append(line);
                }
                break;
            }

        }
        return data;
    }

    /** maximal number of characters to read from error stream. */
    private static final int MAX_ERROR_OUTPUT = 512;
    /** time to wait before checking process errors. */
    private static final int PROC_ERROR_TIME = 500;

    /**
     * Clean up, optionally preparing the tool for the next use.
     * 
     * @param c
     *            the cleanup option
     */
    public synchronized void cleanup(final Cleanup c) {
        StringBuilder error = null;
        if (process != null) {
            InputStream errorStream = process.getErrorStream();
            try {
                if (c == Cleanup.ERROR) {
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

            // terminate the Libavoid process if requested
            if (c == Cleanup.ERROR || c == Cleanup.STOP) {
                try {
                    process.getOutputStream().close();
                    process.getInputStream().close();
                } catch (IOException exception) {
                    // ignore exception
                }
                process.destroy();
                process = null;

                if (tempFile != null) {
                    tempFile.delete();
                    tempFile = null;
                }
            }
        }

        synchronized (nextJob) {
            // reset the stream to indicate that the job is done
            libavoidStream = null;
            if (watchdog != null) {
                Watchdog myWatchdog = watchdog;
                // if requested, reset the watcher to indicate that it should terminate
                if (c == Cleanup.ERROR || c == Cleanup.STOP) {
                    watchdog = null;
                }
                // wake the watcher if it is still waiting for timeout
                myWatchdog.interrupt();
            }
        }

        if (error != null && error.length() > 0) {
            // an error output could be read from Libavoid, so display that to the user
            throw new LibavoidServerException("Libavoid error: " + error.toString());
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
        error.append('.');
    }

    /** interrupt Server operation */
    public void cancelProcess() {
        synchronized (nextJob) {
            Process myProcess = process;
            if (myProcess != null) {
                libavoidStream = null;
                myProcess.destroy();
            }
        }
    }

    /** default timeout for waiting for the server to give some output. */
    public static final int PROCESS_DEF_TIMEOUT = 10000;

    /** synchronization object between the main thread and the watcher thread. */
    private Object nextJob = new Object();

    /**
     * A watcher thread that takes action when a timeout occurs.
     */
    private class Watchdog extends Thread {

        /**
         * {@inheritDoc}
         */
        @Override
        public void run() {
            do {
                synchronized (nextJob) {
                    // the watcher starts working as soon as a stream is made visible
                    while (libavoidStream == null) {
                        try {
                            // wait for notification by the main thread
                            nextJob.wait();
                        } catch (InterruptedException ex) {
                            // the watchdog thread is interrupted: shutdown is requested
                            if (watchdog != this) {
                                return;
                            }
                        }
                    }
                }

                boolean interrupted = false;
                try {
                    Thread.sleep(processTimeout);
                } catch (InterruptedException ex) {
                    // this means the main thread has done a cleanup before the timeout occurred
                    interrupted = true;
                }

                if (!interrupted) {
                    cancelProcess();
                }

            } while (watchdog == this);
        }

    }
}
