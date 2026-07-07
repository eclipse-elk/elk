/*******************************************************************************
 * Copyright (c) 2019 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.core.debug.actions;

import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Iterator;
import java.util.List;
import java.util.stream.Collectors;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import org.eclipse.core.runtime.IProgressMonitor;
import org.eclipse.core.runtime.IStatus;
import org.eclipse.core.runtime.Status;
import org.eclipse.core.runtime.jobs.Job;
import org.eclipse.elk.core.debug.ElkDebugPlugin;
import org.eclipse.elk.core.debug.model.ExecutionInfo;
import org.eclipse.elk.core.debug.views.AbstractLayoutDebugView;
import org.eclipse.jface.action.Action;
import org.eclipse.swt.SWT;
import org.eclipse.swt.widgets.Display;
import org.eclipse.swt.widgets.FileDialog;

/**
 * Action that allows users to produce a zip file containing the logs of selected layout runs.
 */
public class CompressLogFolderAction extends Action {

    /** identifier string for this action. */
    private static final String ACTION_ID = "org.eclipse.elk.debug.actions.compressLogFolder";
    
    /** The layout graph view associated with this action. */
    private AbstractLayoutDebugView view;

    public CompressLogFolderAction(final AbstractLayoutDebugView theview) {
        setId(ACTION_ID);
        setText("&Compress Log Folder");
        setToolTipText("Produces a zip archive for the selected log folder.");
        
        this.view = theview;
    }
    
    /**
     * Updates the enabled property of this action.
     */
    public void updateEnablement() {
        List<ExecutionInfo> selectedLayoutRuns = getSelectedLayoutRunsWithDebugFolder();
        setEnabled(!selectedLayoutRuns.isEmpty());
    }

    @Override
    public void run() {
        final List<ExecutionInfo> selectedInfos = getSelectedLayoutRunsWithDebugFolder();
        if (!selectedInfos.isEmpty()) {
            // Ask the user for the name and location of the target file
            final String fileName = getFileName(selectedInfos, view.getSite().getShell().getDisplay());
            if (fileName != null) {
                // This job will take care of producing the zip file
                Job saveJob = new Job("Compress Log Folder") {
                    protected IStatus run(IProgressMonitor monitor) {
                        return doZip(selectedInfos, fileName, monitor);
                    }
                };
                
                // Start the job
                IProgressMonitor monitor = Job.getJobManager().createProgressGroup();
                saveJob.setProgressGroup(monitor, selectedInfos.size());
                saveJob.setPriority(Job.LONG);
                saveJob.setUser(true);
                saveJob.schedule();
            }
        }
    }
    
    private List<ExecutionInfo> getSelectedLayoutRunsWithDebugFolder() {
        // Iterate over the execution infos and filter out those that are not top-level elements
        return view.getSelectedExecutionInfos().stream()
                .filter(info -> info.getParent() == null)
                .filter(info -> info.getDebugFolder() != null)
                .filter(info -> Files.isDirectory(info.getDebugFolder()))
                .collect(Collectors.toList());
    }
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // File Stuff

    /**
     * Opens a file dialog for the user to choose an output file for the
     * exported image.
     * 
     * @param display display to use as parent for the file dialog
     * @return file name, or null if no file was selected
     */
    private static String getFileName(final List<ExecutionInfo> executionInfos, final Display display) {
        // the suggested file name consists of the debug folder names
        String suggestedFileName = executionInfos.stream()
            .map(info -> info.getDebugFolder())
            .filter(folder -> folder != null)
            .map(folder -> folder.getFileName().toString())
            .collect(Collectors.joining(", "))
            + ".zip";
        
        // create and configure file dialog
        String[] extensions = new String[] {"*.zip"};
        String[] names = new String[] {"ZIP archives"};
        
        FileDialog fileDialog = new FileDialog(display.getActiveShell(), SWT.SAVE);
        fileDialog.setFilterExtensions(extensions);
        fileDialog.setFilterNames(names);
        fileDialog.setOverwrite(true);
        fileDialog.setText("Select Output File");
        fileDialog.setFileName(suggestedFileName);

        // open the file dialog and check the output
        String fileName = completeFileName(fileDialog.open());
        
        return fileName;
    }

    /**
     * Checks whether the extension of the given file name is compatible with
     * the PNG format and changes it if needed.
     * 
     * @param fileName file name given by the user
     * @return file name with .png as extension
     */
    private static String completeFileName(final String fileName) {
        if (fileName == null) {
            return null;
        } else {
            String nameCopy = new String(fileName).toLowerCase();
            if (nameCopy.endsWith(".zip")) {
                return fileName;
            } else if (nameCopy.endsWith(".zi")) {
                return fileName + "p";
            } else if (nameCopy.endsWith(".z")) {
                return fileName + "ip";
            } else if (nameCopy.endsWith(".")) {
                return fileName + "zip";
            } else {
                return fileName + ".zip";
            }
        }
    }
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Zip Code (hrhr)
    
    /**
     * Zips the folders that belong to the given layout runs and saves the result in the ZIP file with the given name.
     * 
     * @param layoutRuns
     *            the layout runs to be compressed.
     * @param zipFileName
     *            name of the resulting zip file.
     * @param monitor
     *            monitor to report progress to.
     * @return a status code.
     */
    private IStatus doZip(final List<ExecutionInfo> layoutRuns, final String zipFileName,
            final IProgressMonitor monitor) {
        
        FileOutputStream outFileStream = null;
        ZipOutputStream zipFileStream = null;
        
        try {
            // Create an output stream for the zip file
            outFileStream = new FileOutputStream(zipFileName);
            zipFileStream = new ZipOutputStream(outFileStream);
            
            monitor.beginTask("Compressing log folders to " + zipFileName, layoutRuns.size());
            
            // Add an entry to represent the zip file's root folder
            zipFileStream.putNextEntry(new ZipEntry("/"));
            zipFileStream.closeEntry();
            
            // Process each execution info
            for (ExecutionInfo layoutRun : layoutRuns) {
                // Check for cancellation
                if (monitor.isCanceled()) {
                    return Status.CANCEL_STATUS;
                } else {
                    doZipLayoutRun(layoutRun, zipFileStream, monitor);
                }
            }
            
            // Return successful status code (note that the output stream is closed in the finally block below)
            monitor.done();
            return Status.OK_STATUS;
            
        } catch (Exception e) {
            return new Status(
                    IStatus.ERROR,
                    ElkDebugPlugin.PLUGIN_ID,
                    "Unable to compress log folders to " + zipFileName + ".",
                    e);
            
        } finally {
            if (zipFileStream != null) {
                try {
                    zipFileStream.close();
                } catch (IOException e) {
                    return new Status(
                            IStatus.ERROR,
                            ElkDebugPlugin.PLUGIN_ID,
                            "Unable to compress log folders to " + zipFileName + ".",
                            e);
                }
            }
        }
    }
    
    /**
     * Adds the content of the given layout run's directory tree to the zip file stream.
     * 
     * @param layoutRun
     *            the layout run whose debug folder to add.
     * @param zipFileStream
     *            the zip file stream to add the content to.
     * @param monitor
     *            progress monitor to report progress to. This method consumes one unit of work.
     * @throws Exception
     *             if anything goes wrong during file operations.
     */
    private void doZipLayoutRun(final ExecutionInfo layoutRun, final ZipOutputStream zipFileStream,
            final IProgressMonitor monitor) throws Exception {
        
        monitor.subTask("Compressing " + layoutRun.getDebugFolder().getFileName());
        
        Path debugFolder = layoutRun.getDebugFolder();
        if (Files.isDirectory(debugFolder)) {
            doZipFolder("/", debugFolder, zipFileStream);
        }
        
        monitor.worked(1);
    }

    /**
     * Adds the content of the given directory tree to the zip file stream. The folder's name is appended to the path
     * prefix to produce zip entry names (which is how we produce path names relative to each layout run's root debug
     * folder).
     * 
     * @param pathPrefix
     *            prefix to prefix the folder's name with when adding it to the zip file.
     * @param folder
     *            root of the directory tree to add.
     * @param zipFileStream
     *            the zip file stream to add the content to.
     * @throws Exception
     *             if anything goes wrong during file operations.
     */
    private void doZipFolder(final String pathPrefix, final Path folder, ZipOutputStream zipFileStream)
            throws Exception {
        
        // Compute the new path prefix for entries (this also doubles as the zip entry's name)
        String newPathPrefix = pathPrefix + folder.getFileName() + "/";
        
        // Add an entry to represent this folder
        zipFileStream.putNextEntry(new ZipEntry(newPathPrefix));
        zipFileStream.closeEntry();
        
        // Iterate over the folder's children and zip them as well
        Iterator<Path> childIterator = Files.list(folder).iterator();
        while (childIterator.hasNext()) {
            Path child = childIterator.next();
            
            if (Files.isDirectory(child)) {
                doZipFolder(newPathPrefix, child, zipFileStream);
                
            } else if (Files.isRegularFile(child)) {
                doZipFile(newPathPrefix, child, zipFileStream);
            }
        }
    }
    
    /**
     * Adds the content of the given fileto the zip file stream. The file's name is appended to the path prefix to
     * produce zip entry names (which is how we produce path names relative to each layout run's root debug folder).
     * 
     * @param pathPrefix
     *            prefix to prefix the folder's name with when adding it to the zip file.
     * @param file
     *            file to add.
     * @param zipFileStream
     *            the zip file stream to add the content to.
     * @throws Exception
     *             if anything goes wrong during file operations.
     */
    private void doZipFile(final String pathPrefix, final Path file, ZipOutputStream zipFileStream) throws Exception {
        // Add an entry to represent this file
        zipFileStream.putNextEntry(new ZipEntry(pathPrefix + file.getFileName()));
        
        // Copy the file's content over to the zip file
        Files.copy(file, zipFileStream);
        
        // Close our zip file's entry
        zipFileStream.closeEntry();
    }

}
