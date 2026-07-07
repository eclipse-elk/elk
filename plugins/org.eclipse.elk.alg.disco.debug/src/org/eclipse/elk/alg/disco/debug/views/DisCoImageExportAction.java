/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.disco.debug.views;

import org.eclipse.core.runtime.IProgressMonitor;
import org.eclipse.core.runtime.IStatus;
import org.eclipse.core.runtime.Status;
import org.eclipse.core.runtime.jobs.Job;
import org.eclipse.elk.alg.disco.debug.DisCoDebugPlugin;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.jface.action.Action;
import org.eclipse.jface.preference.IPreferenceStore;
import org.eclipse.swt.SWT;
import org.eclipse.swt.SWTException;
import org.eclipse.swt.graphics.GC;
import org.eclipse.swt.graphics.Image;
import org.eclipse.swt.graphics.ImageData;
import org.eclipse.swt.graphics.ImageLoader;
import org.eclipse.swt.graphics.Rectangle;
import org.eclipse.swt.widgets.Display;
import org.eclipse.swt.widgets.FileDialog;

/**
 * Action that exports the currently displayed layout graph into an image file.
 */
public class DisCoImageExportAction extends Action {

    /** identifier string for this action. */
    private static final String ACTION_ID = "org.eclipse.elk.alg.disco.debug.imageExport";
    /** relative path to the icon to use for this action. */
    private static final String ICON_PATH = "icons/export.gif";
    /** preference identifier for the last used file name. */
    private static final String LAST_FILE_NAME_PREF = "disCoImageExportAction.lastImageFile";

    /** the layout graph view associated with this action. */
    private DisCoDebugView view;

    /**
     * Creates an image export action for a given layout graph view.
     * 
     * @param theview
     *            layout graph view that created this action
     */
    public DisCoImageExportAction(final DisCoDebugView theview) {
        this.view = theview;
        setId(ACTION_ID);
        setText("&Export PNG");
        setToolTipText("Export the current graph as a PNG image file.");
        setImageDescriptor(DisCoDebugPlugin.imageDescriptorFromPlugin(DisCoDebugPlugin.PLUGIN_ID, ICON_PATH));
    }

    @Override
    public void run() {
        final DisCoGraphRenderingCanvas canvas = view.getCanvas();
        if (canvas != null) {
            // let the user select an output file
            final String fileName = getFileName(canvas.getDisplay());

            if (fileName != null) {
                // create a job for painting and exporting the image
                Job saveJob = new Job("Export PNG Image") {
                    protected IStatus run(final IProgressMonitor monitor) {
                        try {
                            monitor.beginTask("Export PNG Image", 2);

                            // paint the layout graph
                            ElkNode graph = canvas.getLayoutGraph();
                            // KShapeLayout graphSize = graph.getData(KShapeLayout.class);
                            double scale = canvas.getRenderer().getScale();
                            Rectangle area = new Rectangle(0, 0, (int) ((graph.getWidth() + 1) * scale),
                                    (int) ((graph.getHeight() + 1) * scale));
                            Image image = new Image(canvas.getDisplay(), area.width, area.height);
                            canvas.getRenderer().markDirty(area);
                            GC dummyGC = new GC(image);
                            canvas.getRenderer().render(graph, dummyGC, area);
                            monitor.worked(1);
                            if (monitor.isCanceled()) {
                                return new Status(IStatus.INFO, DisCoDebugPlugin.PLUGIN_ID, 0, "Aborted", null);
                            }

                            // save the image into the selected file
                            ImageLoader imageLoader = new ImageLoader();
                            ImageData[] imageData = new ImageData[] { image.getImageData() };
                            // CHECKSTYLEOFF MagicNumber
                            imageLoader.compression = 3;
                            // CHECKSTYLEON MagicNumber

                            imageLoader.data = imageData;
                            imageLoader.save(fileName, SWT.IMAGE_PNG);
                            monitor.worked(1);

                            return new Status(IStatus.INFO, DisCoDebugPlugin.PLUGIN_ID, 0, "OK", null);
                        } catch (SWTException exception) {
                            return new Status(IStatus.ERROR, DisCoDebugPlugin.PLUGIN_ID, exception.code,
                                    "Could not save the selected PNG file.", exception);
                        } finally {
                            monitor.done();
                        }
                    }
                };

                // process the image export job
                IProgressMonitor monitor = Job.getJobManager().createProgressGroup();
                saveJob.setProgressGroup(monitor, 2);
                saveJob.setPriority(Job.SHORT);
                saveJob.setUser(true);
                saveJob.schedule();
            }
        }
    }

    /**
     * Opens a file dialog for the user to choose an output file for the exported
     * image.
     * 
     * @param display
     *            display to use as parent for the file dialog
     * @return file name, or null if no file was selected
     */
    private static String getFileName(final Display display) {
        IPreferenceStore preferenceStore = DisCoDebugPlugin.getDefault().getPreferenceStore();

        // create and configure file dialog
        FileDialog fileDialog = new FileDialog(display.getActiveShell(), SWT.SAVE);
        String[] extensions = new String[] { "*.png" };
        fileDialog.setFilterExtensions(extensions);
        String[] names = new String[] { "Portable Network Graphics files" };
        fileDialog.setFilterNames(names);
        fileDialog.setOverwrite(true);
        fileDialog.setText("Select Output File");
        fileDialog.setFileName(preferenceStore.getString(LAST_FILE_NAME_PREF));

        // open the file dialog and check the output
        String fileName = checkFileName(fileDialog.open());

        if (fileName != null) {
            preferenceStore.setValue(LAST_FILE_NAME_PREF, fileName);
        }
        return fileName;
    }

    /**
     * Checks whether the extension of the given file name is compatible with the
     * PNG format and changes it if needed.
     * 
     * @param fileName
     *            file name given by the user
     * @return file name with .png as extension
     */
    private static String checkFileName(final String fileName) {
        if (fileName == null) {
            return null;
        } else {
            String nameCopy = new String(fileName).toLowerCase();
            if (nameCopy.endsWith(".png")) {
                return fileName;
            } else if (nameCopy.endsWith(".pn")) {
                return fileName + "g";
            } else if (nameCopy.endsWith(".p")) {
                return fileName + "ng";
            } else if (nameCopy.endsWith(".")) {
                return fileName + "png";
            } else {
                return fileName + ".png";
            }
        }
    }

}
