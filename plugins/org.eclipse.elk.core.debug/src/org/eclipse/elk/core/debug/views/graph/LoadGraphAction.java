/*******************************************************************************
 * Copyright (c) 2011, 2019 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.core.debug.views.graph;

import java.io.IOException;

import org.eclipse.elk.core.IGraphLayoutEngine;
import org.eclipse.elk.core.RecursiveGraphLayoutEngine;
import org.eclipse.elk.core.debug.ElkDebugPlugin;
import org.eclipse.elk.core.debug.model.ExecutionInfo;
import org.eclipse.elk.core.debug.model.ExecutionInfoModel;
import org.eclipse.elk.core.service.DiagramLayoutEngine;
import org.eclipse.elk.core.service.ElkServicePlugin;
import org.eclipse.elk.core.util.BasicProgressMonitor;
import org.eclipse.elk.core.util.IElkProgressMonitor;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.emf.common.util.URI;
import org.eclipse.emf.common.util.WrappedException;
import org.eclipse.emf.ecore.resource.Resource;
import org.eclipse.emf.ecore.resource.ResourceSet;
import org.eclipse.emf.ecore.resource.impl.ResourceSetImpl;
import org.eclipse.jface.action.Action;
import org.eclipse.jface.preference.IPreferenceStore;
import org.eclipse.swt.SWT;
import org.eclipse.swt.widgets.Display;
import org.eclipse.swt.widgets.FileDialog;

import com.google.common.collect.Lists;

/**
 * An action for loading an ELK graph file and performing layout on it.
 */
public class LoadGraphAction extends Action {

    /** identifier string for this action. */
    private static final String ACTION_ID = "org.eclipse.elk.debug.actions.loadGraph";

    /** relative path to the icon to use for this action. */
    private static final String ICON_PATH = "icons/import.png";
    /** preference identifier for the last used file name. */
    private static final String LAST_FILE_NAME_PREF = "loadGraphAction.lastGraphFile";

    /** The layout graph view associated with this action. */
    private LayoutGraphView view;
    
    public LoadGraphAction(final LayoutGraphView theview) {
        setId(ACTION_ID);
        setText("Load Graph");
        setToolTipText("Load, layout, and display a saved ELK graph.");
        setImageDescriptor(ElkDebugPlugin.imageDescriptorFromPlugin(ElkDebugPlugin.PLUGIN_ID, ICON_PATH));
        
        this.view = theview;
    }

    @Override
    public void run() {
        IPreferenceStore prefStore = ElkDebugPlugin.getDefault().getPreferenceStore();

        // Create and configure file dialog
        FileDialog fileDialog = new FileDialog(Display.getDefault().getActiveShell(), SWT.OPEN);
        fileDialog.setFilterExtensions(new String[] { "*.elkg", "*.elkt", "*.xmi", "*.*" });
        fileDialog.setText("Select Graph File");
        fileDialog.setFileName(prefStore.getString(LAST_FILE_NAME_PREF));

        // Open the file dialog and wait for file name
        String fileName = fileDialog.open();
        run(fileName, view);
    }
    
    /**
     * Load an ELK graph file and performing layout on it.
     * 
     * @param fullFilePath The full path of the file to load.
     * @param layoutGraphView The Layout Graph view
     */
    public static void run(final String fullFilePath, final LayoutGraphView layoutGraphView) {
        IPreferenceStore prefStore = ElkDebugPlugin.getDefault().getPreferenceStore();
        if (fullFilePath != null) {
            prefStore.setValue(LAST_FILE_NAME_PREF, fullFilePath);

            // Load the file content
            try {
                ElkNode content = loadFromFile(fullFilePath);
                ExecutionInfo info = layout(fullFilePath, content);
                ElkDebugPlugin.getDefault().getModel().addExecution(info);
                // Select the new loaded element (made in async to let the view refresh before)
                layoutGraphView.getSite().getShell().getDisplay().asyncExec(new Runnable() {
                    @Override
                    public void run() {
                        layoutGraphView.setSelectedExecutionInfos(Lists.newArrayList(info));
                    }
                });
            } catch (IOException exception) {
                throw new WrappedException(exception);
            }
        }
    }

    /**
     * Loads an ELK graph from the given file.
     * 
     * @param fileName
     *            path to the file to load.
     * @return {@link ElkNode} loaded from the file.
     * @throws IOException
     *             if anything goes wrong.
     */
    static ElkNode loadFromFile(final String fileName) throws IOException {
        // Load the file content
        ResourceSet resourceSet = new ResourceSetImpl();
        URI uri = URI.createFileURI(fileName);

        Resource resource = resourceSet.createResource(uri);
        resource.load(null);

        return (ElkNode) resource.getContents().get(0);
    }

    /**
     * Perform layout, if activated, and return an {@link ExecutionInfo} that contains information about the layout run.
     * That object is not added to the {@link ExecutionInfoModel} by this method.
     */
    static ExecutionInfo layout(final String fileName, final ElkNode graph) {
        return layout(fileName, LayoutUponLoadSettingAction.shouldLayoutUponLoad(), graph);
    }

    /**
     * Perform layout, if requested, and return an {@link ExecutionInfo} that contains information about the layout run.
     * That object is not added to the {@link ExecutionInfoModel} by this method.
     */
    static ExecutionInfo layout(final String fileName, final boolean performLayout, final ElkNode graph) {
        IPreferenceStore prefStore = ElkServicePlugin.getInstance().getPreferenceStore();
        IElkProgressMonitor monitor =
                new BasicProgressMonitor()
                        .withMaxHierarchyLevels(0)
                        .withLogging(prefStore.getBoolean(DiagramLayoutEngine.PREF_DEBUG_LOGGING))
                        .withLogPersistence(prefStore.getBoolean(DiagramLayoutEngine.PREF_DEBUG_STORE))
                        .withExecutionTimeMeasurement(prefStore.getBoolean(DiagramLayoutEngine.PREF_DEBUG_EXEC_TIME));
        monitor.begin(fileName, 1);

        // Perform layout using a graph layout engine, if enabled
        if (performLayout) {
            IGraphLayoutEngine layoutEngine = new RecursiveGraphLayoutEngine();
            layoutEngine.layout(graph, monitor.subTask(1));
        }

        monitor.done();

        // We're not going through the DiagramLayoutEngine, but directly through the RecursiveGraphLayoutEngine, which
        // means that not layout events will be fired. We'll have to update our model manually.
        monitor.logGraph(graph, "Result");

        return ExecutionInfo.fromProgressMonitorAndFile(monitor, fileName, performLayout);
    }

}
