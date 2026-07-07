/*******************************************************************************
 * Copyright (c) 2019 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.core.debug.views.graph;

import java.io.IOException;
import java.util.List;

import org.eclipse.elk.core.debug.ElkDebugPlugin;
import org.eclipse.elk.core.debug.model.ExecutionInfo;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.emf.common.util.WrappedException;
import org.eclipse.jface.action.Action;

import com.google.common.collect.Lists;

/**
 * Reloads a graph from the file it came from.
 */
public class ReloadFromFileAction extends Action {

    /** identifier string for this action. */
    private static final String ACTION_ID = "org.eclipse.elk.debug.actions.reloadFromFile";
    
    /** relative path to the icon to use for this action. */
    private static final String ICON_PATH = "icons/refresh.gif";

    /** The layout graph view associated with this action. */
    private LayoutGraphView view;

    /**
     * Creates a reload action for a given layout graph view.
     */
    public ReloadFromFileAction(final LayoutGraphView theview) {
        setId(ACTION_ID);
        setText("&Reload from File");
        setToolTipText("Reloads a selected graph from .");
        setImageDescriptor(ElkDebugPlugin.imageDescriptorFromPlugin(ElkDebugPlugin.PLUGIN_ID, ICON_PATH));
        
        this.view = theview;
    }

    /**
     * Updates the enabled property of this action.
     */
    public void updateEnablement() {
        List<ExecutionInfo> selection = view.getSelectedExecutionInfos();
        this.setEnabled(selection.size() == 1 && selection.get(0).isLoadedFromFile());
    }

    @Override
    public void run() {
        List<ExecutionInfo> selection = view.getSelectedExecutionInfos();
        
        // We don't do anything unless a single info is selected which was loaded from a file
        if (selection.size() != 1 || !selection.get(0).isLoadedFromFile()) {
            return;
        }
        
        ExecutionInfo oldInfo = selection.get(0);
        
        // Load the graph again, if possible
        try {
            ElkNode graph = LoadGraphAction.loadFromFile(oldInfo.getFileName());
            ExecutionInfo newInfo = LoadGraphAction.layout(oldInfo.getFileName(), oldInfo.isLaidOutAfterLoad(), graph);
            ElkDebugPlugin.getDefault().getModel().replaceExecution(oldInfo, newInfo);
            // We change selection in async to let the view refresh before
            view.getSite().getShell().getDisplay().asyncExec(new Runnable() {
                @Override
                public void run() {
                    view.setSelectedExecutionInfos(Lists.newArrayList(newInfo));
                }
            });
            
        } catch (IOException e) {
            throw new WrappedException(e);
        }
    }

}
