/*******************************************************************************
 * Copyright (c) 2009, 2019 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.core.debug.views.time;

import java.util.function.Predicate;

import org.eclipse.elk.core.debug.model.ExecutionInfo;
import org.eclipse.elk.core.debug.views.AbstractLayoutDebugView;
import org.eclipse.elk.core.service.DiagramLayoutEngine;
import org.eclipse.jface.dialogs.IDialogSettings;
import org.eclipse.jface.viewers.DelegatingStyledCellLabelProvider;
import org.eclipse.jface.viewers.TreeViewer;
import org.eclipse.jface.viewers.TreeViewerColumn;
import org.eclipse.swt.SWT;

/**
 * View that can be used to display execution results of an algorithm.
 */
public class ExecutionTimeView extends AbstractLayoutDebugView {

    /** The view identifier. */
    public static final String VIEW_ID = "org.eclipse.elk.debug.executionView";
    
    // UI Controls
    private TreeViewerColumn nameColumn;
    private TreeViewerColumn timeColumn;
    private TreeViewerColumn localTimeColumn;
    
    /**
     * Creates an execution view.
     */
    public ExecutionTimeView() {
        super(VIEW_ID, DiagramLayoutEngine.PREF_DEBUG_EXEC_TIME);
    }
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Execution Tree Handling

    @Override
    protected Predicate<ExecutionInfo> getTreeFilter() {
        return info -> info.isExecutionTimeMeasured();
    }
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // UI Setup

    @Override
    protected void customizeTreeViewer(TreeViewer treeViewer) {
        treeViewer.getTree().setHeaderVisible(true);
        
        nameColumn = new TreeViewerColumn(treeViewer, SWT.NONE);
        nameColumn.setLabelProvider(new DelegatingStyledCellLabelProvider(
                new ExecutionTimeLabelProvider(ExecutionTimeLabelProvider.DisplayMode.NAME)));
        nameColumn.getColumn().setText("Name");
        nameColumn.getColumn().setWidth(500);
        
        timeColumn = new TreeViewerColumn(treeViewer, SWT.NONE);
        timeColumn.setLabelProvider(new DelegatingStyledCellLabelProvider(
                new ExecutionTimeLabelProvider(ExecutionTimeLabelProvider.DisplayMode.TIME_TOTAL)));
        timeColumn.getColumn().setAlignment(SWT.RIGHT);
        timeColumn.getColumn().setText("Time [ms]");
        timeColumn.getColumn().setWidth(100);
        
        localTimeColumn = new TreeViewerColumn(treeViewer, SWT.NONE);
        localTimeColumn.setLabelProvider(new DelegatingStyledCellLabelProvider(
                new ExecutionTimeLabelProvider(ExecutionTimeLabelProvider.DisplayMode.TIME_LOCAL)));
        localTimeColumn.getColumn().setAlignment(SWT.RIGHT);
        localTimeColumn.getColumn().setText("Local Time [ms]");
        localTimeColumn.getColumn().setWidth(100);
    }
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // State Persistence
    
    private static final String COL_NAME_WIDTH_KEY = "nameColumn.width";
    private static final String COL_TIME_WIDTH_KEY = "timeColumn.width";
    private static final String COL_LOCAL_TIME_WIDTH_KEY = "localTimeColumn.width";
    
    @Override
    protected void loadState(IDialogSettings viewSettings) {
        super.loadState(viewSettings);
        
        loadColumnWidth(viewSettings, nameColumn, COL_NAME_WIDTH_KEY);
        loadColumnWidth(viewSettings, timeColumn, COL_TIME_WIDTH_KEY);
        loadColumnWidth(viewSettings, localTimeColumn, COL_LOCAL_TIME_WIDTH_KEY);
    }
    
    private void loadColumnWidth(IDialogSettings viewSettings, TreeViewerColumn column, String settingKey) {
        if (viewSettings.get(settingKey) != null) {
            column.getColumn().setWidth(viewSettings.getInt(settingKey));
        }
    }

    @Override
    protected void persistState(IDialogSettings viewSettings) {
        super.persistState(viewSettings);
        
        viewSettings.put(COL_NAME_WIDTH_KEY, nameColumn.getColumn().getWidth());
        viewSettings.put(COL_TIME_WIDTH_KEY, timeColumn.getColumn().getWidth());
        viewSettings.put(COL_LOCAL_TIME_WIDTH_KEY, localTimeColumn.getColumn().getWidth());
    }

}
