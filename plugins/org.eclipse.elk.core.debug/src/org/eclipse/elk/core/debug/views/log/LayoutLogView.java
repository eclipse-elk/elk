/*******************************************************************************
 * Copyright (c) 2019 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.core.debug.views.log;

import java.util.function.Predicate;

import org.eclipse.elk.core.debug.model.ExecutionInfo;
import org.eclipse.elk.core.debug.views.AbstractLayoutDebugView;
import org.eclipse.elk.core.service.DiagramLayoutEngine;
import org.eclipse.jface.resource.JFaceResources;
import org.eclipse.jface.viewers.DelegatingStyledCellLabelProvider;
import org.eclipse.jface.viewers.TreeViewer;
import org.eclipse.swt.SWT;
import org.eclipse.swt.widgets.Composite;
import org.eclipse.swt.widgets.Text;

/**
 * This view displays the debug logs of an layout algorithm.
 */
public class LayoutLogView extends AbstractLayoutDebugView {

    /** The view identifier. */
    public static final String VIEW_ID = "org.eclipse.elk.core.debug.logView";

    // UI Controls
    private Text logViewer;
    
    public LayoutLogView() {
        super(VIEW_ID, DiagramLayoutEngine.PREF_DEBUG_LOGGING);
    }

    @Override
    protected void treeSelectionChanged() {
        super.treeSelectionChanged();
        
        String text = "";

        for (ExecutionInfo info : getSelectedExecutionInfos()) {
            if (!info.getLogMessages().isEmpty()) {
                text += String.join("\n", info.getLogMessages()) + "\n\n";
            }
        }

        logViewer.setText(text);
    }

    @Override
    protected Predicate<ExecutionInfo> getTreeFilter() {
        return info -> info.hasLogMessages() || info.hasDescendantsWithLogMessages();
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    // User Interface

    @Override
    protected void customizeTreeViewer(TreeViewer treeViewer) {
        treeViewer.setLabelProvider(new DelegatingStyledCellLabelProvider(new LogLabelProvider()));
    }

    @Override
    protected void setupRemainingControls(Composite parent) {
        logViewer = new Text(parent, SWT.H_SCROLL | SWT.V_SCROLL);
        logViewer.setFont(JFaceResources.getFont(JFaceResources.TEXT_FONT));
        logViewer.setEditable(false);
    }

}
