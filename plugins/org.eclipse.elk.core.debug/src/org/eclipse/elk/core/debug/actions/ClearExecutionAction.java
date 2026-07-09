/*******************************************************************************
 * Copyright (c) 2019 Kiel University and others.
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
package org.eclipse.elk.core.debug.actions;

import org.eclipse.elk.core.debug.ElkDebugPlugin;
import org.eclipse.elk.core.debug.model.ExecutionInfo;
import org.eclipse.elk.core.debug.views.AbstractLayoutDebugView;
import org.eclipse.jface.action.Action;
import org.eclipse.ui.ISharedImages;
import org.eclipse.ui.PlatformUI;

/**
 * Action that clears the currently selected execution infos. The debug views are automatically updated since they
 * listen for model changes.
 */
public class ClearExecutionAction extends Action {

    /** Identifier string for this action. */
    private static final String ACTION_ID = "org.eclipse.elk.debug.actions.clearExecution";
    
    /** The debug view associated with this action. */
    private final AbstractLayoutDebugView view;
    
    public ClearExecutionAction(final AbstractLayoutDebugView view) {
        setId(ACTION_ID);
        setText("&Remove");
        setToolTipText("Removes selected layout executions.");
        setImageDescriptor(PlatformUI.getWorkbench().getSharedImages().getImageDescriptor(
                ISharedImages.IMG_ELCL_REMOVE));
        
        this.view = view;
    }
    
    /**
     * Updates the enabled property of this action.
     */
    public void updateEnablement() {
        setEnabled(!view.getSelectedExecutionInfos().isEmpty());
    }

    @Override
    public void run() {
        ExecutionInfo[] selectedInfos = view.getSelectedExecutionInfos().toArray(new ExecutionInfo[0]);
        ElkDebugPlugin.getDefault().getModel().removeExecutionInfos(selectedInfos);
    }

}
