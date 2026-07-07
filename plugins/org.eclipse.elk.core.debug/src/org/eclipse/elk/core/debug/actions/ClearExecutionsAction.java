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

import org.eclipse.elk.core.debug.ElkDebugPlugin;
import org.eclipse.jface.action.Action;
import org.eclipse.ui.ISharedImages;
import org.eclipse.ui.PlatformUI;

/**
 * Action that clears all execution infos. The debug views are automatically updated since they listen for model
 * changes.
 */
public class ClearExecutionsAction extends Action {

    /** Identifier string for this action. */
    private static final String ACTION_ID = "org.eclipse.elk.debug.actions.clearExecutions";
    
    public ClearExecutionsAction() {
        setId(ACTION_ID);
        setText("&Remove All");
        setToolTipText("Removes all layout executions.");
        setImageDescriptor(PlatformUI.getWorkbench().getSharedImages().getImageDescriptor(
                ISharedImages.IMG_ELCL_REMOVEALL));
    }
    
    /**
     * Updates the enabled property of this action.
     */
    public void updateEnablement() {
        setEnabled(!ElkDebugPlugin.getDefault().getModel().getExecutionInfos().isEmpty());
    }

    @Override
    public void run() {
        ElkDebugPlugin.getDefault().getModel().removeAllExecutionInfos();
    }

}
