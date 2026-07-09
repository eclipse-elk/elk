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
package org.eclipse.elk.core.debug.views.graph;

import org.eclipse.elk.core.debug.ElkDebugPlugin;
import org.eclipse.jface.action.Action;
import org.eclipse.jface.action.IAction;
import org.eclipse.jface.dialogs.IDialogSettings;

/**
 * Action that controls whether layout is applied when a graph is loaded.
 */
public class LayoutUponLoadSettingAction extends Action {

    /** identifier string for this action. */
    private static final String ACTION_ID = "org.eclipse.elk.debug.actions.layoutUponLoad";

    public LayoutUponLoadSettingAction() {
        super("&Apply layout to loaded graph files", IAction.AS_CHECK_BOX);
        setId(ACTION_ID);
        setToolTipText("When active, automatic layout is applied to graphs loaded from files. "
                + "Otherwise, the graphs are drawn as they are.");
        
        // Initialize from dialog settings
        setChecked(shouldLayoutUponLoad());
    }

    @Override
    public void run() {
        // Save to dialog settings
        ElkDebugPlugin.getDefault().getDialogSettings().put(ACTION_ID, isChecked());
    }
    
    /**
     * Returns the current value of the setting. Defaults to {@code true}.
     */
    public static boolean shouldLayoutUponLoad() {
        IDialogSettings dialogSettings = ElkDebugPlugin.getDefault().getDialogSettings();
        if (dialogSettings.get(ACTION_ID) == null) {
            return true;
        } else {
            return dialogSettings.getBoolean(ACTION_ID);
        }
    }

}
