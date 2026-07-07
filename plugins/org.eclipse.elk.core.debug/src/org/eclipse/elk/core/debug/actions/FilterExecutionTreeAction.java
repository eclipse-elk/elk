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
import org.eclipse.elk.core.debug.views.AbstractLayoutDebugView;
import org.eclipse.jface.action.Action;
import org.eclipse.jface.action.IAction;

/**
 * Check box action that toggles between filtered and unfiltered logs in the layout log view.
 */
public class FilterExecutionTreeAction extends Action {
    
    /** Identifier string for this action. */
    private static final String ACTION_ID = "org.eclipse.elk.debug.actions.filter";
    /** Relative path to the icon to use for this action. */
    private static final String ICON_PATH = "icons/filter.png";
    
    /** The debug view associated with this action. */
    private final AbstractLayoutDebugView view;
    
    public FilterExecutionTreeAction(final AbstractLayoutDebugView view) {
        super("&Filter", IAction.AS_CHECK_BOX);
        setToolTipText("Display only those executions that have relevant data associated with them.");
        setId(ACTION_ID);
        setImageDescriptor(ElkDebugPlugin.imageDescriptorFromPlugin(ElkDebugPlugin.PLUGIN_ID, ICON_PATH));
        
        this.view = view;
    }
    
    
    @Override
    public void run() {
        view.setFilterTree(isChecked());
    }
}
