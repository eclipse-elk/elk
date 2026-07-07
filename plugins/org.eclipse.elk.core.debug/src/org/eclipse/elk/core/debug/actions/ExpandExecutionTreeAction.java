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

/**
 * Action for expanding all elements in the tree view of the debug log view. 
 */
public class ExpandExecutionTreeAction extends Action {
    
    /** Identifier string for this action. */
    private static final String ACTION_ID = "org.eclipse.elk.debug.actions.expandAll";
    /** Relative path to the icon to use for this action. */
    private static final String ICON_PATH = "icons/expand_all.gif";
    
    /** The debug view associated with this action. */
    private final AbstractLayoutDebugView view;
    
    public ExpandExecutionTreeAction(final AbstractLayoutDebugView view) {
        setId(ACTION_ID);
        setText("&Expand All");
        setToolTipText("Expands every element of the tree viewer.");
        setImageDescriptor(ElkDebugPlugin.imageDescriptorFromPlugin(ElkDebugPlugin.PLUGIN_ID, ICON_PATH));
        
        this.view = view;
    }

    @Override
    public void run() {
        view.expandAllTreeViewerElements();
    }

}
