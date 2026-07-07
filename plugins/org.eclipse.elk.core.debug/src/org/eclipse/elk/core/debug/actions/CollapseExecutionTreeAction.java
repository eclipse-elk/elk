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

import org.eclipse.elk.core.debug.views.AbstractLayoutDebugView;
import org.eclipse.jface.action.Action;
import org.eclipse.ui.ISharedImages;
import org.eclipse.ui.PlatformUI;

/**
 * Action for collapsing all elements of the tree view in the debug log view.
 */
public class CollapseExecutionTreeAction extends Action {
    
    /** Identifier string for this action. */
    private static final String ACTION_ID = "org.eclipse.elk.debug.actions.collapseAll";
    
    /** The debug view associated with this action. */
    private final AbstractLayoutDebugView view;

    public CollapseExecutionTreeAction(final AbstractLayoutDebugView view) {
        setId(ACTION_ID);
        setText("&Collapse All");
        setToolTipText("Collapses every element of the tree viewer.");
        setImageDescriptor(PlatformUI.getWorkbench().getSharedImages().getImageDescriptor(
                ISharedImages.IMG_ELCL_COLLAPSEALL));
        
        this.view = view;
    }

    @Override
    public void run() {
        view.collapseAllTreeViewerElements();
    }
}
