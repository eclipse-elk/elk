/*******************************************************************************
 * Copyright (c) 2008, 2015 Kiel University and others.
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
package org.eclipse.elk.core.ui;

import org.eclipse.core.commands.AbstractHandler;
import org.eclipse.core.commands.ExecutionEvent;
import org.eclipse.core.commands.ExecutionException;
import org.eclipse.elk.core.service.DiagramLayoutEngine;
import org.eclipse.jface.preference.IPreferenceStore;
import org.eclipse.jface.viewers.ISelection;
import org.eclipse.jface.viewers.IStructuredSelection;
import org.eclipse.ui.IEditorPart;
import org.eclipse.ui.handlers.HandlerUtil;

/**
 * Handler for execution of layout triggered by menu, toolbar, or keyboard command.
 * 
 * @author msp
 */
public class LayoutHandler extends AbstractHandler {
    
    /** preference identifier for animation of layout. */
    public static final String PREF_ANIMATION = "org.eclipse.elk.animation";
    /** preference identifier for zoom-to-fit after layout. */
    public static final String PREF_ZOOM = "org.eclipse.elk.zoomToFit";
    /** preference identifier for progress dialog. */
    public static final String PREF_PROGRESS = "org.eclipse.elk.progressDialog";
    /** parameter identifier for the scope of automatic layout. */
    public static final String PARAM_LAYOUT_SCOPE = "org.eclipse.elk.core.ui.layoutScope";
    /** value for diagram scope. */
    public static final String VAL_DIAGRAM = "diagram";
    /** value for selection scope. */
    public static final String VAL_SELECTION = "selection";

    @Override
    public Object execute(final ExecutionEvent event) throws ExecutionException {
        ISelection selection = null;

        // check parameter for layout scope, default is diagram scope
        String layoutScope = event.getParameter(PARAM_LAYOUT_SCOPE);
        if (layoutScope != null && layoutScope.equals(VAL_SELECTION)) {
            selection = HandlerUtil.getCurrentSelection(event);
        }
        
        // fetch general settings from preferences
        IPreferenceStore preferenceStore = ElkUiPlugin.getInstance().getPreferenceStore();
        boolean animation = preferenceStore.getBoolean(PREF_ANIMATION);
        boolean zoomToFit = preferenceStore.getBoolean(PREF_ZOOM);
        boolean progressDialog = preferenceStore.getBoolean(PREF_PROGRESS);

        // get the active editor, which is expected to contain the diagram for applying layout
        IEditorPart editorPart = HandlerUtil.getActiveEditor(event);
        Object diagramPart = null;
        if (selection instanceof IStructuredSelection && !selection.isEmpty()) {
            // perform layout with the given selection (only the first element is considered)
            IStructuredSelection structuredSelection = (IStructuredSelection) selection;
            if (structuredSelection.size() == 1) {
                diagramPart = structuredSelection.getFirstElement();
            } else {
                diagramPart = structuredSelection.toList();
            }
        }
        DiagramLayoutEngine.invokeLayout(editorPart, diagramPart, animation, progressDialog,
                false, zoomToFit);
        return null;
    }

}
