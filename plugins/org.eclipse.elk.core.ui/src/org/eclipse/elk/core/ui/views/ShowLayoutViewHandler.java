/*******************************************************************************
 * Copyright (c) 2011, 2015 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.core.ui.views;

import org.eclipse.core.commands.AbstractHandler;
import org.eclipse.core.commands.ExecutionEvent;
import org.eclipse.core.commands.ExecutionException;
import org.eclipse.core.runtime.IStatus;
import org.eclipse.core.runtime.Status;
import org.eclipse.elk.core.ui.ElkUiPlugin;
import org.eclipse.ui.IWorkbenchWindow;
import org.eclipse.ui.PartInitException;
import org.eclipse.ui.handlers.HandlerUtil;
import org.eclipse.ui.statushandlers.StatusManager;

/**
 * A command handler for displaying the layout view.
 *
 * @author msp
 */
public class ShowLayoutViewHandler extends AbstractHandler {

    @Override
    public Object execute(final ExecutionEvent event) throws ExecutionException {
        IWorkbenchWindow workbenchWindow = HandlerUtil.getActiveWorkbenchWindow(event);
        if (workbenchWindow != null) {
            try {
                workbenchWindow.getActivePage().showView(LayoutViewPart.VIEW_ID);
            } catch (PartInitException exception) {
                IStatus status = new Status(IStatus.ERROR, ElkUiPlugin.PLUGIN_ID,
                        "Could not open Layout View.", exception);
                StatusManager.getManager().handle(status, StatusManager.SHOW);
            }
        }
        return null;
    }

}
