/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
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
package org.eclipse.elk.alg.disco.debug.views;

import org.eclipse.jface.action.Action;

/**
 * A button for showing/hiding different aspects of the drawing in the
 * {@link DisCoDebugView}, controlled via the {@link DisCoGraphRenderer}.
 */
public class SelectivePaintToggleAction extends Action {
    /** identifier string for this action. */
    private static final String ACTION_ID = "org.eclipse.elk.debug.disco.selectivepainttoggle";
    /** the canvas to draw on. Used for enforcing a redraw */
    private DisCoGraphRenderingCanvas canvas;
    /**
     * the renderer providing a method for enabling/disabling the visibility of
     * certain layers of the drawing.
     */
    private DisCoGraphRenderer renderer;
    /**
     * a command for the method
     * {@link DisCoGraphRenderer#setVisible(String, boolean)} (see method
     * description for supported commands).
     */
    private String command;
    /**
     * Indicating whether the button press is supposed to show (or hide) layers of
     * the drawing.
     */
    private boolean val;

    /**
     * Makes a new button for executing a specific toggle command.
     * 
     * @param command
     *            a command for the method
     *            {@link DisCoGraphRenderer#setVisible(String, boolean)} (see method
     *            description for supported commands)
     * @param title
     *            text visible on the button in the ui
     * @param tooltip
     *            text appearing when hovering over the button for some time
     * @param canvas
     *            canvas to influence
     */
    public SelectivePaintToggleAction(final String command, final String title, final String tooltip,
            final DisCoGraphRenderingCanvas canvas) {
        super(title, Action.AS_CHECK_BOX);
        this.canvas = canvas;
        this.renderer = canvas.getRenderer();
        this.command = command;
        this.val = false;
        setChecked(true);
        setId(ACTION_ID + command);
        setToolTipText(tooltip);
    }

    @Override
    public void run() {
        renderer.setVisible(command, val);
        // toggle between true and false
        setChecked(val);
        val = !val;
        this.canvas.redraw();
    }
}
