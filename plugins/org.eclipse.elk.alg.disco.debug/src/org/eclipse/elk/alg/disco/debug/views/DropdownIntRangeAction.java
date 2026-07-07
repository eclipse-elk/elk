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

import java.util.Observable;
import java.util.Observer;

import org.eclipse.jface.action.Action;
import org.eclipse.jface.action.ActionContributionItem;
import org.eclipse.jface.action.IMenuCreator;
import org.eclipse.swt.widgets.Control;
import org.eclipse.swt.widgets.Menu;

/**
 * Shows a drop down menu for either setting the lowest or the highest hierarchy
 * level of a graph to be drawn.
 */
public class DropdownIntRangeAction extends Action implements IMenuCreator, Observer {
    private Menu menu;
    private int lb, ub;
    private String label;
    private DropdownIntRangeAction that;
    private State bound;
    private ConfigState state;
    private DisCoGraphRenderingCanvas canvas;

    /**
     * Constructor for this drop down menu.
     * 
     * @param observeable
     *            State to be observered
     * @param label
     *            Text do be displayed alongside the menu
     * @param low
     *            Lowest option
     * @param high
     *            Highest option
     * @param initial
     *            Initially chosen option
     * @param bound
     *            Feature of the drawing to be controlled
     * @param canvas
     *            The canvas displaying the debug drawing
     */
    public DropdownIntRangeAction(final ConfigState observeable, final String label, final int low, final int high,
            final int initial, final State bound, final DisCoGraphRenderingCanvas canvas) {
        super(label + initial, Action.AS_DROP_DOWN_MENU);
        lb = low;
        ub = high;
        that = this;
        this.label = label;
        this.bound = bound;
        setMenuCreator(this);
        observeable.addObserver(this);
        state = observeable;
        this.canvas = canvas;
    }

    @Override
    public void dispose() {
        if (menu != null) {
            menu.dispose();
            menu = null;
        }
    }

    @Override
    public Menu getMenu(final Control parent) {
        if (menu != null) {
            menu.dispose();
        }

        menu = new Menu(parent);
        for (int i = lb; i <= ub; i++) {
            NumAction numAction = new NumAction("" + i);
            numAction.setCanvas(canvas);
            addToMenu(menu, numAction);
        }

        return menu;
    }

    void updateLabel(final String add) {
        this.setText(label + add);
    }

    @Override
    public Menu getMenu(final Menu parent) {
        return null;
    }

    private void addToMenu(final Menu parent, final Action action) {
        ActionContributionItem item = new ActionContributionItem(action);
        item.fill(parent, -1);
    }

    @Override
    public void update(final Observable o, final Object arg) {

        if (State.LOWER.equals(arg) && bound.equals(State.LOWER)) {
            updateLabel("" + state.getLowerLvl());
        } else {

            if (State.UPPER.equals(arg) && bound.equals(State.UPPER)) {
                updateLabel("" + state.getUpperLvl());
            } else {

                if (State.DEPTH.equals(arg)) {
                    this.ub = state.getMaxDepth() - 1;
                }
            }
        }

    }

    /**
     * An Action for each item of the drow down menu.
     * 
     * @author mic
     *
     */
    private class NumAction extends Action {
        private DisCoGraphRenderingCanvas canvas;

        NumAction(final String text) {
            super(text);
        }

        public void run() {
            String number = this.getText();
            int num = Integer.parseInt(number);
            that.updateLabel(number);
            if (that.bound.equals(State.LOWER)) {
                that.state.setLow(num);
            } else {
                that.state.setHigh(num);
            }
            canvas.redraw();

        }

        public void setCanvas(final DisCoGraphRenderingCanvas canvas) {
            this.canvas = canvas;
        }
    }
}
