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

/**
 * Main class for configuring the visuals of the debug view.
 */
public class ConfigState extends Observable {

    // Adjustable via UI buttons in the running instance
    private int lowerLvl = 0;
    private int upperLvl = 0;

    // Max hierarchy level of current graph
    private Integer maxDepth;

    // For some properties for png export and specific graphs, there is no
    // corresponding ui-element to set them, just change these constants. Most of
    // them are intended for better legibility when using debug output in a paper.
    // UI might be implemented in the future, if needed.
    private static final double SCALE = 1.0;
    private static final boolean GHOST_PARENT = false;
    private static final BaseGridVisibility GRID = BaseGridVisibility.SIMPLE;
    private static final boolean MARK_CENTER = false;
    private static final boolean DRAW_POLYLINES_BLACK = false;
    private static final boolean DRAW_LABELS = true;
    private static final boolean MAKE_PATTERNS_SOLID = false;
    private static final double LINE_THICKNESS_SCALE = 1.0;
    private static final int PATTERN_SCALE = 1;
    private static final boolean REMOVE_LEVEL_NUMBER = false;

    /**
     * Sets the lowest hierarchy level to be drawn by the view.
     * 
     * @param low
     *            Lowest hierarchy level to be drawn
     */
    void setLow(final int low) {
        lowerLvl = low;
        setChanged();
        notifyObservers(State.LOWER);
        if (lowerLvl > upperLvl) {
            upperLvl = lowerLvl;
            setChanged();
            notifyObservers(State.UPPER);
        }
    }

    /**
     * Sets the highest hierarchy level to be drawn by the view.
     * 
     * @param high
     *            Highest hierarchy level to be drawn
     */
    void setHigh(final int high) {
        upperLvl = high;
        setChanged();
        notifyObservers(State.UPPER);
        if (lowerLvl > upperLvl) {
            lowerLvl = upperLvl;
            setChanged();
            notifyObservers(State.LOWER);
        }
    }

    /**
     * Sets the max hierarchy level of the current graph.
     * 
     * @param depth
     *            max hierarchy level of the current graph.
     */
    void setDepth(final int depth) {
        maxDepth = depth;
        setChanged();
        notifyObservers(State.DEPTH);
    }

    /**
     * Returns the lowest hierarchy level to be drawn.
     * 
     * @return the lowerLvl
     */
    public int getLowerLvl() {
        return lowerLvl;
    }

    /**
     * Returns the highest hierarchy level to be drawn.
     * 
     * @return the upperLvl
     */
    public int getUpperLvl() {
        return upperLvl;
    }

    /**
     * Returns the max hierarchy level of the current graph.
     * 
     * @return the maxDepth
     */
    public Integer getMaxDepth() {
        return maxDepth;
    }

    /**
     * Returns the linear scaling factor the image is to be drawn with.
     * 
     * @return the scale
     */
    public double getScale() {
        return SCALE;
    }

    /**
     * Returns whether the highest level of the graph just below the root is to be
     * drawn as a reference (only a line drawing without filling).
     * 
     * @return true, draw it; false, otherwise
     */
    public boolean paintParentTransparently() {
        return drawGhostParent();
    }

    /**
     * Returns the option for the base grid visibilty.
     * 
     * @return visibility
     */
    public BaseGridVisibility howToPaintGrid() {
        return GRID;
    }

    /**
     * Returns whether the centers of polyominoes' are to be hihglightet in the
     * view.
     * 
     * @return true, if yes; false, otherwise
     */
    public boolean markTheCenter() {
        return MARK_CENTER;
    }

    /**
     * Should the lines of a polyomino be drawn in solid black (good for printing on
     * paper).
     * 
     * @return true, yes; false, otherwise.
     */
    public boolean drawPolyLinesBlack() {
        return DRAW_POLYLINES_BLACK;
    }

    /**
     * Should each element of the view be drawn with its debug id as a label?
     * 
     * @return true, yes; false, otherwise.
     */
    public boolean drawLabels() {
        return DRAW_LABELS;
    }

    /**
     * Should extension patterns be drawn in a solid pattern (good for printing on
     * paper).
     * 
     * @return true, yes; false, otherwise.
     */
    public boolean makeSolid() {
        return MAKE_PATTERNS_SOLID;
    }

    /**
     * Returns the thickness of drawn lines in pixels.
     * 
     * @return the thickness.
     */
    public double getThicknessScale() {
        return LINE_THICKNESS_SCALE;
    }

    /**
     * Returns the scaling factor for the fill pattern of extenstions.
     * 
     * @return the scaling factor.
     */
    public int getPatternScale() {
        return PATTERN_SCALE;
    }

    /**
     * /** Returns whether the highest level of the graph just below the root is to
     * be drawn as a reference (only a line drawing without filling).
     * 
     * @return true, draw it; false, otherwise
     */
    public boolean drawGhostParent() {
        return GHOST_PARENT;
    }

    /**
     * Returns whether hierarchy level numbers are to be omitted when drawing
     * labels.
     * 
     * @return true, omit them; false, include them.
     */
    public boolean removeLvl() {
        return REMOVE_LEVEL_NUMBER;
    }
}
