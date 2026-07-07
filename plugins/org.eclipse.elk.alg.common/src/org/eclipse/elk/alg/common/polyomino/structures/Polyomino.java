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
package org.eclipse.elk.alg.common.polyomino.structures;

import java.util.List;

import org.eclipse.elk.alg.common.utils.UniqueTriple;

import com.google.common.collect.Lists;

/**
 * This class provides a low resolution grid for representing polyominoes.
 */
public class Polyomino extends PlanarGrid {

    ///////////////////////////////////////////////////////////////////////////////
    // Variables

    /**
     * If this polyomino is added to another grid structure, the coordinates of its upper left corner relative to the
     * underlying grid, are set here. This is the x-coordinate.
     */
    private int x = 0;
    /**
     * If this polyomino is added to another grid structure, the coordinates of its upper left corner relative to the
     * underlying grid, are set here. This is the y-coordinate.
     */
    private int y = 0;

    private List<UniqueTriple<Direction, Integer, Integer>> polyominoExtensions;

    ///////////////////////////////////////////////////////////////////////////////
    // Constructors

    /**
     * Creates an empty polyomino of a specified size. Attaches the given extensions to the boundaries of the polyomino.
     * Cells have to be manually set to blocked and weakly blocked afterwards, though.
     * 
     * @param width
     *            width of the polyomino in discrete cells
     * @param height
     *            height of the polyomino in discrete cells
     * @param extensions
     *            extensions defined by a cardinal direction, a horizontal or vertical offset from the origin of the
     *            grid depending on the given direction, and the width of the extension in discrete cells.
     */
    public Polyomino(final int width, final int height,
            final List<UniqueTriple<Direction, Integer, Integer>> extensions) {
        super(width, height);
        polyominoExtensions = extensions;
    }

    /**
     * Creates an empty polyomino of a specified size.
     * 
     * @param width
     *            width of the polyomino in discrete cells
     * @param height
     *            height of the polyomino in discrete cells
     */
    public Polyomino(final int width, final int height) {
        this(width, height, Lists.newArrayList());
    }

    /**
     * Creates a polyomino of zero width and height.
     */
    public Polyomino() {
        this(0, 0);
    }

    ///////////////////////////////////////////////////////////////////////////////
    // Public methods

    /**
     * Gets the x-coordinate of the origin of this polyomino (potentially relative to another grid data structure).
     * 
     * @return X-coordinate
     */
    public int getX() {
        return x;
    }

    /**
     * Gets the y-coordinate of the origin of this polyomino (potentially relative to another grid data structure).
     * 
     * @return Y-coordinate
     */
    public int getY() {
        return y;
    }

    /**
     * Sets the x-coordinate of the origin of this polyomino.
     * 
     * @param val
     *            New value
     */
    public void setX(final int val) {
        x = val;
    }

    /**
     * Sets the y-coordinate of the origin of this polyomino.
     * 
     * @param val
     *            New value
     */
    public void setY(final int val) {
        y = val;
    }

    /**
     * @return the polyominoExtensions
     */
    public List<UniqueTriple<Direction, Integer, Integer>> getPolyominoExtensions() {
        return polyominoExtensions;
    }

    /**
     * Attaches the given extension to the boundaries of the polyomino.
     * 
     * @param dir
     *            cardinal direction of the extension
     * @param offset
     *            a horizontal or vertical offset from the origin of the grid depending on the given direction
     * @param width
     *            the width of the extension in discrete cells
     */
    public void addExtension(final Direction dir, final int offset, final int width) {
        polyominoExtensions.add(new UniqueTriple<Direction, Integer, Integer>(dir, offset, width));
    }

}
