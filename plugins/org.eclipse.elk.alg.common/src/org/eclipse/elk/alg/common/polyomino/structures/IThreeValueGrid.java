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

/**
 * Basic data structure for representing a discrete ℕ x ℕ grid, albeit bounded by a given width and height. Each grid
 * cell can be free or filled.
 */
public interface IThreeValueGrid {

    /**
     * Gets the width of the grid.
     * 
     * @return The width
     */
    int getWidth();

    /**
     * Gets the height of the grid.
     * 
     * @return The height
     */
    int getHeight();
    
    /**
     * Returns whether a given cell is empty or not.
     * 
     * @param x
     *            X-Coordinate of grid cell to be tested
     * @param y
     *            Y-Coordinate of grid cell to be tested
     * @return true: cell is empty; false: otherwise
     * @throws IndexOutOfBoundsException
     *             Coordinates are outside of the grid.
     */
    boolean isEmpty(int x, int y) throws IndexOutOfBoundsException;

    /**
     * Returns whether a given cell is blocked or not.
     * 
     * @param x
     *            X-Coordinate of grid cell to be tested
     * @param y
     *            Y-Coordinate of grid cell to be tested
     * @return true: cell is blocked; false: otherwise
     * @throws IndexOutOfBoundsException
     *             Coordinates are outside of the grid.
     */
    boolean isBlocked(int x, int y) throws IndexOutOfBoundsException;
    
    /**
     * Returns whether a given cell is weakly blocked or not.
     * 
     * @param x
     *            X-Coordinate of grid cell to be tested
     * @param y
     *            Y-Coordinate of grid cell to be tested
     * @return true: cell is weakly blocked; false: otherwise
     * @throws IndexOutOfBoundsException
     *             Coordinates are outside of the grid.
     */
    boolean isWeaklyBlocked(int x, int y) throws IndexOutOfBoundsException;

    /**
     * Returns whether given point (x,y) is inside the grid.
     * 
     * @param x
     *            X-coordinate
     * @param y
     *            Y-coordinate
     * @return true: point is inside the grid; false: otherwise
     */
    boolean inBounds(int x, int y);

    /**
     * reinitializes this grid with new given dimensions and empty cells.
     * 
     * @param width
     *            New width
     * @param height
     *            New Height
     */
    void reinitialize(int width, int height);
    
    /**
     * Empties the cell at position (x,y).
     * 
     * @param x
     *            X-coordinate of cell to empty
     * @param y
     *            Y-coordinate of cell to empty
     * @throws IndexOutOfBoundsException
     *             Coordinates are outside of the grid.
     */
    void setEmpty(int x, int y) throws IndexOutOfBoundsException;

    /**
     * Marks the cell at position (x,y) as blocked.
     * 
     * @param x
     *            X-coordinate of cell to mark as blocked
     * @param y
     *            Y-coordinate of cell to mark as blocked
     * @throws IndexOutOfBoundsException
     *             Coordinates are outside of the grid.
     */
    void setBlocked(int x, int y) throws IndexOutOfBoundsException;
    
    /**
     * Marks the cell at position (x,y) as weakly blocked.
     * 
     * @param x
     *            X-coordinate of cell to mark as weakly blocked
     * @param y
     *            Y-coordinate of cell to mark as weakly blocked
     * @throws IndexOutOfBoundsException
     *             Coordinates are outside of the grid.
     */
    void setWeaklyBlocked(int x, int y) throws IndexOutOfBoundsException;

}
