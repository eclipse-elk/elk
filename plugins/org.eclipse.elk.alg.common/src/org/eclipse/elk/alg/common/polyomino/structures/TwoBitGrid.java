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
 * Basic data structure for representing a discrete ℕ x ℕ grid, albeit bounded by a given width and height. The data
 * structure is based on a two dimensional long array, making it possible to address each bit individually to save
 * space.
 */
public class TwoBitGrid implements IThreeValueGrid {
    ///////////////////////////////////////////////////////////////////////////////
    // Constants

    private static final long LSB_MASK = 0x01;
    private static final long TWO_LSBS_MASK = 0x03;
    private static final long EMPTY = 0x00;
    private static final long BLOCKED = 0x01;
    private static final long WEAKLY_BLOCKED = 0x02;
    private static final double HALF_WORD = 32.0;
    private static final int REST_MASK = 0x1F;
    private static final int RIGHT_SHIFT = 5; // 2^5 = 32

    ///////////////////////////////////////////////////////////////////////////////
    // Variables

    /** The implementation uses an encapsulated byte array. */
    private long[][] grid;
    /** The dimensions of the current grid. */
    private int xSize, ySize;

    ///////////////////////////////////////////////////////////////////////////////
    // Package private Constructors

    /**
     * Constructs an empty grid of size 0 * 0.
     */
    TwoBitGrid() {
        this(0, 0);
    }

    /**
     * Constructs a grid of size width * height.
     * 
     * @param width
     *            The width
     * @param height
     *            The height
     */
    public TwoBitGrid(final int width, final int height) {
        grid = new long[height][(int) Math.ceil(width / HALF_WORD)];
        xSize = width;
        ySize = height;
    }

    ///////////////////////////////////////////////////////////////////////////////
    // Public methods

    @Override
    public int getWidth() {
        return xSize;
    }

    @Override
    public int getHeight() {
        return ySize;
    }

    @Override
    public boolean isEmpty(final int x, final int y) throws IndexOutOfBoundsException {
        try {
            return retrieve(x, y) == EMPTY ? true : false;
        } catch (ArrayIndexOutOfBoundsException e) {
            throw new IndexOutOfBoundsException("Grid is only of size " + xSize + "*" + ySize + ". Requested point ("
                    + x + ", " + y + ") is out of bounds.");
        }
    }

    @Override
    public boolean isBlocked(final int x, final int y) throws IndexOutOfBoundsException {
        try {
            return retrieve(x, y) == BLOCKED ? true : false;
        } catch (ArrayIndexOutOfBoundsException e) {
            throw new IndexOutOfBoundsException("Grid is only of size " + xSize + "*" + ySize + ". Requested point ("
                    + x + ", " + y + ") is out of bounds.");
        }
    }

    @Override
    public boolean isWeaklyBlocked(final int x, final int y) throws IndexOutOfBoundsException {
        try {
            return retrieve(x, y) == WEAKLY_BLOCKED ? true : false;
        } catch (ArrayIndexOutOfBoundsException e) {
            throw new IndexOutOfBoundsException("Grid is only of size " + xSize + "*" + ySize + ". Requested point ("
                    + x + ", " + y + ") is out of bounds.");
        }
    }

    @Override
    public boolean inBounds(final int x, final int y) {
        return x >= 0 && y >= 0 && x < xSize && y < ySize;
    }

    @Override
    public void reinitialize(final int width, final int height) {
        grid = new long[height][(int) Math.ceil(width / HALF_WORD)];
        xSize = width;
        ySize = height;

    }

    @Override
    public void setEmpty(final int x, final int y) throws IndexOutOfBoundsException {
        set(x, y, false, false);

    }

    @Override
    public void setBlocked(final int x, final int y) throws IndexOutOfBoundsException {
        set(x, y, false, true);
    }

    @Override
    public void setWeaklyBlocked(final int x, final int y) throws IndexOutOfBoundsException {
        set(x, y, true, false);
    }

    @Override
    public String toString() {
        String output = " ";
        // A basic coordinate system
        Integer count = 0;
        for (int x = 0; x < xSize; x++) {
            output += count.toString();
            count = incModTen(count);
        }
        output += "\n";
        count = 0;
        for (int y = 0; y < ySize; y++) {
            output += count.toString();
            count = incModTen(count);
            for (int x = 0; x < xSize; x++) {
                long item = retrieve(x, y);
                if (item == EMPTY) {
                    output += "_";
                } else if (item == BLOCKED) {
                    output += "X";
                } else {
                    output += "0";
                }
            }
            output += "\n";
        }
        return output.substring(0, output.length() - 1);
    }

    ///////////////////////////////////////////////////////////////////////////////
    // Private methods

    /**
     * Returns the current state of a given cell.
     * 
     * @param x
     *            X-Coordinate of grid cell to be tested
     * @param y
     *            Y-Coordinate of grid cell to be tested
     * @return {code BLOCKED}: cell is blocked; {code WEAKLY_BLOCKED}: cell is weakly blocked; {@code EMPTY}: otherwise
     * @throws ArrayIndexOutOfBoundsException
     *             Coordinates are outside of the grid.
     */
    private long retrieve(final int x, final int y) throws ArrayIndexOutOfBoundsException {
        int xWord = x >> RIGHT_SHIFT;
        long xRest = x & REST_MASK;
        long value = (grid[y][xWord] >>> (xRest << 1)) & TWO_LSBS_MASK;
        return value;
    }

    /**
     * Sets the two bits of the cell at position (x,y).
     * 
     * @param x
     *            X-coordinate of cell
     * @param y
     *            Y-coordinate of cell
     * @param msb
     *            true: most significant bit of the cell will be set to 1, otherwise: 0
     * @param lsb
     *            true: least significant bit of the cell will be set to 1, otherwise: 0
     * @throws IndexOutOfBoundsException
     *             Coordinates are outside of the grid.
     */
    private void set(final int x, final int y, final boolean msb, final boolean lsb) throws IndexOutOfBoundsException {
        try {
            if (x >= xSize) {
                throw new ArrayIndexOutOfBoundsException();
            }

            int xWord = x >> RIGHT_SHIFT;
            long xRest = x & REST_MASK;
            long mask = LSB_MASK << (xRest << 1);
            if (lsb) {
                grid[y][xWord] = grid[y][xWord] | mask;
            } else {
                grid[y][xWord] = grid[y][xWord] & ~mask;
            }
            mask <<= 1;
            if (msb) {
                grid[y][xWord] = grid[y][xWord] | mask;
            } else {
                grid[y][xWord] = grid[y][xWord] & ~mask;
            }

        } catch (ArrayIndexOutOfBoundsException e) {
            throw new IndexOutOfBoundsException("Grid is only of size " + xSize + "*" + ySize + ". Requested point ("
                    + x + ", " + y + ") is out of bounds.");
        }
    }

    /**
     * Increments a number by 1 if it is between 0 and 8, returns 0 otherwise.
     * 
     * @param num
     *            number to increment
     * @return number incremented by 1 if it is between 0 and 8, 0 otherwise
     */
    private int incModTen(final int num) {
        final int eight = 8;
        if (num > eight) {
            return 0;
        }
        return num + 1;
    }

}
