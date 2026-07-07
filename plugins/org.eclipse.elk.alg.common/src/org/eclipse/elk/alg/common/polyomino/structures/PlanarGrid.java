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

import org.eclipse.elk.alg.common.utils.UniqueTriple;
import org.eclipse.elk.core.util.Quadruple;

/**
 * Binary grid structure with some extra operations, like addressing coordinates of the grid relative to its center,
 * testing for intersections with other planar grids and adding filled cells from other planar grids.
 */
public class PlanarGrid extends TwoBitGrid {

    ///////////////////////////////////////////////////////////////////////////////
    // Variables

    /** Indicates the x-coordinate of the center of the grid (if the center is not integral, value is truncated). */
    private int xCenter;
    /** Indicates the y-coordinate of the center of the grid (if the center is not integral, value is truncated). */
    private int yCenter;

    ///////////////////////////////////////////////////////////////////////////////
    // Constructors

    /**
     * Constructs an empty planar grid of size 0 * 0.
     */
    PlanarGrid() {
        super(0, 0);
        xCenter = 0;
        yCenter = 0;
    }

    /**
     * Constructs a planar grid of size width * height. Computes the center coordinates of this grid. Use positive
     * integer values only.
     * 
     * @param width
     *            The width
     * @param height
     *            The height
     */
    public PlanarGrid(final int width, final int height) {
        super(width, height);
        xCenter = (width - 1) >> 1;
        yCenter = (height - 1) >> 1;

    }

    ///////////////////////////////////////////////////////////////////////////////
    // Public methods

    /**
     * Returns whether a given cell is blocked or not. Coordinates are given relative to the center of the grid.
     * 
     * @param x
     *            X-Coordinate of grid cell to be tested
     * @param y
     *            Y-Coordinate of grid cell to be tested
     * @return true: cell is filled; false: otherwise
     * @throws IndexOutOfBoundsException
     *             Coordinates are outside of the grid.
     */
    public boolean isEmptyCenterBased(final int x, final int y) throws IndexOutOfBoundsException {
        try {
            return isEmpty(x + xCenter, y + yCenter);
        } catch (IndexOutOfBoundsException e) {
            throw new IndexOutOfBoundsException(
                    e.getLocalizedMessage() + " Given center based coordinates were (" + x + ", " + y + ").");
        }
    }

    /**
     * Returns whether a given cell is blocked or not. Coordinates are given relative to the center of the grid.
     * 
     * @param x
     *            X-Coordinate of grid cell to be tested
     * @param y
     *            Y-Coordinate of grid cell to be tested
     * @return true: cell is filled; false: otherwise
     * @throws IndexOutOfBoundsException
     *             Coordinates are outside of the grid.
     */
    public boolean isBlockedCenterBased(final int x, final int y) throws IndexOutOfBoundsException {
        try {
            return isBlocked(x + xCenter, y + yCenter);
        } catch (IndexOutOfBoundsException e) {
            throw new IndexOutOfBoundsException(
                    e.getLocalizedMessage() + " Given center based coordinates were (" + x + ", " + y + ").");
        }
    }

    /**
     * Returns whether a given cell is weakly blocked or not. Coordinates are given relative to the center of the grid.
     * 
     * @param x
     *            X-Coordinate of grid cell to be tested
     * @param y
     *            Y-Coordinate of grid cell to be tested
     * @return true: cell is weakly blocked; false: otherwise
     * @throws IndexOutOfBoundsException
     *             Coordinates are outside of the grid.
     */
    public boolean isWeaklyBlockedCenterBased(final int x, final int y) throws IndexOutOfBoundsException {
        try {
            return isWeaklyBlocked(x + xCenter, y + yCenter);
        } catch (IndexOutOfBoundsException e) {
            throw new IndexOutOfBoundsException(
                    e.getLocalizedMessage() + " Given center based coordinates were (" + x + ", " + y + ").");
        }
    }

    /**
     * Returns whether given point (x,y) is inside the grid. Coordinates are given relative to the center of the grid.
     * 
     * @param x
     *            X-coordinate
     * @param y
     *            Y-coordinate
     * @return true: point is inside the grid; false: otherwise
     */
    public boolean inBoundsCenterBased(final int x, final int y) {
        int xt = x;
        int yt = y;
        xt += xCenter;
        yt += yCenter;
        return inBounds(xt, yt);
    }

    /**
     * Tests whether the filled cells of two grids intersect with which each other, whereby the center of both grids are
     * placed upon each other (unless there is an offset). Weakly blocked cells are allowed to overlap each other
     * without being determined intersected, but blocked cells instersect with weakly blocked cells.
     * 
     * @param other
     *            Other grid to be tested for intersection with this grid
     * @param xOffset
     *            Offset of the center of other relative to the center of this grid in x-direction
     * @param yOffset
     *            Offset of the center of other relative to the center of this grid in y-direction
     * @return true: both grids overlap/intersect with each other; false: otherwise
     */
    public <G extends PlanarGrid> boolean intersectsWithCenterBased(final G other, final int xOffset,
            final int yOffset) {
        for (int x = 0; x < other.getWidth(); x++) {
            int xTranslated = x - other.getCenterX() + xOffset;
            for (int y = 0; y < other.getHeight(); y++) {
                int yTranslated = y - other.getCenterY() + yOffset;
                if (this.inBoundsCenterBased(xTranslated, yTranslated)
                        && ((!(other.isEmpty(x, y)) && this.isBlockedCenterBased(xTranslated, yTranslated))
                                || (other.isBlocked(x, y) && !(this.isEmptyCenterBased(xTranslated, yTranslated))))) {
                    return true;
                }
            }
        }
        return false;
    }

    /**
     * Tests whether the filled cells of two grids intersect with which each other, whereby the center of both grids are
     * placed upon each other (unless there is an offset). Weakly blocked cells are allowed to overlap each other
     * without being determined intersected, but blocked cells intersect with weakly blocked cells.
     * 
     * @param other
     *            Other grid to be tested for intersection with this grid
     * @param xOffset
     *            Offset of the center of other relative to the center of this grid in x-direction
     * @param yOffset
     *            Offset of the center of other relative to the center of this grid in y-direction
     * @return true: both grids overlap/intersect with each other; false: otherwise
     */
    public <G extends PlanarGrid> boolean intersectsWithCenterBased(final Polyomino other, final int xOffset,
            final int yOffset) {
        if (this.intersectsWithCenterBased((PlanarGrid) other, xOffset, yOffset)) {
            return true;
        } else {
            for (UniqueTriple<Direction, Integer, Integer> ext : other.getPolyominoExtensions()) {
                boolean intersects = false;

                // Some transformations of the coordinates, as center based coordinates aren't that practical for
                // checking the extensions
                int leftX = getCenterX() - other.getCenterX() + xOffset;
                int rightX = leftX + other.getWidth();
                int topY = getCenterY() - other.getCenterY() + yOffset;
                int bottomY = topY + other.getHeight();

                switch (ext.getFirst()) {
                case NORTH:
                    intersects = weaklyIntersectsArea(leftX + ext.getSecond(), 0, leftX + ext.getThird(), topY - 1);
                    break;
                case EAST:
                    intersects =
                            weaklyIntersectsArea(rightX, topY + ext.getSecond(), getWidth() - 1, topY + ext.getThird());
                    break;
                case SOUTH:
                    intersects = weaklyIntersectsArea(leftX + ext.getSecond(), bottomY, leftX + ext.getThird(),
                            getHeight() - 1);
                    break;
                default: // WEST
                    intersects = weaklyIntersectsArea(0, topY + ext.getSecond(), leftX - 1, topY + ext.getThird());
                }
                if (intersects) {
                    return true;
                }
            }
        }
        return false;
    }

    /**
     * Adds blocked and weakly blocked cells from another grid to this one, whereby the center of both grids are placed
     * upon each other (unless there is an offset). Existing blocked cells won't be overwritten.
     * 
     * @param other
     *            Grid to add cells from
     * @param xOffset
     *            Offset of the center of other relative to the center of this grid in x-direction
     * @param yOffset
     *            Offset of the center of other relative to the center of this grid in y-direction
     * @throws IndexOutOfBoundsException
     *             not all cells of other fit into this grid (considering the given offset).
     */
    public <G extends PlanarGrid> void addFilledCellsFrom(final G other, final int xOffset, final int yOffset)
            throws IndexOutOfBoundsException {
        for (int x = 0; x < other.getWidth(); x++) {
            int xTranslated = x - other.getCenterX() + xOffset;
            for (int y = 0; y < other.getHeight(); y++) {
                int yTranslated = y - other.getCenterY() + yOffset;
                if (other.isBlocked(x, y)) {
                    if (!this.isWeaklyBlockedCenterBased(xTranslated, yTranslated)) {
                        // Might throw exception.
                        this.setBlockedCenterBased(xTranslated, yTranslated);
                    }
                } else {
                    // (priorities for overwriting: blocked > weakly_blocked > empty)
                    // Might throw exception.
                    if (other.isWeaklyBlocked(x, y)) {
                        if (!this.isBlockedCenterBased(xTranslated, yTranslated)) {
                            this.setWeaklyBlockedCenterBased(xTranslated, yTranslated);
                        }
                    }
                }
            }
        }
    }

    /**
     * Adds filled cells from a {@link Polyomino} to this one, whereby the center of the {@link Polyomino} is placed on
     * (unless there is an offset). Additionally sets fields in the {@link Polyomino} containing its upper left corner
     * coordinates relative to the upper left corner of this grid, considering the offset.
     * 
     * @param other
     *            Grid to add cells from
     * @param xOffset
     *            Offset of the center of other relative to the center of this grid in x-direction
     * @param yOffset
     *            Offset of the center of other relative to the center of this grid in y-direction
     * @throws IndexOutOfBoundsException
     *             not all cells of other fit into this grid (considering the given offset).
     */
    public void addFilledCellsFrom(final Polyomino other, final int xOffset, final int yOffset)
            throws IndexOutOfBoundsException {
        // Might throw exception.
        addFilledCellsFrom((PlanarGrid) other, xOffset, yOffset);
        other.setX(this.xCenter - other.getCenterX() + xOffset);
        other.setY(this.yCenter - other.getCenterY() + yOffset);

        for (UniqueTriple<Direction, Integer, Integer> ext : other.getPolyominoExtensions()) {
            switch (ext.getFirst()) {
            case NORTH:
                weaklyBlockArea(other.getX() + ext.getSecond(), 0, other.getX() + ext.getThird(), other.getY() - 1);
                break;
            case EAST:
                weaklyBlockArea(other.getX() + other.getWidth(), other.getY() + ext.getSecond(), getWidth() - 1,
                        other.getY() + ext.getThird());
                break;
            case SOUTH:
                weaklyBlockArea(other.getX() + ext.getSecond(), other.getY() + other.getHeight(),
                        other.getX() + ext.getThird(), getHeight() - 1);
                break;
            default: // WEST
                weaklyBlockArea(0, other.getY() + ext.getSecond(), other.getX() - 1, other.getY() + ext.getThird());
            }
        }
    }

    /**
     * Gets the bounds of the rectangle bounding the filled cells of this grid.
     * 
     * @return (x, y, width, height) of the bounding rectangle. (x, y) is the position of the upper left corner of the
     *         rectangle raltive to the upper left corner of this grid.
     */
    public Quadruple<Integer, Integer, Integer, Integer> getFilledBounds() {
        int gridWidth = getWidth(), gridHeight = getHeight();
        int minX = Integer.MAX_VALUE;
        int maxX = Integer.MIN_VALUE;
        int minY = Integer.MAX_VALUE;
        int maxY = Integer.MIN_VALUE;
        for (int xi = 0; xi < gridWidth; ++xi) {
            for (int yi = 0; yi < gridHeight; ++yi) {
                if (isBlocked(xi, yi)) {
                    minX = Math.min(minX, xi);
                    maxX = Math.max(maxX, xi);
                    minY = Math.min(minY, yi);
                    maxY = Math.max(maxY, yi);
                }
            }
        }
        int width = maxX - minX + 1;
        int height = maxY - minY + 1;

        return new Quadruple<Integer, Integer, Integer, Integer>(minX, minY, width, height);
    }

    /**
     * Weakly blocks all cells in a given rectangular area given in coordinates with the upper left corner as origin,
     * doesn't overwrite blocked cells.
     * 
     * @param xUpperLeft
     *            x-coordinate of upper left corner
     * @param yUpperleft
     *            y-coordinate of upper left corner
     * @param xBottomRight
     *            x-coordinate of bottom right corner
     * @param yBottomRight
     *            y-coordinate of bottom right corner
     * @throws IndexOutOfBoundsException
     *             One of given corner points is out of bounds.
     */
    public void weaklyBlockArea(final int xUpperLeft, final int yUpperleft, final int xBottomRight,
            final int yBottomRight) {
        for (int yi = yUpperleft; yi <= yBottomRight; yi++) {
            for (int xi = xUpperLeft; xi <= xBottomRight; xi++) {
                if (!isBlocked(xi, yi)) {
                    // Might throw IndexOutOfBoundsException
                    setWeaklyBlocked(xi, yi);
                }
            }
        }
    }

    /**
     * Returns whether a given rectangular area given in coordinates with the upper left corner as origin contains at
     * least one blocked cell.
     * 
     * @param xUpperLeft
     *            x-coordinate of upper left corner
     * @param yUpperleft
     *            y-coordinate of upper left corner
     * @param xBottomRight
     *            x-coordinate of bottom right corner
     * @param yBottomRight
     *            y-coordinate of bottom right corner
     * @return true, area contains blocked cell(s), false, otherwise
     * @throws IndexOutOfBoundsException
     *             One of given corner points is out of bounds.
     */
    public boolean weaklyIntersectsArea(final int xUpperLeft, final int yUpperleft, final int xBottomRight,
            final int yBottomRight) {
        for (int yi = yUpperleft; yi <= yBottomRight; yi++) {
            for (int xi = xUpperLeft; xi <= xBottomRight; xi++) {
                // Might throw IndexOutOfBoundsException
                if (isBlocked(xi, yi)) {
                    return true;
                }
            }
        }
        return false;

    }

    /**
     * Gets the x-coordinate of the (truncated) center of this grid.
     * 
     * @return X-coordinate of the center of the grid
     */
    public int getCenterX() {
        return xCenter;
    }

    /**
     * Gets the y-coordinate of the (truncated) center of this grid.
     * 
     * @return Y-coordinate of the center of the grid
     */
    public int getCenterY() {
        return yCenter;
    }

    ///////////////////////////////////////////////////////////////////////////////
    // Package private methods

    @Override
    public void reinitialize(final int width, final int height) {
        super.reinitialize(width, height);
        xCenter = (width - 1) >> 1;
        yCenter = (height - 1) >> 1;
    }

    /**
     * Blocks the cell at position (x,y). Coordinates are given relative to the center of the grid.
     * 
     * @param x
     *            X-coordinate of cell to fill
     * @param y
     *            Y-coordinate of cell to fill
     * @throws IndexOutOfBoundsException
     *             Coordinates are outside of the grid.
     */
    void setBlockedCenterBased(final int x, final int y) throws IndexOutOfBoundsException {
        try {
            setBlocked(x + xCenter, y + yCenter);
        } catch (IndexOutOfBoundsException e) {
            throw new IndexOutOfBoundsException(
                    e.getLocalizedMessage() + " Given center based coordinates were (" + x + ", " + y + ").");
        }
    }

    /**
     * Empties the cell at position (x,y). Coordinates are given relative to the center of the grid.
     * 
     * @param x
     *            X-coordinate of cell to empty
     * @param y
     *            Y-coordinate of cell to empty
     * @throws IndexOutOfBoundsException
     *             Coordinates are outside of the grid.
     */
    void setEmptyCenterBased(final int x, final int y) throws IndexOutOfBoundsException {
        try {
            setEmpty(x + xCenter, y + yCenter);
        } catch (IndexOutOfBoundsException e) {
            throw new IndexOutOfBoundsException(
                    e.getLocalizedMessage() + " Given center based coordinates were (" + x + ", " + y + ").");
        }
    }

    /**
     * Fills the cell at position (x,y). Coordinates are given relative to the center of the grid.
     * 
     * @param x
     *            X-coordinate of cell to weakly block
     * @param y
     *            Y-coordinate of cell to weakly block
     * @throws IndexOutOfBoundsException
     *             Coordinates are outside of the grid.
     */
    void setWeaklyBlockedCenterBased(final int x, final int y) throws IndexOutOfBoundsException {
        try {
            setWeaklyBlocked(x + xCenter, y + yCenter);
        } catch (IndexOutOfBoundsException e) {
            throw new IndexOutOfBoundsException(
                    e.getLocalizedMessage() + " Given center based coordinates were (" + x + ", " + y + ").");
        }
    }

}
