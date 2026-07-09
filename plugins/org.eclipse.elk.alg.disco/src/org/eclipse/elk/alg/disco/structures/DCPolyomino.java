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
package org.eclipse.elk.alg.disco.structures;

import java.util.List;

import org.eclipse.elk.alg.common.polyomino.structures.Direction;
import org.eclipse.elk.alg.common.polyomino.structures.Polyomino;
import org.eclipse.elk.alg.disco.graph.DCComponent;
import org.eclipse.elk.alg.disco.graph.DCDirection;
import org.eclipse.elk.alg.disco.graph.DCElement;
import org.eclipse.elk.alg.disco.graph.DCExtension;
import org.eclipse.elk.core.math.ElkRectangle;
import org.eclipse.elk.core.math.KVector;

/**
 * This class provides a low resolution binary grid for representing {@link DCComponent DCComponents}.
 */
public class DCPolyomino extends Polyomino {

    ///////////////////////////////////////////////////////////////////////////////
    // Variables

    /** {@link DCComponent} represented by this polyomino. */
    private DCComponent representee;
    /** The width of this polyomino. */
    private int pWidth;
    /** The height of this polyomino. */
    private int pHeight;
    /** Width of one polyomino cell in the {@link DCComponent DCComponent's} coordinate system. */
    private double cellSizeX;
    /** Height of one polyomino cell in the {@link DCComponent DCComponent's} coordinate system. */
    private double cellSizeY;

    ///////////////////////////////////////////////////////////////////////////////
    // Constructor

    /**
     * Constructs a Polyomino out of a given {@link DCComponent}.
     * 
     * @param comp
     *            {@link DCComponent} to be transformed
     * @param csX
     *            Step size for rasterizing the {@link DCComponent}
     * @param csY
     *            Step size for rasterizing the {@link DCComponent}
     */
    public DCPolyomino(final DCComponent comp, final double csX, final double csY) {
        super();
        cellSizeX = csX;
        cellSizeY = csY;
        representee = comp;
        KVector compDims = comp.getDimensionsOfBoundingRectangle();

        pWidth = computeLowResDimension(compDims.x, cellSizeX);
        pHeight = computeLowResDimension(compDims.y, cellSizeY);

        reinitialize(pWidth, pHeight);

        fillCells();

        for (DCElement elem : representee.getElements()) {
            if (elem.getExtensions().size() > 0) {
                addExtensionsToPoly(elem);
            }
        }
    }

    ///////////////////////////////////////////////////////////////////////////////
    // Public methods

    /**
     * Gets the offset of the upper left corner this polyomino relative to the upper left corner of the
     * {@link DCComponent} it represents.
     * 
     * @return The offset
     */
    public KVector getOffset() {
        final double halfNeg = -0.5d;
        return representee.getDimensionsOfBoundingRectangle().clone().sub(pWidth * cellSizeX, pHeight * cellSizeY)
                .scale(halfNeg);
    }

    // The following public methods are intended for DEBUGGING PURPOSES (s. DisCoDebugView).

    /**
     * Computes the supposed position of the upper left corner of the polyomino, if it where drawn in the same
     * coordinate system as the {@link DCComponent} it represents.
     * 
     * @return Position the polyomino would have, if it were drawn to scale on the same canvas as its associated
     *         {@link DCComponent}
     */
    public KVector getMinCornerOnCanvas() {
        return representee.getMinCorner().clone().sub(getOffset());
    }

    /**
     * Gets the length of a side of a grid cell of the polyomino if translated to the original coordinate system of its
     * respective {@link DCComponent}.
     * 
     * @return the cellSize
     */
    public double getCellSizeX() {
        return cellSizeX;
    }

    /**
     * Gets the length of a side of a grid cell of the polyomino if translated to the original coordinate system of its
     * respective {@link DCComponent}.
     * 
     * @return the cellSize
     */
    public double getCellSizeY() {
        return cellSizeY;
    }

    /**
     * Gets the id of this polyomino.
     * 
     * @return The id
     */
    public int getId() {
        return representee.getId();
    }

    /**
     * Sets the id of this polyomino.
     * 
     * @param id
     *            The id
     */
    public void setId(final int id) {
        representee.setId(id);
    }

    /**
     * Gets the {@link DCComponent} associated with this polyomino.
     * 
     * @return The associated component
     */
    public DCComponent getRepresentee() {
        return representee;
    }

    ///////////////////////////////////////////////////////////////////////////////
    // Private methods

    /**
     * Computes one dimension of the polyomino if given the width or height of its associated component.
     * 
     * @param dim
     *            Width oder height of this polyomino {@link DCComponent}.
     * @return Integral low resolution width or height of this polyomino
     */
    private int computeLowResDimension(final double dim, final double cellSize) {
        double cellFit = dim / cellSize;
        int fitTruncated = (int) cellFit;
        if (cellFit > fitTruncated) {
            fitTruncated++;
        }
        return fitTruncated;
    }

    /**
     * Fills the (assumedly) empty Polyomino grid cells, if they intersect with the associated {@link DCComponent}.
     */
    private void fillCells() {

        KVector compCorner = representee.getMinCorner();
        KVector polyoOffset = getOffset();

        double baseX = compCorner.x - polyoOffset.x;
        double curX;
        double curY = compCorner.y - polyoOffset.y;

        // TODO Option for filling enclosed spaces, rectilinear convex/concave hull

        for (int y = 0; y < pHeight; y++) {
            curX = baseX;
            for (int x = 0; x < pWidth; x++) {
                if (representee.intersects(new ElkRectangle(curX, curY, cellSizeX, cellSizeY))) {
                    setBlocked(x, y);
                }
                curX += cellSizeX;
            }
            curY += cellSizeY;
        }
    }

    /**
     * Adds extenstions to this polyomino taken from a {@link DCElement} of the component it represents.
     * 
     * @param elem
     *            an element of the component this polyomino represents
     */
    private void addExtensionsToPoly(final DCElement elem) {
        List<DCExtension> extensions = elem.getExtensions();

        KVector compCorner = representee.getMinCorner();
        KVector polyoOffset = getOffset();

        double baseX = compCorner.x - polyoOffset.x;
        double baseY = compCorner.y - polyoOffset.y;

        ElkRectangle elemPos = elem.getBounds();

        baseX = elemPos.x - baseX;
        baseY = elemPos.y - baseY;

        for (DCExtension extension : extensions) {
            KVector pos = extension.getOffset();
            double xe = baseX + pos.x;
            double ye = baseY + pos.y;

            int xp = (int) (xe / cellSizeX);
            int yp = (int) (ye / cellSizeY);

            DCDirection dir = extension.getDirection();
            Direction polyDir;
            switch (dir) {
            case NORTH:
                polyDir = Direction.NORTH;
                break;
            case EAST:
                polyDir = Direction.EAST;
                break;
            case SOUTH:
                polyDir = Direction.SOUTH;
                break;
            default:
                polyDir = Direction.WEST;
            }

            if (dir.isHorizontal()) {
                int ypPlusWidth = (int) ((ye + extension.getWidth()) / cellSizeY);
                addExtension(polyDir, yp, ypPlusWidth);
                if (dir.equals(DCDirection.WEST)) {
                    weaklyBlockArea(0, yp, xp, ypPlusWidth);
                } else { // direction is EAST
                    weaklyBlockArea(xp, yp, pWidth - 1, ypPlusWidth);
                }

            } else { // direction is vertical
                int xpPlusWidth = (int) ((xe + extension.getWidth()) / cellSizeX);
                addExtension(polyDir, xp, xpPlusWidth);
                if (dir.equals(DCDirection.NORTH)) {
                    weaklyBlockArea(xp, 0, xpPlusWidth, yp);
                } else { // direction is SOUTH
                    weaklyBlockArea(xp, yp, xpPlusWidth, pHeight - 1);
                }
            }
        }
    }

}
