/*******************************************************************************
 * Copyright (c) 2018, 2020 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.rectpacking.p1widthapproximation;

import java.util.List;

import org.eclipse.elk.alg.rectpacking.util.DrawingData;
import org.eclipse.elk.alg.rectpacking.util.DrawingDataDescriptor;
import org.eclipse.elk.graph.ElkNode;

/**
 * Class that offers methods for small operations in the first iteration.
 * 
 * @see AreaApproximation
 */
public final class Calculations {

    private Calculations() {
    }

    /**
     * Calculates the drawing width for placement options LPR and LPB.
     * 
     * @param drawingWidth The width of drawing before placement of new rectangle.
     * @param x The potential x-coordinate of new rectangle.
     * @param width The width of new rectangle.
     * @return The width of drawing including new rectangle in case LPR and LPB.
     */
    protected static double getWidthLPRorLPB(final double drawingWidth, final double x,
            final double width) {
        return Math.max(drawingWidth, x + width);
    }

    /**
     * Calculates the drawing height for placement options LPR and LPB.
     * 
     * @param drawingHeight The height of drawing before placement of new rectangle.
     * @param y The potential y-coordinate of new rectangle.
     * @param height The height of new rectangle.
     * @return The height of drawing including new rectangle in case LPR and LPB.
     */
    protected static double getHeightLPRorLPB(final double drawingHeight, final double y,
            final double height) {
        return Math.max(drawingHeight, y + height);
    }

    /**
     * Calculates the y-coordinate after the shift of the rectangle to be placed right of lastPlaced.
     * 
     * @param x The x-coordinate of the rectangle to be placed.
     * @param placedRects A list of already placed rectangles.
     * @param lastPlaced The rectangle that was placed most recently.
     * @return The y-coordinate after the shift of the rectangle to be placed right of lastPlaced
     */
    protected static double calculateYforLPR(final double x, final List<ElkNode> placedRects,
            final ElkNode lastPlaced, final double nodeNodeSpacing) {
        ElkNode closestUpperNeighbor = null;
        double closestNeighborBottomBorder = 0;
        // find neighbors that lay between the upper and lower border of the rectangle to be placed.
        for (ElkNode placedRect : placedRects) {
            double placedRectBottomBorder = placedRect.getY() + placedRect.getHeight();
            if (verticalOrderConstraint(placedRect, x, nodeNodeSpacing)) {
                // is closest neighbor?
                if (closestUpperNeighbor == null) {
                    closestUpperNeighbor = placedRect;
                } else if (lastPlaced.getY() - placedRectBottomBorder < lastPlaced.getY() - closestNeighborBottomBorder) {
                    closestUpperNeighbor = placedRect;
                }
                closestNeighborBottomBorder = closestUpperNeighbor.getY() + closestUpperNeighbor.getHeight();
            }
        }

        // no neighbor yet
        if (closestUpperNeighbor == null) {
            return 0;
        } else {
            // else, choose closest neighbors bottom border
            return closestNeighborBottomBorder + nodeNodeSpacing;
        }
    }

    /**
     * Calculates the x-coordinate after shift of the rectangle to be placed below the last placed rectangle.
     * 
     * @param placedRects A list of already placed rectangles.
     * @param y The y-coordinate of the rectangle to be placed.
     * @param lastPlaced The most recently placed rectangle.
     * @return The x-coordinate after shift of the rectangle to be placed below lastPlaced.
     */
    protected static double calculateXforLPB(final double y, final List<ElkNode> placedRects,
            final ElkNode lastPlaced, final double nodeNodeSpacing) {
        ElkNode closestLeftNeighbour = null;
        double closestNeighborRightBorder = 0;
        // Find neighbors that lay in between the height of the rectangle to be placed.
        for (ElkNode placedRect : placedRects) {
            double placedRectRightBorder = placedRect.getX() + placedRect.getWidth();
            if (horizontalOrderConstraint(placedRect, y, nodeNodeSpacing)) {
                // Is closest neighbor?
                if (closestLeftNeighbour == null) {
                    closestLeftNeighbour = placedRect;
                } else if (lastPlaced.getX() - placedRectRightBorder < lastPlaced.getX() - closestNeighborRightBorder) {
                    closestLeftNeighbour = placedRect;
                }
                closestNeighborRightBorder = closestLeftNeighbour.getX() + closestLeftNeighbour.getWidth();
            }
        }

        // No neighbor yet.
        if (closestLeftNeighbour == null) {
            return 0;
        } else {
            // Else, choose closest neighbors right border
            return closestNeighborRightBorder + nodeNodeSpacing;
        }
    }

    /**
     * Calculates the area taken up by the last placed rectangle and the rectangle to place for the
     * {@link DrawingDataDescriptor} last placed right.
     * 
     * @param lastPlaced The last placed rectangle.
     * @param toPlace The rectangle to be placed.
     * @param lprOpt The drawing data containing the coordinates of the rectangle to be placed.
     * @return The area taken up by the two given rectangles.
     */
    protected static double calculateAreaLPR(final ElkNode lastPlaced, final ElkNode toPlace,
            final DrawingData lprOpt) {
        double lastPlacedBottomBorder = lastPlaced.getY() + lastPlaced.getHeight();
        double toPlaceBottomBorder = lprOpt.getNextYcoordinate() + toPlace.getHeight();
        double maxYLPR = Math.max(lastPlacedBottomBorder, toPlaceBottomBorder);

        double heightLPR = maxYLPR - Math.min(lastPlaced.getY(), lprOpt.getNextYcoordinate());
        double widthLPR = lprOpt.getNextXcoordinate() + toPlace.getWidth() - lastPlaced.getX();

        return widthLPR * heightLPR;
    }

    /**
     * Calculates the area taken up by the last placed rectangle and the rectangle to place for the
     * {@link DrawingDataDescriptor} last placed below.
     * 
     * @param lastPlaced The last placed rectangle.
     * @param toPlace The rectangle to be placed.
     * @param lpbOpt The drawing data containing the coordinates of the rectangle to be placed.
     * @return The area taken up by the two given rectangles.
     */
    protected static double calculateAreaLPB(final ElkNode lastPlaced, final ElkNode toPlace,
            final DrawingData lpbOpt) {
        double lastPlacedRightBorder = lastPlaced.getX() + lastPlaced.getWidth();
        double toPlaceRightBorder = lpbOpt.getNextXcoordinate() + toPlace.getWidth();
        double maxXLPB = Math.max(lastPlacedRightBorder, toPlaceRightBorder);

        double widthLPB = maxXLPB - Math.min(lastPlaced.getX(), lpbOpt.getNextXcoordinate());
        double heightLPB = lpbOpt.getNextYcoordinate() + toPlace.getHeight() - lastPlaced.getY();

        return widthLPB * heightLPB;
    }

    /**
     * Checks whether the placedRect produces an vertical order constraint regarding the order. If toPlace is placed
     * left of placedRect, toPlace can at most be placed at the bottom border of placedRect.
     * 
     * @param placedRect The already placed rectangle.
     * @param x The x-coordinate of the rectangle to be placed.
     * @return True, if the placedRect produces a constraint. False otherwise.
     */
    private static boolean verticalOrderConstraint(final ElkNode placedRect, final double x, final double nodeNodeSpacing) {
        return x < placedRect.getX() + placedRect.getWidth() + nodeNodeSpacing;
    }

    /**
     * Checks whether the placedRect produces an horizontal order constraint regarding the order. If toPlace is placed
     * above of placedRect, toPlace can at most be placed at the right border of placedRect.
     * 
     * @param placedRect The already placed rectangle.
     * @param yCoordRectToPlace The y-coordinate of the rectangle to be placed.
     * @return True, if the placedRect produces a constraint. False otherwise.
     */
    private static boolean horizontalOrderConstraint(final ElkNode placedRect, final double yCoordRectToPlace, final double nodeNodeSpacing) {
        return yCoordRectToPlace < placedRect.getY() + placedRect.getHeight() + nodeNodeSpacing;
    }
}
