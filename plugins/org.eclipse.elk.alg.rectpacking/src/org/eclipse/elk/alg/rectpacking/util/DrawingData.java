/*******************************************************************************
 * Copyright (c) 2018, 2020 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.alg.rectpacking.util;

/**
 * Class that offers instances to save information about a possible drawing. Can be used to return multiple values that
 * are calculated for a possible drawing in a method for example.
 */
public class DrawingData {
    //////////////////////////////////////////////////////////////////
    // FIELDS
    /** Scale measure of the given drawing. */
    private double scaleMeasure;
    /** Width of the drawing. */
    private double drawingWidth;
    /** Height of the drawing. */
    private double drawingHeight;
    /** Drawing area. */
    private double area;
    /** Aspect ratio of the drawing. */
    private double aspectRatio;
    /** Desired aspect ratio. */
    private double dar;
    /** Indicates what placement option this data belongs to or whether it belongs to a whole drawing. */
    private DrawingDataDescriptor placementOption;
    /**
     * If a this object contains information about a possible drawing, this field contains the potential x-coordinate
     * for the new rectangle.
     */
    private double nextXcoordinate;
    /**
     * If a this object contains information about a possible drawing, this field contains the potential y-coordinate
     * for the new rectangle.
     */
    private double nextYcoordinate;

    //////////////////////////////////////////////////////////////////
    // Constructors.
    /**
     * Creates a new object saving parameters of a drawing.
     * 
     * @param dar
     *            desired aspect ratio.
     * @param drawingWidth
     *            drawing width.
     * @param drawingHeight
     *            drawing height.
     * @param placementOption
     *            placement option.
     */
    public DrawingData(final double dar, final double drawingWidth, final double drawingHeight,
            final DrawingDataDescriptor placementOption) {
        this(dar, drawingWidth, drawingHeight, placementOption, 0.0, 0.0);
    }

    /**
     * Creates a new object saving parameters of a drawing.
     * 
     * @param dar
     *            desired aspect ratio.
     * @param drawingWidth
     *            drawing width.
     * @param drawingHeight
     *            drawing height.
     * @param placementOption
     *            placement option.
     * @param nextXcoord
     *            x-coordinate for rectangle to place.
     * @param nextYcoord
     *            y-coordinate for rectangle to place.
     */
    public DrawingData(final double dar, final double drawingWidth, final double drawingHeight,
            final DrawingDataDescriptor placementOption, final double nextXcoord, final double nextYcoord) {
        this.dar = dar;
        this.drawingWidth = drawingWidth;
        this.drawingHeight = drawingHeight;
        this.placementOption = placementOption;
        this.nextXcoordinate = nextXcoord;
        this.nextYcoordinate = nextYcoord;
        calcAreaAspectRatioScaleMeasure();
    }

    //////////////////////////////////////////////////////////////////
    // Helper methods.

    /**
     * Calculates the area, aspect ratio, and scale measure. Sets the respective fields accordingly.
     */
    private void calcAreaAspectRatioScaleMeasure() {
        if (this.drawingWidth > 0 && this.drawingHeight > 0) {
            this.area = this.drawingWidth * this.drawingHeight;
            this.aspectRatio = this.drawingWidth / this.drawingHeight;
            this.scaleMeasure = DrawingUtil.computeScaleMeasure(this.drawingWidth, this.drawingHeight, this.dar);
        }
    }

    //////////////////////////////////////////////////////////////////
    // Getters and setters.

    /**
     * Gets the drawing width.
     */
    public double getDrawingWidth() {
        return drawingWidth;
    }

    /**
     * Sets the drawing width.
     */
    public void setDrawingWidth(final double drawingWidth) {
        this.drawingWidth = drawingWidth;
        calcAreaAspectRatioScaleMeasure();
    }

    /**
     * Gets the drawing height.
     */
    public double getDrawingHeight() {
        return drawingHeight;
    }

    /**
     * Sets the drawing height.
     */
    public void setDrawingHeight(final double drawingHeight) {
        this.drawingHeight = drawingHeight;
        calcAreaAspectRatioScaleMeasure();
    }

    /**
     * Returns the scale measure.
     */
    public double getScaleMeasure() {
        return scaleMeasure;
    }

    /**
     * Gets the area.
     */
    public double getArea() {
        return area;
    }

    /**
     * Gets the aspect ratio.
     */
    public double getAspectRatio() {
        return aspectRatio;
    }

    /**
     * Gets placement option.
     */
    public DrawingDataDescriptor getPlacementOption() {
        return placementOption;
    }

    /**
     * Sets placement option.
     */
    public void setPlacementOption(final DrawingDataDescriptor placementOption) {
        this.placementOption = placementOption;
    }

    /**
     * Gets potential x-coordinate of the rectangle to place.
     */
    public double getNextXcoordinate() {
        return nextXcoordinate;
    }

    /**
     * Sets potential x-coordinate of the rectangle to place.
     */
    public void setNextXcoordinate(final double potentialX) {
        this.nextXcoordinate = potentialX;
    }

    /**
     * Gets potential y-coordinate of the rectangle to place.
     */
    public double getNextYcoordinate() {
        return nextYcoordinate;
    }

    /**
     * Sets potential y-coordinate of rectangle to place.
     */
    public void setNextYcoordinate(final double potentialY) {
        this.nextYcoordinate = potentialY;
    }

    /**
     * @return the desired aspect ratio
     */
    public double getDesiredAspectRatio() {
        return dar;
    }
}
