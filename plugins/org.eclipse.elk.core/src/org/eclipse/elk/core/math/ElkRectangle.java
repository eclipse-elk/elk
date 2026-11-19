/*******************************************************************************
 * Copyright (c) 2014, 2017 Kiel University and others.
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
package org.eclipse.elk.core.math;

import java.util.Objects;

/**
 * Class resembles basic functionality of {@link java.awt.geom.Rectangle2D}. This way it is possible
 * to avoid awt dependencies in ELK Layered's code.
 */
public class ElkRectangle {

    // CHECKSTYLEOFF VisibilityModifier

    /** The X coordinate of this <code>Rectangle</code>. */
    public double x;

    /** The Y coordinate of this <code>Rectangle</code>. */
    public double y;

    /** The width of this <code>Rectangle</code>. */
    public double width;

    /** The height of this <code>Rectangle</code>. */
    public double height;

    /**
     * Constructs a new <code>Rectangle</code>, initialized to location (0,&nbsp;0) and size
     * (0,&nbsp;0).
     */
    public ElkRectangle() {
        this(0, 0, 0, 0);
    }

    /**
     * Constructs and initializes a <code>Rectangle</code> from the specified <code>double</code>
     * coordinates.
     * 
     * @param x
     *            the x coordinate of the upper-left corner of the newly constructed
     *            <code>Rectangle</code>
     * @param y
     *            the y coordinate of the upper-left corner of the newly constructed
     *            <code>Rectangle</code>
     * @param w
     *            the width of the newly constructed <code>Rectangle</code>
     * @param h
     *            the height of the newly constructed <code>Rectangle</code>
     */
    public ElkRectangle(final double x, final double y, final double w, final double h) {
        this.x = x;
        this.y = y;
        this.width = w;
        this.height = h;
    }

    /**
     * Constructs and initializes a <code>Rectangle</code> from the specified instance.
     * 
     * @param rect
     *            the existing rectangle whose values to copy.
     */
    public ElkRectangle(final ElkRectangle rect) {
        this.x = rect.x;
        this.y = rect.y;
        this.width = rect.width;
        this.height = rect.height;
    }

    /**
     * Sets the location and size of this <code>Rectangle</code> to the specified
     * <code>double</code> values.
     * 
     * @param nx
     *            the X coordinate of the upper-left corner of this <code>Rectangle</code>
     * @param ny
     *            the Y coordinate of the upper-left corner of this <code>Rectangle</code>
     * @param nw
     *            the width of this <code>Rectangle</code>
     * @param nh
     *            the height of this <code>Rectangle</code>
     */
    public void setRect(final double nx, final double ny, final double nw, final double nh) {
        this.x = nx;
        this.y = ny;
        this.width = nw;
        this.height = nh;
    }

    /**
     * @return the (x,y) position of this rectangle.
     */
    public KVector getPosition() {
        return new KVector(x, y);
    }
    
    /**
     * @return the top left coordinate (x,y).
     */
    public KVector getTopLeft() {
        return getPosition();
    }

    /**
     * @return the top right coordinate (x+w,y).
     */
    public KVector getTopRight() {
        return new KVector(x + width, y);
    }

    /**
     * @return the bottom left coordinate (x,y+h).
     */
    public KVector getBottomLeft() {
        return new KVector(x, y + height);
    }

    /**
     * @return the bottom right coordinate (x+w,y+h).
     */
    public KVector getBottomRight() {
        return new KVector(x + width, y + height);
    }
    
    /**
     * @return the center point of a rectangle.
     */
    public KVector getCenter() {
        return new KVector(x + width / 2, y + height / 2);
    }
    
    /**
     * Unions the receiver and the given <code>Rectangle</code> objects and puts the result into the
     * receiver.
     * 
     * @param other
     *            the <code>Rectangle</code> to be combined with this instance
     */
    public void union(final ElkRectangle other) {
        double x1 = Math.min(this.x, other.x);
        double y1 = Math.min(this.y, other.y);
        double x2 = Math.max(this.x + this.width, other.x + other.width);
        double y2 = Math.max(this.y + this.height, other.y + other.height);
        if (x2 < x1) {
            double t = x1;
            x1 = x2;
            x2 = t;
        }
        if (y2 < y1) {
            double t = y1;
            y1 = y2;
            y2 = t;
        }
        setRect(x1, y1, x2 - x1, y2 - y1);
    }
    
    /**
     * Moves the rectangle by the given offset.
     * 
     * @param offset
     *            the offset to move the rectangle by.
     */
    public void move(final KVector offset) {
        x += offset.x;
        y += offset.y;
    }
    
    @Override
    public String toString() {
        return "Rect[x=" + x + ",y=" + y + ",w=" + width + ",h=" + height + "]";
    }

    @Override
    public boolean equals(final Object obj) {
        if (obj == null || !(obj instanceof ElkRectangle)) {
            return false;
        }
        ElkRectangle other = (ElkRectangle) obj;
        return Objects.equals(x, other.x) && Objects.equals(y, other.y)
                && Objects.equals(width, other.width) && Objects.equals(height, other.height);
    }
    
    @Override
    public int hashCode() {
        return Objects.hash(x, y, width, height);
    }

    /**
    * Returns the largest X coordinate of the rectangle in double precision.
    */
    public double getMaxX() {
        return this.x + this.width;
    }
    
    /**
    * Returns the largest Y coordinate of the rectangle in double precision.
    */
    public double getMaxY() {
        return this.y + this.height;
    }

    /**
     * Tests if the interior of this rect intersects the interior of a specified rect.
     * 
     * @param rect the other rectangle.
     */
    public boolean intersects(final ElkRectangle rect) {
        double r1x1 = this.x;
        double r1y1 = this.y;
        double r1x2 = this.x + this.width;
        double r1y2 = this.y + this.height;
        double r2x1 = rect.x;
        double r2y1 = rect.y;
        double r2x2 = rect.x + rect.width;
        double r2y2 = rect.y + rect.height;

        return r1x1 < r2x2 && r1x2 > r2x1 && r1y2 > r2y1 && r1y1 < r2y2;
    }

}
