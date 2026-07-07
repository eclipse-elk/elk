/*******************************************************************************
 * Copyright (c) 2014, 2015 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.layered.p5edges.splines;

import java.text.DecimalFormat;
import java.util.Collection;
import java.util.Iterator;

import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LShape;
import org.eclipse.elk.core.math.ElkMargin;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.options.PortSide;

import com.google.common.collect.Iterables;

/**
 * This class is representing a rectangle and provides some comfort methods to work with the
 * rectangle. Special about this implementation is that there is no width and height stored. This
 * improves the performance of unions and accesses to the right/top side.
 * 
 * @author tit
 */
public final class Rectangle {
    
    // The rectangle is defined by the four outer coordinates.
    /** Minimum y-value. */
    private double top = Double.MAX_VALUE; 
    /** Maximum y-value. */
    private double bottom = -Double.MAX_VALUE; 
    /** Maximum x-value. */
    private double left = Double.MAX_VALUE;
    /** Minimum x-value. */
    private double right = -Double.MAX_VALUE; 
    
    /**
     * Creates a new rectangle with the given values.
     * 
     * @param left The minimum x-value of the rectangle.
     * @param top The maximum y-value of the rectangle.
     * @param right The maximum x-value of the rectangle.
     * @param bottom The minimum y-value of the rectangle.
     */
    public Rectangle(final double left, final double top, final double right, final double bottom) {
        if (top > bottom) {
            throw new IllegalArgumentException("Top must be smaller or equal to bottom.");
        } else if (left > right) {
            throw new IllegalArgumentException("Left must be smaller or equal to right.");
        }
        this.top = top;
        this.right = right;
        this.bottom = bottom;
        this.left = left;
    }

    /**
     * A copy constructor.
     * 
     * @param rectangle The rectangle to copy.
     */
    public Rectangle(final Rectangle rectangle) {
        top = rectangle.top;
        right = rectangle.right;
        bottom = rectangle.bottom;
        left = rectangle.left;
    }
    
    /**
     * Creates a new Rectangle for the given KVector, resulting in a rectangle with height = width = 0.
     * 
     * @param position The position of the rectangle.
     */
    public Rectangle(final KVector position) {
        top = position.y;
        bottom = position.y;
        right = position.x;
        left = position.x;
    }
    
    /**
     * Creates a rectangle around the given shape.
     * 
     * @param shape The shape to form the rectangle around.
     */
    public Rectangle(final LShape shape) {
        final KVector position = shape.getPosition();
        final KVector extend = position.clone().add(shape.getSize());
        
        top = Math.min(position.y, extend.y);
        bottom = Math.max(position.y, extend.y);
        left = Math.min(position.x, extend.x);
        right = Math.max(position.x, extend.x);
    }

    /**
     * Creates a rectangle around all given coordinates.
     * 
     * @param vectors The coordinates to span the rectangle around.
     */
    public Rectangle(final Iterable<KVector> vectors) {
        if (Iterables.isEmpty(vectors)) {
            throw new IllegalArgumentException("The list of vectors may not be empty.");
        }
            
        for (final KVector vector : vectors) {
            top = Math.min(top, vector.y); 
            right = Math.max(right, vector.x); 
            bottom = Math.max(bottom, vector.y); 
            left = Math.min(left, vector.x); 
        }
    }
    
    /**
     * Creates a new rectangle including all coordinates of the given vectors.
     * 
     * @param vectors The vectors to span a rectangle around.
     */
    public Rectangle(final KVector... vectors) {
        if (vectors.length == 0) {
            throw new IllegalArgumentException("The list of vectors may not be empty.");
        }
            
        for (final KVector vector : vectors) {
            top = Math.min(top, vector.y); 
            right = Math.max(right, vector.x); 
            bottom = Math.max(bottom, vector.y); 
            left = Math.min(left, vector.x); 
        }
    }
    
    /**
     * Enlarges this rectangle to include the given coordinate.
     * 
     * @param vector The vector to the coordinate to include.
     */
    public void union(final KVector vector) {
        top = Math.min(top, vector.y); 
        right = Math.max(right, vector.x); 
        bottom = Math.max(bottom, vector.y); 
        left = Math.min(left, vector.x); 
    }

    /**
     * Enlarges this rectangle to include the given shape.
     * 
     * @param shape The shape to include.
     */
    public void union(final LShape shape) {
        this.union(shape.getPosition());
        this.union(shape.getPosition().clone().add(shape.getSize().x, shape.getSize().y));
    }

    /**
     * Enlarges this rectangle to include the given rectangle.
     * 
     * @param rectangle The rectangle to include.
     */
    public void union(final Rectangle rectangle) {
        top = Math.min(top, rectangle.top); 
        right = Math.max(right, rectangle.right); 
        bottom = Math.max(bottom, rectangle.bottom); 
        left = Math.min(left, rectangle.left); 
    }

    /**
     * Creates a new rectangle, spanning around the two given rectangles.
     * 
     * @param rectangle1 The first rectangle to include.
     * @param rectangle2 The second rectangle to include.
     * @return The minimum rectangle, that includes both input rectangles.
     */
    public static Rectangle union(final Rectangle rectangle1, final Rectangle rectangle2) {
        final Rectangle retVal = new Rectangle(rectangle1);
        retVal.union(rectangle2);
        return retVal;
    }

    /**
     * Creates a new rectangle, spanning around the given rectangles.
     * 
     * @param rectangles The rectangles to include.
     * @return The minimum rectangle, that includes all input rectangles.
     */
    public static Rectangle union(final Collection<Rectangle> rectangles) {
        final Iterator<Rectangle> iter = rectangles.iterator();
        if (iter.hasNext()) {
            final Rectangle rectangle = iter.next();
            final Rectangle retVal = new Rectangle(rectangle);
            while (iter.hasNext()) {
                rectangle.union(iter.next());
            }
            return retVal;
        } else {
            throw new IllegalArgumentException("The list of vectors may not be null.");
        }
    }
    
    /**
     * Enlarges this rectangle in all four directions by the given value.
     * 
     * @param enlargement The value this rectangle will be enlarged in all four directions.
     */
    public void enlarge(final double enlargement) {
        top = top - enlargement;
        left = left - enlargement;
        right = right + enlargement;
        bottom = bottom + enlargement;
    }

    /**
     * Calculates and returns the height of this rectangle.
     * 
     * @return The height.
     */
    public double getHeight() {
        return bottom - top;
    }

    /**
     * Calculates and returns the width of this rectangle.
     * 
     * @return The width.
     */
    public double getWidth() {
        return right - left;
    }
    
    /**
     * Returns the top edge of this rectangle.
     * 
     * @return The top edge.
     */
    public double getTop() {
        return top;
    }

    /**
     * Returns the right edge of this rectangle.
     * 
     * @return The right edge.
     */
    public double getRight() {
        return right;
    }

    /**
     * Returns the bottom edge of this rectangle.
     * 
     * @return The bottom edge.
     */
    public double getBottom() {
        return bottom;
    }

    /**
     * Returns the left edge of this rectangle.
     * 
     * @return The left edge.
     */
    public double getLeft() {
        return left;
    }
    
    /**
     * Returns the coordinate of the edge of this rectangle that is corresponding to the 
     * given {@link PortSide}.
     * 
     * @param side The {@link PortSide} for which to return the edge of this rectangle.
     * @return The position of the specified edge.
     */
    public double getFromPortSide(final PortSide side) {
        switch (side) {
        case NORTH:
            return top;
        case EAST:
            return right;
        case SOUTH:
            return bottom;
        case WEST:
            return left;
        default:
            return 0.0;
        }
    }

    /**
     * Creates a margin that is relative to a {@link LShape}, so that this {@link Rectangle} completely
     * lays inside the margins and is of minimal size. The margins of any side of the shape will be at 
     * least 0. If the given shape is a {@link LNode} it's margin is neither modified nor taken into 
     * account. 
     * 
     * @param shape The shape the margins shall lay around.
     * @return A new {@link Margins} laying around the shape.
     */
    public ElkMargin toNodeMargins(final LShape shape) {
        final ElkMargin retVal =  new ElkMargin();

        final Rectangle shapeRectangle = new Rectangle(0, 0, shape.getSize().x, shape.getSize().y);
        retVal.top = Math.max(0.0, shapeRectangle.top - top);
        retVal.left = Math.max(0.0, shapeRectangle.left - left);
        retVal.bottom = Math.max(0.0, bottom - shapeRectangle.bottom);
        retVal.right = Math.max(0.0, right - shapeRectangle.right);
        
        return retVal;
    }
    
    // elkjs-exclude-start
    // DecimalFormat not available in gwt
    /** Format for the values shown on toString(). */
    private static final DecimalFormat DEC_FORMAT = new DecimalFormat("#0.0");

    /**
     * Converts this rectangle to a readable string with rounded values.
     * 
     * @return The readable string.
     */
    public String toString() {
        return "[" 
        + "top= " + DEC_FORMAT.format(top) 
        + ",left= " + DEC_FORMAT.format(left) 
        + ",bottom= " + DEC_FORMAT.format(bottom) 
        + ",right= " + DEC_FORMAT.format(right) + "]";
    }
    // elkjs-exclude-end
}
