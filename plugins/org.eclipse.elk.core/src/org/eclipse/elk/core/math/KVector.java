/*******************************************************************************
 * Copyright (c) 2010, 2015 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.core.math;

import java.util.Random;

import org.eclipse.elk.core.util.IDataObject;

/**
 * A simple 2D vector class which supports translation, scaling, normalization etc.
 *
 * @author uru
 * @author owo
 * @author cds
 */
public final class KVector implements IDataObject, Cloneable {

    /** the default fuzzyness used when comparing two vectors fuzzily. */
    private static final double DEFAULT_FUZZYNESS = 0.05;

    /** the serial version UID. */
    private static final long serialVersionUID = -4780985519832787684L;

    // CHECKSTYLEOFF VisibilityModifier
    /** x coordinate. */
    public double x;
    /** y coordinate. */
    public double y;
    // CHECKSTYLEON VisibilityModifier

    /**
     * Create vector with default coordinates (0,0).
     */
    public KVector() {
        this.x = 0.0;
        this.y = 0.0;
    }

    /**
     * Constructs a new vector from given values.
     *
     * @param thex
     *            x value
     * @param they
     *            y value
     */
    public KVector(final double thex, final double they) {
        this.x = thex;
        this.y = they;
    }

    /**
     * Creates an exact copy of a given vector v.
     *
     * @param v
     *            existing vector
     */
    public KVector(final KVector v) {
        this.x = v.x;
        this.y = v.y;
    }
    
    /**
     * Creates a new vector that points from the given start to the end vector.
     * 
     * @param start
     *            start vector.
     * @param end
     *            end vector.
     */
    public KVector(final KVector start, final KVector end) {
        this.x = end.x - start.x;
        this.y = end.y - start.y;
    }

    /**
     * Creates a normalized vector for the passed angle in radians.
     *
     * @param angle
     *            angle in radians.
     */
    public KVector(final double angle) {
        this.x = Math.cos(angle);
        this.y = Math.sin(angle);
    }

    /**
     * Returns an exact copy of this vector.
     *
     * @return identical vector
     */
    // elkjs-exclude-start
    @Override
    // elkjs-exclude-end
    public KVector clone() {
        return new KVector(x, y);
    }

    @Override
    public String toString() {
        return "(" + x + "," + y + ")";
    }

    @Override
    public boolean equals(final Object obj) {
        if (obj instanceof KVector) {
            KVector other = (KVector) obj;
            return this.x == other.x && this.y == other.y;
        } else {
            return false;
        }
    }

    /**
     * Calls {@link #equalsFuzzily(KVector, double)} with a default fuzzyness.
     *
     * @param other
     *            the vector to compare this vector to.
     * @return {@code true} if the vectors are approximately equal, {@code false} otherwise.
     */
    public boolean equalsFuzzily(final KVector other) {
        return equalsFuzzily(other, DEFAULT_FUZZYNESS);
    }

    /**
     * Compares if this and the given vector are approximately equal. What <i>approximately</i>
     * means is defined by the fuzzyness: for both x and y coordinate, the two vectors may only
     * differ by at most the fuzzyness to still be considered equal.
     *
     * @param other
     *            the vector to compare this vector to.
     * @param fuzzyness
     *            the maximum difference per dimension that is still considered equal.
     * @return {@code true} if the vectors are approximately equal, {@code false} otherwise.
     */
    public boolean equalsFuzzily(final KVector other, final double fuzzyness) {
        return Math.abs(this.x - other.x) <= fuzzyness
                && Math.abs(this.y - other.y) <= fuzzyness;
    }

    @Override
    public int hashCode() {
        return Double.valueOf(x).hashCode() + Integer.reverse(Double.valueOf(y).hashCode());
    }

    /**
     * returns this vector's length.
     *
     * @return Math.sqrt(x*x + y*y)
     */
    public double length() {
        return Math.sqrt(x * x + y * y);
    }

    /**
     * returns square length of this vector.
     *
     * @return x*x + y*y
     */
    public double squareLength() {
        return x * x + y * y;
    }

    /**
     * Set vector to (0,0).
     *
     * @return {@code this}
     */
    public KVector reset() {
        this.x = 0.0;
        this.y = 0.0;
        return this;
    }

    /**
     * Resets this vector to the value of the other vector.
     *
     * @param other the vector whose values to copy.
     * @return {@code this}
     */
    public KVector set(final KVector other) {
        this.x = other.x;
        this.y = other.y;
        return this;
    }

    /**
     * Resets this vector to the given values.
     *
     * @param newX new x value.
     * @param newY new y value.
     * @return {@code this}
     */
    public KVector set(final double newX, final double newY) {
        this.x = newX;
        this.y = newY;
        return this;
    }

    /**
     * Vector addition.
     *
     * @param v
     *            vector to add
     * @return <code>this + v</code>
     */
    public KVector add(final KVector v) {
        this.x += v.x;
        this.y += v.y;
        return this;
    }

    /**
     * Translate the vector by adding the given amount.
     *
     * @param dx
     *            the x offset
     * @param dy
     *            the y offset
     * @return {@code this}
     */
    public KVector add(final double dx, final double dy) {
        this.x += dx;
        this.y += dy;
        return this;
    }

    /**
     * Returns the sum of arbitrarily many vectors as a new vector instance.
     *
     * @param vs vectors to be added
     * @return a new vector containing the sum of given vectors
     */
    public static KVector sum(final KVector... vs) {
        KVector sum = new KVector();
        for (KVector v : vs) {
            sum.x += v.x;
            sum.y += v.y;
        }
        return sum;
    }

    /**
     * Vector subtraction.
     *
     * @param v
     *            vector to subtract
     * @return {@code this}
     */
    public KVector sub(final KVector v) {
        this.x -= v.x;
        this.y -= v.y;
        return this;
    }

    /**
     * Translate the vector by subtracting the given amount.
     *
     * @param dx
     *            the x offset
     * @param dy
     *            the y offset
     * @return {@code this}
     */
    public KVector sub(final double dx, final double dy) {
        this.x -= dx;
        this.y -= dy;
        return this;
    }

    /**
     * Returns the difference of two vectors as a new vector instance.
     *
     * @param v1
     *            the minuend
     * @param v2
     *            the subtrahend
     * @return a new vector containing the difference of given vectors
     */
    public static KVector diff(final KVector v1, final KVector v2) {
        return new KVector(v1.x - v2.x, v1.y - v2.y);
    }

    /**
     * Scale the vector.
     *
     * @param scale
     *            scaling factor
     * @return {@code this}
     */
    public KVector scale(final double scale) {
        this.x *= scale;
        this.y *= scale;
        return this;
    }

    /**
     * Scale the vector with different values for X and Y coordinate.
     *
     * @param scalex
     *            the x scaling factor
     * @param scaley
     *            the y scaling factor
     * @return {@code this}
     */
    public KVector scale(final double scalex, final double scaley) {
        this.x *= scalex;
        this.y *= scaley;
        return this;
    }

    /**
     * Normalize the vector.
     *
     * @return {@code this}
     */
    public KVector normalize() {
        double length = this.length();
        if (length > 0) {
            this.x /= length;
            this.y /= length;
        }
        return this;
    }

    /**
     * scales this vector to the passed length.
     *
     * @param length
     *            length to scale to
     * @return {@code this}
     */
    public KVector scaleToLength(final double length) {
        this.normalize();
        this.scale(length);
        return this;
    }

    /**
     * Negate the vector.
     *
     * @return {@code this}
     */
    public KVector negate() {
        this.x = -this.x;
        this.y = -this.y;
        return this;
    }

    /**
     * Returns angle representation of this vector in degree. The length of the vector must not be 0.
     *
     * @return value within [0,360)
     */
    public double toDegrees() {
        return Math.toDegrees(toRadians());
    }

    /**
     * Returns angle representation of this vector in radians. The length of the vector must not be 0.
     *
     * @return value within [0,2*pi)
     */
    public double toRadians() {
        double length = this.length();
        assert length > 0;

        if (x >= 0 && y >= 0) {  // 1st quadrant
            return Math.asin(y / length);
        } else if (x < 0) {      // 2nd or 3rd quadrant
            return Math.PI - Math.asin(y / length);
        } else {                 // 4th quadrant
            return 2 * Math.PI + Math.asin(y / length);
        }
    }

    /**
     * Add some "noise" to this vector.
     *
     * @param random
     *            the random number generator
     * @param amount
     *            the amount of noise to add
     */
    public void wiggle(final Random random, final double amount) {
        this.x += random.nextDouble() * amount - (amount / 2);
        this.y += random.nextDouble() * amount - (amount / 2);
    }

    /**
     * Returns the distance between two vectors.
     *
     * @param v2
     *            second vector
     * @return distance between this and second vector
     */
    public double distance(final KVector v2) {
        double dx = this.x - v2.x;
        double dy = this.y - v2.y;
        return Math.sqrt((dx * dx) + (dy * dy));
    }

    /**
     * Returns the dot product of the two given vectors.
     *
     * @param v2
     *            second vector
     * @return this.x * v2.x + this.y * v2.y
     */
    public double dotProduct(final KVector v2) {
        return this.x * v2.x + this.y * v2.y;
    }
    
    /**
     * Calculates the cross product of two vectors v and w.
     * @param v
     * @param w
     * @return the cross product of v and w
     */
    public static double crossProduct(final KVector v, final KVector w) {
        return v.x * w.y - v.y * w.x;
    }
    
    /**
     * Rotates the given vector v by the angle provided in radians. Rotation is anti-clockwise.
     * @param v
     * @param angle
     * @return the rotated vector
     */
    public KVector rotate(final double angle) {
        double newX = this.x * Math.cos(angle) - this.y * Math.sin(angle);
        this.y = this.x * Math.sin(angle) + this.y * Math.cos(angle);
        this.x = newX;
        return this;
    }
    
    /**
     * Returns the angle between this vector and another given vector in radians.
     * @param other
     * @return angle between vectors
     */
    public double angle(KVector other) {
        return Math.acos(this.dotProduct(other) / (this.length() * other.length()));
    }

    /**
     * Apply the given bounds to this vector.
     *
     * @param lowx
     *            the lower bound for x coordinate
     * @param lowy
     *            the lower bound for y coordinate
     * @param highx
     *            the upper bound for x coordinate
     * @param highy
     *            the upper bound for y coordinate
     * @return {@code this}
     * @throws IllegalArgumentException
     *             if highx < lowx or highy < lowy
     */
    public KVector bound(final double lowx, final double lowy, final double highx,
            final double highy) {
        if (highx < lowx || highy < lowy) {
            throw new IllegalArgumentException(
                    "The highx must be bigger then lowx and the highy must be bigger then lowy");
        }
        if (x < lowx) {
            x = lowx;
        } else if (x > highx) {
            x = highx;
        }
        if (y < lowy) {
            y = lowy;
        } else if (y > highy) {
            y = highy;
        }
        return this;
    }

    /**
     * Determine whether any of the two values are NaN.
     *
     * @return true if x is NaN or y is NaN
     */
    public boolean isNaN() {
        return Double.isNaN(x) || Double.isNaN(y);
    }

    /**
     * Determine whether any of the two values are infinite.
     *
     * @return true if x is infinite or y is infinite
     */
    public boolean isInfinite() {
        return Double.isInfinite(x) || Double.isInfinite(y);
    }

    @Override
    public void parse(final String string) {
        int start = 0;
        while (start < string.length() && isdelim(string.charAt(start), "([{\"' \t\r\n")) {
            start++;
        }
        int end = string.length();
        while (end > 0 && isdelim(string.charAt(end - 1), ")]}\"' \t\r\n")) {
            end--;
        }
        if (start >= end) {
            throw new IllegalArgumentException("The given string does not contain any numbers.");
        }
        String[] tokens = string.substring(start, end).split(",|;|\r|\n");
        if (tokens.length != 2) {
            throw new IllegalArgumentException("Exactly two numbers are expected, "
                    + tokens.length + " were found.");
        }
        try {
            x = Double.parseDouble(tokens[0].trim());
            y = Double.parseDouble(tokens[1].trim());
        } catch (NumberFormatException exception) {
            throw new IllegalArgumentException(
                    "The given string contains parts that cannot be parsed as numbers." + exception);
        }
    }

    /**
     * Determine whether the given character is a delimiter.
     *
     * @param c
     *            a character
     * @param delims
     *            a string of possible delimiters
     * @return true if {@code c} is one of the characters in {@code delims}
     */
    private static boolean isdelim(final char c, final String delims) {
        for (int i = 0; i < delims.length(); i++) {
            if (c == delims.charAt(i)) {
                return true;
            }
        }
        return false;
    }

}
