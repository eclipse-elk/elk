/*******************************************************************************
 * Copyright (c) 2014, 2015 Kiel University and others.
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

import org.eclipse.elk.core.util.IDataObject;

/**
 * Stores the spacing of an object in {@code double} precision.
 * 
 * @author cds
 * @author uru
 */
public abstract class Spacing implements IDataObject, Cloneable {
    // Allow public fields in these utility classes.
    // CHECKSTYLEOFF VisibilityModifier

    /** The serial version UID. */
    private static final long serialVersionUID = 4358555478195088364L;

    /** The spacing from the top. */
    public double top = 0.0;
    /** The spacing from the bottom. */
    public double bottom = 0.0;
    /** The spacing from the left. */
    public double left = 0.0;
    /** The spacing from the right. */
    public double right = 0.0;

    // CHECKSTYLEON VisibilityModifier
    
    /**
     * Creates a new instance with all fields set to {@code 0.0}.
     */
    protected Spacing() { }

    /**
     * Creates a new instance initialized with the given values.
     * 
     * @param top
     *            the spacing from the top.
     * @param right
     *            the spacing from the right.
     * @param bottom
     *            the spacing from the bottom.
     * @param left
     *            the spacing from the left.
     */
    protected Spacing(final double top, final double right, final double bottom, final double left) {
        set(top, right, bottom, left);
    }
    
    /**
     * Sets all spacings to the ones of the given {@link Spacing} object.
     * 
     * @param spacing
     *            the spacing object to set this instance's values to.
     */
    public void set(final Spacing spacing) {
        set(spacing.top, spacing.right, spacing.bottom, spacing.left);
    }

    /**
     * Sets all four spacings at once.
     * 
     * @param newTop
     *            the spacing from the top.
     * @param newLeft
     *            the spacing from the left.
     * @param newBottom
     *            the spacing from the bottom.
     * @param newRight
     *            the spacing from the right.
     */
    public void set(final double newTop, final double newRight, final double newBottom, final double newLeft) {
        this.top = newTop;
        this.right = newRight;
        this.bottom = newBottom;
        this.left = newLeft;
    }
    
    /**
     * @return the top
     */
    public double getTop() {
        return top;
    }
    
    /**
     * @param top the top to set
     */
    public void setTop(final double top) {
        this.top = top;
    }
    
    /**
     * @return the right
     */
    public double getRight() {
        return right;
    }
    
    /**
     * @param right the right to set
     */
    public void setRight(final double right) {
        this.right = right;
    }
    
    /**
     * @return the bottom
     */
    public double getBottom() {
        return bottom;
    }
    
    /**
     * @param bottom the bottom to set
     */
    public void setBottom(final double bottom) {
        this.bottom = bottom;
    }
    
    /**
     * @return the left
     */
    public double getLeft() {
        return left;
    }
    
    /**
     * @param left the left to set
     */
    public void setLeft(final double left) {
        this.left = left;
    }
    
    /**
     * Sets both the left and the right spacing to the same value.
     * 
     * @param val the new value.
     */
    public void setLeftRight(final double val) {
        this.left = val;
        this.right = val;
    }
    
    /**
     * Sets both the top and the bottom spacing to the same value.
     * 
     * @param val the new value.
     */
    public void setTopBottom(final double val) {
        this.top = val;
        this.bottom = val;
    }
    
    /** 
     * @return the combined horizontal padding, i.e. {@link #getLeft()} + {@link #getRight()}.   
     */ 
    public double getHorizontal() { 
        return this.getLeft() + this.getRight(); 
    } 
     
    /** 
     * @return the combined vertical padding, i.e. {@link #getTop()} + {@link #getBottom()}.   
     */ 
    public double getVertical() { 
        return this.getTop() + this.getBottom(); 
    } 

    @Override
    public boolean equals(final Object obj) {
        if (obj instanceof Spacing) {
            Spacing other = (Spacing) obj;
            return this.top == other.top && this.bottom == other.bottom && this.left == other.left
                    && this.right == other.right;
        } else {
            return false;
        }
    }

    @Override
    public int hashCode() {
        // I don't care very much to define explicit constants for this calculation.
        // CHECKSTYLEOFF MagicNumber

        int code1 = java.lang.Double.valueOf(left).hashCode() << 16;
        code1 |= java.lang.Double.valueOf(bottom).hashCode() & 0xffff;

        int code2 = java.lang.Double.valueOf(right).hashCode() << 16;
        code2 |= java.lang.Double.valueOf(top).hashCode() & 0xffff;

        return code1 ^ code2;

        // CHECKSTYLEON MagicNumber
    }

    @Override
    public String toString() {
        return "[top=" + top + ",left=" + left + ",bottom=" + bottom + ",right=" + right + "]";
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
        if (start < end) {
            String[] tokens = string.substring(start, end).split(",|;");
            try {
                for (String token : tokens) {
                    String[] keyandvalue = token.split("=");
                    if (keyandvalue.length != 2) {
                        throw new IllegalArgumentException("Expecting a list of key-value pairs.");
                    }
                    String key = keyandvalue[0].trim();
                    double value = Double.parseDouble(keyandvalue[1].trim());
                    if (key.equals("top")) {
                        top = value;
                    } else if (key.equals("left")) {
                        left = value;
                    } else if (key.equals("bottom")) {
                        bottom = value;
                    } else if (key.equals("right")) {
                        right = value;
                    }
                }
            } catch (NumberFormatException exception) {
                throw new IllegalArgumentException(
                        "The given string contains parts that cannot be parsed as numbers." + exception);
            }
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

    /**
     * Copy spacings values from another double valued spacings.
     * 
     * @param other
     *            another spacings
     * @return this instance.
     */
    public Spacing copy(final Spacing other) {
        this.left = other.left;
        this.right = other.right;
        this.top = other.top;
        this.bottom = other.bottom;
        return this;
    }
    
    /**
     * Adds spacings values from another double valued spacings.
     * 
     * @param other
     *            another spacings.
     * @return this instance.
     */
    public Spacing add(final Spacing other) {
        this.left += other.left;
        this.right += other.right;
        this.top += other.top;
        this.bottom += other.bottom;
        return this;
    }

    /**
     * Return a new instance of the concrete class. 
     */
    public abstract Spacing clone();
}
