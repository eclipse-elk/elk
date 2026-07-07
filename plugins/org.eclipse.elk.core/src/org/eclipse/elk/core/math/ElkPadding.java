/*******************************************************************************
 * Copyright (c) 2016 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.core.math;

/**
 * Stores the padding of an element. The padding is spacings from the border of an element to
 * other internal elements.
 */
public class ElkPadding extends Spacing {

    /** The serial version UID. */
    private static final long serialVersionUID = -2159860709896900657L;

    /**
     * Creates a new instance with all fields set to {@code 0.0}.
     */
    public ElkPadding() {
        super();
    }
    
    /**
     * Create new padding using the same value for every side.
     * 
     * @param any
     *            padding value for every side.
     */
    public ElkPadding(final double any) {
        super(any, any, any, any);
    }

    /**
     * Create new padding using the same padding for both sides of a common dimension.
     * 
     * @param leftRight
     *            value used for the left and right padding.
     * @param topBottom
     *            value used for the top and bottom padding.
     */
    public ElkPadding(final double leftRight, final double topBottom) {
        super(topBottom, leftRight, topBottom, leftRight);
    }

    /**
     * Creates a new instance initialized with the given values.
     * 
     * @param top
     *            the padding from the top.
     * @param right
     *            the padding from the right.
     * @param bottom
     *            the padding from the bottom.
     * @param left
     *            the padding from the left.
     */
    public ElkPadding(final double top, final double right, final double bottom, final double left) {
        super(top, right, bottom, left);
    }
    
    /**
     * Creates a new instance with all fields set to the value of {@code other}.
     * 
     * @param other
     *            paddings object from which to copy the values.
     */
    public ElkPadding(final ElkPadding other) {
        super(other.top, other.right, other.bottom, other.left);
    }
    
    // elkjs-exclude-start
    @Override
    // elkjs-exclude-end
    public ElkPadding clone() {
        return new ElkPadding(this);
    }
}
