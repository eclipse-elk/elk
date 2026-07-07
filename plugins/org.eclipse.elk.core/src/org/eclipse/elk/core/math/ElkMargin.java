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
 * Stores the margins of an element. The margin is the area around the border of an element that has to be kept free
 * of any other elements.
 */
public class ElkMargin extends Spacing {

    /** The serial version UID. */
    private static final long serialVersionUID = 7465583871643915474L;

    /**
     * Creates a new instance with all fields set to {@code 0.0}.
     */
    public ElkMargin() {
        super();
    }
    
    /**
     * Create new margin using the same value for every side.
     * 
     * @param any
     *            margin value for every side.
     */
    public ElkMargin(final double any) {
        super(any, any, any, any);
    }

    /**
     * Create new margin using the same margin for both sides of a common dimension.
     * 
     * @param leftRight
     *            value used for the left and right margin.
     * @param topBottom
     *            value used for the top and bottom margin.
     */
    public ElkMargin(final double leftRight, final double topBottom) {
        super(topBottom, leftRight, topBottom, leftRight);
    }

    /**
     * Creates a new instance initialized with the given values.
     * 
     * @param top
     *            the margin from the top.
     * @param right
     *            the margin from the right.
     * @param bottom
     *            the margin from the bottom.
     * @param left
     *            the margin from the left.
     */
    public ElkMargin(final double top, final double right, final double bottom, final double left) {
        super(top, right, bottom, left);
    }

    /**
     * Creates a new instance with all fields set to the value of {@code other}.
     * 
     * @param other
     *            margins object from which to copy the values.
     */
    public ElkMargin(final ElkMargin other) {
        super(other.top, other.right, other.bottom, other.left);
    }
    
    // elkjs-exclude-start
    @Override
    // elkjs-exclude-end
    public ElkMargin clone() {
        return new ElkMargin(this);
    }
    
}
