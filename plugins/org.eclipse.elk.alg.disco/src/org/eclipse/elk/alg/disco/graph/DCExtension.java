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
package org.eclipse.elk.alg.disco.graph;

import org.eclipse.elk.core.math.ElkRectangle;
import org.eclipse.elk.core.math.KVector;

/**
 * Models short hierarchical edges conceptually as semi-indefinite strips ending in a {@link DCElement}.
 */
public class DCExtension {
    private DCElement parent;
    private DCDirection direction;
    // The offset is not centered in the extension but denotes the top-left corner in the respective coordinate system
    // (origin in top-left corner of the display)
    private KVector offset;
    private double width;

    /**
     * Creates a beam-like extension.
     * 
     * @param parent
     *            Element this extension is going to be attached to
     * @param direction
     *            The cardinal direction this extension points to
     * @param middlePos
     *            The exact point the extension attaches to its element. If the extensions has a greater width than 1
     *            unit, the point cuts the beam exactly in half.
     * @param width
     *            Width of the extension
     */
    public DCExtension(final DCElement parent, final DCDirection direction, final KVector middlePos,
            final double width) {
        this.parent = parent;
        this.direction = direction;
        this.width = width;

        ElkRectangle bounds = parent.getBounds();
        setOffset(new KVector(-bounds.x, -bounds.y));
        getOffset().add(middlePos);
        // adjust position so that the Vector describes the left bound of the extension when viewed clock wise;
        double halfWidth = width / 2;
        // If horizontal extension ...
        if (direction.isHorizontal()) {
            // ... subtract half the width vertically
            getOffset().sub(0, halfWidth);
        } else {
            // other case is analog
            getOffset().sub(halfWidth, 0);
        }

        parent.addExtension(this);
    }

    /**
     * Returns the element this extension is attached to.
     * 
     * @return the parent Element this extension is attached to
     */
    public DCElement getParent() {
        return parent;
    }

    /**
     * Returns the cardinal direction of this extension.
     * 
     * @return the direction The cardinal direction of this extension.
     */
    public DCDirection getDirection() {
        return direction;
    }

    /**
     * Retruns the offset of the extension relative to the upper left corner of its element.
     * 
     * @return the offset offset relative to the element of the extension
     */
    public KVector getOffset() {
        return offset;
    }

    /**
     * Returns the width of the extensions.
     * 
     * @return the width
     */
    public double getWidth() {
        return width;
    }

    /**
     * Explicitly sets the offset of the extension relative to its element.
     * 
     * @param offset
     *            the offset to set
     */
    private void setOffset(final KVector offset) {
        this.offset = offset;
    }

}
