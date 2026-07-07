/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.disco.graph;

import java.util.Collection;
import java.util.Set;

import org.eclipse.elk.alg.disco.ICompactor;
import org.eclipse.elk.core.math.ElkRectangle;
import org.eclipse.elk.core.math.KVector;

import com.google.common.collect.Sets;

/**
 * Represents a connected component of a {@link DCGraph graph}.
 */
public class DCComponent {

    ///////////////////////////////////////////////////////////////////////////////
    // Variables

    /** Offset from original position, might be set by a class implementing {@link ICompactor}. */
    private KVector offset;
    /** All elements that are part of this component (e.g. nodes, edges, ...) */
    private Set<DCElement> shapes;
    /** Indicates whether some fields need to be updated. */
    private boolean changed = true;
    /** dimensions of the bounding rectangle of this component. */
    private KVector bounds;
    /** Upper left corner of the bounding rectangle of this component. */
    private KVector minCornerOfBoundingRectangle;
    /** For debugging purposes only. */
    private int id = -1;

    ///////////////////////////////////////////////////////////////////////////////
    // Package private constructor

    /**
     * Creates a new component with an offset of (0,0) and and no associated {@link DCElement DCElements}.
     */
    DCComponent() {
        offset = new KVector(0, 0);
        shapes = Sets.newHashSet();
    }

    ///////////////////////////////////////////////////////////////////////////////
    // Public methods

    /**
     * Sets the offset of the {@link DCComponent}.
     * 
     * @param offset
     *            New offset to add to the position of this {@link DCComponent}
     */
    public void setOffset(final KVector offset) {
        changed = true;
        this.offset = offset;
    }

    /**
     * Gets the offset of the {@link DCComponent}.
     * 
     * @return Offset to add to the position of this {@link DCComponent}
     */
    public KVector getOffset() {
        return offset;
    }

    /**
     * Gets all {@link DCElement DCElements} associated with this {@link DCComponent}.
     * 
     * @return All {@link DCElement DCElements} of this {@link DCComponent}
     */
    public Set<DCElement> getElements() {
        return shapes;
    }

    /**
     * Gets dimensions of the bounding rectangle of this {@link DCComponent}.
     * 
     * @return The bounding rectangle of this component
     */
    public KVector getDimensionsOfBoundingRectangle() {
        if (changed) {
            update();
        }
        return bounds;
    }

    /**
     * Gets the position of the upper left corner of this {@link DCComponent}.
     * 
     * @return The position of the upper left corner of this component
     */
    public KVector getMinCorner() {
        if (changed) {
            update();
        }
        return minCornerOfBoundingRectangle;
    }

    /**
     * Tests whether this {@link DCComponent} intersects with a rectangular area given by the four parameters of this
     * method.
     * 
     * @param rect
     * @return true: rectangle intersects with this {@link DCComponent}; false: otherwise.
     */
    public boolean intersects(final ElkRectangle rect) {
        for (DCElement elem : shapes) {
            if (elem.intersects(rect)) {
                return true;
            }
        }
        return false;
    }

    ///////////////////////////////////////////////////////////////////////////////
    // Package private methods for constructing the component

    /**
     * Adds a {@link DCElement} to this {@link DCComponent}.
     * 
     * @param element
     *            The element to add to this {@link DCComponent}
     */
    void addElement(final DCElement element) {
        changed = true;
        shapes.add(element);
        element.setComponent(this);
    }

    /**
     * Adds {@link DCElement DCElements} to this {@link DCComponent}.
     * 
     * @param elements
     *            The elements to add to this {@link DCComponent}
     */
    <C extends Collection<DCElement>> void addElements(final C elements) {
        for (DCElement elem : elements) {
            addElement(elem);
        }
    }

    ///////////////////////////////////////////////////////////////////////////////
    // Private method for updating internal state

    /**
     * Updates the position and the dimensions of the bounding rectangle of this {@link DCComponent}. To be called
     * whenever the component has been modified and properties of the bounding rectangle are being requested.
     */
    private void update() {
        changed = false;
        double minX = Double.POSITIVE_INFINITY;
        double maxX = Double.NEGATIVE_INFINITY;
        double minY = Double.POSITIVE_INFINITY;
        double maxY = Double.NEGATIVE_INFINITY;

        for (DCElement elem : shapes) {
            ElkRectangle elemBounds = elem.getBounds();

            minX = Math.min(minX, elemBounds.x);
            maxX = Math.max(maxX, elemBounds.x + elemBounds.width);
            minY = Math.min(minY, elemBounds.y);
            maxY = Math.max(maxY, elemBounds.y + elemBounds.height);

            for (DCExtension ext : elem.getExtensions()) {
                DCDirection dir = ext.getDirection();
                double minPos, maxPos;

                if (dir.isHorizontal()) {
                    minPos = elemBounds.y + ext.getOffset().y;
                    maxPos = minPos + ext.getWidth();
                    minY = Math.min(minY, minPos);
                    maxY = Math.max(maxY, maxPos);
                } else {
                    minPos = elemBounds.x + ext.getOffset().x;
                    maxPos = minPos + ext.getWidth();
                    minX = Math.min(minX, minPos);
                    maxX = Math.max(maxX, maxPos);
                }
            }
        }

        bounds = new KVector(maxX - minX, maxY - minY);
        minCornerOfBoundingRectangle = new KVector(minX + offset.x, minY + offset.y);
    }

    ///////////////////////////////////////////////////////////////////////////////
    // For debugging purposes only

    /**
     * Gets the integer id of this {@link DCComponent}.
     * 
     * @return the id
     */
    public int getId() {
        return id;
    }

    /**
     * Sets the integer id of this {@link DCComponent}.
     * 
     * @param id
     *            the id to set
     */
    public void setId(final int id) {
        this.id = id;
    }
}
