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

import java.util.Iterator;
import java.util.List;

import org.eclipse.elk.alg.disco.ICompactor;
import org.eclipse.elk.alg.disco.transform.IGraphTransformer;
import org.eclipse.elk.core.math.ElkMath;
import org.eclipse.elk.core.math.ElkRectangle;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.math.KVectorChain;
import org.eclipse.elk.graph.properties.MapPropertyHolder;

import com.google.common.collect.Lists;

/**
 * Represents an element of a {@link DCGraph}. Expresses common graph features like nodes, edges, labels, ... simply as
 * polygons, as this abstraction is sufficient for the compaction of different connected components, whereas short
 * hierarchical edges are represented by {@link DCExtension DCExtensions}.
 */
public class DCElement extends MapPropertyHolder {
    /** the serial version UID. */
    private static final long serialVersionUID = 58841189716563145L;

    ///////////////////////////////////////////////////////////////////////////////
    // Variables

    /** Path representing this elements shape based on the positions of polygon vertices . */
    private final KVectorChain shape;
    /** The bounding box of the {@link #shape}. */
    private final ElkRectangle bounds;
    /** The {@link DCComponent} this element belongs to. */
    private DCComponent cp;
    /** Maintained for drawing purposes in the DisCoGraphRenderer of the DisCoDebugView. */
    private double[] coords;
    /** Maintained for drawing purposes in the DisCoGraphRenderer of the DisCoDebugView. */
    private KVector parentCoords = null;
    /** List holding all short hierarchical edges to or from the parent node (if this DCElement represents a node). */
    private List<DCExtension> extensions = Lists.newArrayList();

    ///////////////////////////////////////////////////////////////////////////////
    // Constructor

    /**
     * Constructs a {@link DCElement} from a list of of points representing the vertices of a closed polygonal path.
     * Typically invoked by a class implementing {@link IGraphTransformer}.
     * 
     * @param polyPath
     *            Ordered list of points expressing the vertices of a Polygon
     */
    public DCElement(final KVectorChain polyPath) {
        this.shape = polyPath;
        double minX = Double.POSITIVE_INFINITY;
        double minY = Double.POSITIVE_INFINITY;
        double maxX = Double.NEGATIVE_INFINITY;
        double maxY = Double.NEGATIVE_INFINITY;
        for (KVector v : polyPath) {
            minX = Math.min(minX, v.x);
            minY = Math.min(minY, v.y);
            maxX = Math.max(maxX, v.x);
            maxY = Math.max(maxY, v.y);
        }
        bounds = new ElkRectangle(minX, minY, maxX - minX, maxY - minY);
    }

    ///////////////////////////////////////////////////////////////////////////////
    // Public methods

    /**
     * Returns the offset which has to be added to the elements original position to reach its new position after
     * component compactions using a class implementing {@link ICompactor}.
     * 
     * @return Offset
     * @throws NullPointerException
     *             Offset hasn't been set, yet.
     */
    public KVector getOffset() throws NullPointerException {
        return cp.getOffset();
    }

    /**
     * Used for debugging purposes in the DisCoDebugView. Returns the coordinates of the vertices of this element in an
     * easily paintable form.
     *
     * @return coordinates in the form {@code (x_0, y_0, x_1, y_1, ..., x_n, y_n)} for points {@code p_i = (x_i, y_i)}
     *         representing a polygon (a closed path) p_0,p_1,...,p_n,p_0
     */
    public double[] getCoords() {
        if (coords != null) {
            return coords;
        }

        coords = new double[shape.size() * 2];
        Iterator<KVector> it = shape.iterator();
        int i = 0;
        while (it.hasNext()) {
            KVector c = it.next();
            coords[i] = c.x;
            coords[i + 1] = c.y;
            i += 2;
        }

        return coords;
    }

    /**
     * Adds an {@link DCExtension} to this element.
     * 
     * @param extension
     *            {@link DCExtension} to add.
     */
    public void addExtension(final DCExtension extension) {
        extensions.add(extension);
    }

    /**
     * Returns the extensions of this element.
     * 
     * @return list of extensions
     */
    public List<DCExtension> getExtensions() {
        return extensions;
    }

    /**
     * Gets the bounding rectangle of this element.
     * 
     * @return The bounding rectangle of this element
     */
    public ElkRectangle getBounds() {
        return bounds;
    }

    /**
     * Sets the coordinates of the surrounding parent of this element, Needed for the debug view, only.
     * 
     * @param pCoords
     *            Coordinates of the parent of the element
     */
    public void setParentCoords(final KVector pCoords) {
        parentCoords = pCoords;
    }

    /**
     * Returns the coordinates of the parent of this element for debugging purposes.
     * 
     * @return The coordinates of the surrounding parent of the element. {@link <code>null</code>} if not set.
     */
    public KVector getParentcoords() {
        return parentCoords;
    }

    ///////////////////////////////////////////////////////////////////////////////
    // Package private methods

    /**
     * Sets the {@link DCComponent} this element belongs to.
     * 
     * @param comp
     *            The {@link DCComponent} this element belongs to.
     */
    void setComponent(final DCComponent comp) {
        this.cp = comp;
    }

    /**
     * Tests whether this {@link DCElement} intersects with, or fully contains, 
     * a rectangular area given by the four parameters.
     * 
     * @param rect
     * @return true: rectangle intersects with this {@link DCElement}; false: otherwise.
     */
    boolean intersects(final ElkRectangle rect) {
        return ElkMath.intersects(rect, shape) || ElkMath.contains(rect, shape); 
    }
}
