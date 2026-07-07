/*******************************************************************************
 * Copyright (c) 2013, 2015 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.mrtree.graph;

import org.eclipse.elk.core.math.KVector;

/**
 * Abstract superclass for {@link TGraphElement}s that can have a position and a size.
 * 
 * @author sor
 * @author sgu
 */
public abstract class TShape extends TGraphElement {

    /** the serial version UID. */
    private static final long serialVersionUID = 1L;

    /** the current position of the element. */
    private final KVector pos = new KVector();
    /** the size of the element. */
    private final KVector size = new KVector();

    /**
     * Creates a shape in the context of the given graph.
     * 
     * @param id the identifier
     */
    public TShape(final int id) {
        super(id);
    }

    /**
     * Implicit super constructor.
     */
    public TShape() {
    }

    /**
     * Returns the element's current position. For nodes this is the coordinates of their center
     * point, for other elements the coordinates of the upper left corner.
     * 
     * @return the position
     */
    public KVector getPosition() {
        return pos;
    }

    /**
     * Returns the element's current size.
     * 
     * @return the size
     */
    public KVector getSize() {
        return size;
    }
}
