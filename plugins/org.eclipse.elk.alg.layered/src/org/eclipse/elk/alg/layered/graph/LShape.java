/*******************************************************************************
 * Copyright (c) 2011, 2019 Kiel University and others.
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
package org.eclipse.elk.alg.layered.graph;

import org.eclipse.elk.core.math.KVector;

/**
 * Abstract superclass for {@link LGraphElement}s that can have a position and a size.
 */
public abstract class LShape extends LGraphElement {
    
    /** the serial version UID. */
    private static final long serialVersionUID = 5111245934175354687L;
    
    /** the current position of the element. */
    private final KVector pos = new KVector();
    /** the size of the element. */
    private final KVector size = new KVector();
    

    /**
     * Returns the element's current position. This is the coordinate of the element's upper
     * left corner. May be modified.
     * 
     * @return the position
     */
    public KVector getPosition() {
        return pos;
    }

    /**
     * Returns the element's current size. May be modified.
     * 
     * @return the size
     */
    public KVector getSize() {
        return size;
    }
    
}
