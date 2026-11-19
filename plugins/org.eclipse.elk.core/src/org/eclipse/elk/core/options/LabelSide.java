/*******************************************************************************
 * Copyright (c) 2014, 2017 Kiel University and others.
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
package org.eclipse.elk.core.options;

import org.eclipse.elk.graph.properties.IProperty;
import org.eclipse.elk.graph.properties.Property;

/**
 * Enumeration for the definition of a side of the edge to place the (edge) label to.
 */
public enum LabelSide {
    /** The label's placement side hasn't been decided yet. */
    UNKNOWN,
    /** The label is placed above the edge. */
    ABOVE,
    /** The label is placed below the edge. */
    BELOW,
    /** The label is placed directly on top of the edge. */
    INLINE;
    

    /**
     * Property set on edge and port labels by layout algorithms depending on which side they decide is
     * appropriate for any given label.
     */
    public static final IProperty<LabelSide> LABEL_SIDE = new Property<LabelSide>(
            "org.eclipse.elk.labelSide", LabelSide.UNKNOWN);
    
    
    /**
     * Returns the side opposite to the one this method is called on. {@link #UNKNOWN} is mapped to itself.
     */
    public LabelSide opposite() {
        switch (this) {
        case ABOVE:
            return BELOW;
        case BELOW:
            return ABOVE;
        case INLINE:
            return INLINE;
        default:
            return UNKNOWN;
        }
    }
}
