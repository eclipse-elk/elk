/*******************************************************************************
 * Copyright (c) 2013, 2015 Kiel University and others.
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
package org.eclipse.elk.alg.mrtree.graph;

import org.eclipse.elk.graph.properties.MapPropertyHolder;

/**
 * Abstract superclass for the layers, nodes, ports, and edges of a tree graph.
 * 
 * @author sor
 * @author sgu
 */
public abstract class TGraphElement extends MapPropertyHolder {

    /** the serial version UID. */
    private static final long serialVersionUID = 1L;

    /** Identifier value, may be arbitrarily used by algorithms. */
    public int id;  // SUPPRESS CHECKSTYLE VisibilityModifier

    /**
     * Create a graph element.
     * 
     * @param id the element identifier
     */
    public TGraphElement(final int id) {
        super();
        this.id = id;
    }

    /**
     * Implicit super constructor.
     * 
     */
    public TGraphElement() {
    }

}
