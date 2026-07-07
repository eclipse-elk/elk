/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.common;

import java.util.Objects;

import org.eclipse.elk.core.math.KVector;

/**
 * An undirected edge modeled as two vectors.
 */
public class TEdge {
    // CHECKSTYLEOFF VisibilityModifier
    /** vertex u. */
    public final KVector u;
    /** vertex v. */
    public final KVector v;
    // CHECKSTYLEON VisibilityModifier
    
    /**
     * Creates an edge (u,v).
     * 
     * @param u
     * @param v
     */
    public TEdge(final KVector u, final KVector v) {
        this.u = u;
        this.v = v;
    }
    
    @Override
    public boolean equals(final Object obj) {
        if (obj instanceof TEdge) {
            TEdge other = (TEdge) obj;
            return (Objects.equals(u, other.u) && Objects.equals(v, other.v)) 
                || (Objects.equals(u, other.v) && Objects.equals(v, other.u));
        } else {
            return false;
        }
    }
    
    @Override
    public int hashCode() {
        return Objects.hashCode(u) + Objects.hashCode(v);
    }
}
