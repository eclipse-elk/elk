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

import java.util.List;
import java.util.Objects;

import org.eclipse.elk.alg.common.spore.InternalProperties;
import org.eclipse.elk.core.math.KVector;

import com.google.common.collect.Lists;
import com.google.common.math.DoubleMath;

/**
 *  Representation of a triangle.
 */
public class TTriangle {
    // CHECKSTYLEOFF VisibilityModifier
    /** vertex a. */
    public final KVector a;
    /** vertex b. */
    public final KVector b;
    /** vertex c. */
    public final KVector c;
    /** A list containing the three edges of the triangle. */
    public final List<TEdge> tEdges;
    /** A list containing the three vertices of the triangle. */
    public final List<KVector> vertices;
    /** he center of the circumscribed circle of the triangle. */
    private KVector circumcenter;
    // CHECKSTYLEON VisibilityModifier
    
    /**
     * Creates a triangle (a,b,c).
     * 
     * @param a
     * @param b
     * @param c
     */
    public TTriangle(final KVector a, final KVector b, final KVector c) {
        this.a = a;
        this.b = b;
        this.c = c;
        tEdges = Lists.newArrayList(new TEdge(a, b), new TEdge(b, c), new TEdge(c, a));
        vertices = Lists.newArrayList(a, b, c);
        circumcenter = calculateCircumcenter();
    }
    
    /**
     * The center of the circumscribed circle of the triangle.
     * 
     * @return the center point of the circumcircle
     */
    public KVector getCircumcenter() {
        return circumcenter;
    }
    
    /**
     * Determines the center of the circumscribed circle of the triangle.
     * 
     * @return the center point of the circumcircle
     */
    private KVector calculateCircumcenter() {
        KVector ab = b.clone().sub(a);
        KVector ac = c.clone().sub(a);
        KVector bc = c.clone().sub(b);
        double e = ab.x * (a.x + b.x) + ab.y * (a.y + b.y);
        double f = ac.x * (a.x + c.x) + ac.y * (a.y + c.y);
        double g = 2 * (ab.x * bc.y - ab.y * bc.x);
        
        double px = (ac.y * e - ab.y * f) / g;
        double py = (ab.x * f - ac.x * e) / g;
        return new KVector(px, py);
    }
    
    /**
     * Check if v is located inside the circumcircle of the triangle.
     * 
     * @param v point to check
     * @return true if v lies in circumcircle
     */
    public boolean inCircumcircle(final KVector v) {
        // Fuzziness prevents non-deterministic behavior caused by double imprecision e.g. let square a,b,c,d with d 
        // not in circumcircle a,b,c then b should not be in circumcircle a,c,d.
        return DoubleMath.fuzzyCompare(circumcenter.distance(v), circumcenter.distance(a), 
                InternalProperties.FUZZINESS) < 0;
    }
    
    /**
     * {@inheritDoc}
     * Two TTriangles are equal, if they are described by the same vertices.
     */
    @Override
    public boolean equals(final Object obj) {
        if (obj instanceof TTriangle) {
            TTriangle other = (TTriangle) obj;
            return contains(other.a) && contains(other.b) && contains(other.c);
        } else {
            return false;
        }
    }
    
    @Override
    public int hashCode() {
        // Note that the hash is independent of the vertex order to guarantee that equal triangles
        // yield identical hashes.
        return Objects.hashCode(a) + Objects.hashCode(b) + Objects.hashCode(c);
    }
    
    /**
     * Checks whether an edge is part of this triangle.
     * @param tEdge the edge to check
     * @return true if tEdge is equal to one of the triangle's edges and false otherwise
     */
    public boolean contains(final TEdge tEdge) {
        return Objects.equals(tEdge, tEdges.get(0))
            || Objects.equals(tEdge, tEdges.get(1))
            || Objects.equals(tEdge, tEdges.get(2));
    }
    
    /**
     * Checks whether a vertex is part of this triangle.
     * @param vertex the vertex to check
     * @return true if vertex is equal to one of the triangle's vertices and false otherwise
     */
    public boolean contains(final KVector vertex) {
        return Objects.equals(vertex, vertices.get(0))
                || Objects.equals(vertex, vertices.get(1))
                || Objects.equals(vertex, vertices.get(2));
    }
}
