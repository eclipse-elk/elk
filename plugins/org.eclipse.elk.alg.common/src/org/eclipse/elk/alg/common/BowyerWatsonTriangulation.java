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

import java.util.Iterator;
import java.util.List;
import java.util.Set;

import org.eclipse.elk.alg.common.utils.SVGImage;
import org.eclipse.elk.core.math.KVector;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;

/**
 * This class creates a Delaunay triangulation for the given points that is represented by a list of edges.
 * The approach was devised independently by Bowyer, Adrian (1981) "Computing Dirichlet tessellations" and 
 * Watson, David F. (1981) "Computing the n-dimensional Delaunay tessellation with application to Voronoi polytopes".
 * <p> precondition: All vertices have to be distinct.</p>
 * <p> postcondition: The returned edges form a connected graph that includes all input points without duplicate 
 * edges.</p>
 */
public final class BowyerWatsonTriangulation {
    
    /** Hidden constructor. */
    private BowyerWatsonTriangulation() { };
    
    /**
     * Triangulates a list of points.
     * 
     * @param vertices the input points
     * @param debugOutputFile file name for debug SVG. Debug output will be deactivated if this is null.
     * @return the edges of the triangulation
     */
    public static Set<TEdge> triangulate(final List<KVector> vertices, final String debugOutputFile) {
        /*d*/SVGImage svg = new SVGImage(debugOutputFile);
        /*d*/svg.addGroups("invalid", "tri", "bndry", "done", "new");
        
        /* preliminaries */
        
        // determine the bounding box of the given points
        KVector topleft = new KVector(Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY);
        KVector bottomright = new KVector(Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY);
        for (KVector v : vertices) {
            topleft.x = Math.min(topleft.x, v.x);
            topleft.y = Math.min(topleft.y, v.y);
            bottomright.x = Math.max(bottomright.x, v.x);
            bottomright.y = Math.max(bottomright.y, v.y);
            // CHECKSTYLEOFF MagicNumber
            /*d*/svg.g("bb").addCircle(v.x, v.y, 18, "stroke=\"black\" stroke-width=\"1\" fill=\"lightgray\"");
        }
        KVector size = new KVector(bottomright.x - topleft.x, bottomright.y - topleft.y);
        /*d*/svg.g("bb").addRect(topleft.x, topleft.y, size.x, size.y, 
                "stroke=\"blue\" stroke-width=\"4\" fill=\"none\"");
        
        // find a super-triangle spanning over all vertices
        final double wiggleroom = 50;
        // This ensures that no vertices are too close to the edges of the super-triangle
        // because that would result in circumcircles that are so large that a double would not
        // be precise enough to compute the circumcircle criterion.
        KVector sa = new KVector(topleft.x - wiggleroom, topleft.y - size.x - wiggleroom);
        KVector sb = new KVector(topleft.x - wiggleroom, bottomright.y + size.x + wiggleroom);
        KVector sc = new KVector(bottomright.x + size.y / 2 + wiggleroom, topleft.y + size.y / 2);
        /*d*/svg.g("bb").addPoly("stroke=\"gray\" stroke-width=\"4\" fill=\"none\" stroke-dasharray=\"20,20\"", 
                sa, sb, sc, sa);
        TTriangle superTriangle = new TTriangle(sa, sb, sc);
        /*d*/svg.setViewBox(sa.x, sa.y, sc.x - sa.x, sb.y - sa.y); // circumcircles will make it gigantic
        /*d*/svg.isave();
        /*d*/svg.removeGroup("bb");
        
        /* Bowyer Watson algorithm*/
        
        Set<TTriangle> triangulation = Sets.newHashSet();
        List<TTriangle> invalidTriangles = Lists.newArrayList();
        List<TEdge> boundary = Lists.newArrayList();
        triangulation.add(superTriangle);
        
        // incrementally adding vertices
        for (KVector vertex : vertices) {
            /*d*/svg.g("done").addCircle(vertex.x, vertex.y, 18, 
                    "stroke=\"black\" stroke-width=\"1\" fill=\"lightgray\"");
            /*d*/svg.g("new").addCircle(vertex.x, vertex.y, 18, "stroke=\"black\" stroke-width=\"1\" fill=\"black\"");
            
            // gather invalid triangles where the new vertex lies inside the circumcircle
            invalidTriangles.clear();
            for (TTriangle triangle : triangulation) {
                /*d*/svg.g("tri").addPoly("stroke=\"black\" fill=\"none\" stroke-width=\"4\"", triangle.a, 
                        triangle.b, triangle.c, triangle.a);
                /*d*/KVector c = triangle.getCircumcenter();
                /*d*/svg.g("invalid").addCircle(c.x, c.y, c.distance(triangle.a), 
                        "stroke=\"orange\" stroke-width=\"4\" fill=\"none\"");
                if (triangle.inCircumcircle(vertex)) {
                    invalidTriangles.add(triangle);
                    /*d*/svg.g("invalid").addPoly("stroke=\"none\" fill=\"red\" opacity=\"0.18\"", triangle.a, 
                            triangle.b, triangle.c, triangle.a);
                }
            }
            /*d*/svg.isave();
            /*d*/svg.clearGroup("invalid");
            
            // calculate boundary of invalid triangles
            boundary.clear();
            for (TTriangle triangle : invalidTriangles) {
                for (TEdge tEdge : triangle.tEdges) {
                    boolean onBoundary = true;
                    // edges that are not shared with other invalid triangles are on the boundary
                    for (TTriangle other : invalidTriangles) {
                        if (other != triangle && other.contains(tEdge)) {
                            onBoundary = false;
                        }
                    }
                    if (onBoundary) {
                        boundary.add(tEdge);
                        /*d*/svg.g("bndry").addLine(tEdge.u.x, tEdge.u.y, tEdge.v.x, tEdge.v.y, 
                                "stroke=\"purple\" stroke-width=\"18\" stroke-dasharray=\"20,20\"");
                    }
                }
            }
            /*d*/svg.isave();
            
            // remove invalid triangles
            triangulation.removeAll(invalidTriangles);
            /*d*/svg.clearGroup("tri");
            /*d*/triangulation.forEach(triangle -> 
                        svg.g("tri").addPoly("stroke=\"black\" fill=\"none\" stroke-width=\"4\"", triangle.a, 
                                triangle.b, triangle.c, triangle.a));
            /*d*/svg.isave();
            
            // triangulate boundary
            for (TEdge tEdge : boundary) {
                triangulation.add(new TTriangle(vertex, tEdge.u, tEdge.v));
                /*d*/svg.g("tri").addPoly("stroke=\"black\" fill=\"none\" stroke-width=\"4\"", vertex, 
                        tEdge.u, tEdge.v, vertex);
            }
            /*d*/svg.isave();
            /*d*/svg.clearGroup("new");
            /*d*/svg.clearGroup("bndry");
            /*d*/svg.clearGroup("tri");
        }
        
        // convert triangulation to set of edges
        // Because it's a set, edges that were redundant in the triangulation exist only once in the set.
        Set<TEdge> tEdges = Sets.newHashSet();
        triangulation.forEach(triangle -> tEdges.addAll(triangle.tEdges));
        
        // remove edges connected to super triangle
        Iterator<TEdge> i = tEdges.iterator();
        while (i.hasNext()) {
            TEdge tEdge = i.next();
            if (superTriangle.contains(tEdge.u) || superTriangle.contains(tEdge.v)) {
                i.remove();
            }
        }
        
        /*d*/tEdges.forEach(tEdge -> 
                svg.addLine(tEdge.u.x, tEdge.u.y, tEdge.v.x, tEdge.v.y, "stroke=\"black\" stroke-width=\"4\""));
        /*d*/svg.isave();
        
        return tEdges;
    }
    
    /**
     * Triangulates a list of points.
     * 
     * @param vertices the input points
     * @return the edges of the triangulation
     */
    public static Set<TEdge> triangulate(final List<KVector> vertices) {
        return triangulate(vertices, null);
    }
}
