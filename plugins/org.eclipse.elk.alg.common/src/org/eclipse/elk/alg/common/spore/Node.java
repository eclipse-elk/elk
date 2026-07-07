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
package org.eclipse.elk.alg.common.spore;

import org.eclipse.elk.alg.common.TEdge;
import org.eclipse.elk.alg.common.utils.Utils;
import org.eclipse.elk.core.math.ElkMath;
import org.eclipse.elk.core.math.ElkRectangle;
import org.eclipse.elk.core.math.KVector;

import com.google.common.math.DoubleMath;

/**
 * Representation of a node as a rectangle and center point.
 */
public class Node {
    // CHECKSTYLEOFF VisibilityModifier
    /** The original center point from the time of creation. Used to match the object this node was derived from. */
    public KVector originalVertex;
    /** A modifiable center position to that a translation can be applied. */
    public KVector vertex;
    /** The bounding box of the represented diagram element. */
    public ElkRectangle rect;
    // CHECKSTYLEON VisibilityModifier
    
    /**
     * Constructor initializing the rectangle and center point.
     * @param v center point, also used for identification
     * @param r rectangle
     */
    public Node(final KVector v, final ElkRectangle r) {
        originalVertex = v;
        vertex = originalVertex.clone();
        rect = new ElkRectangle(r);
    }
    
    /**
     * Translate the rectangle and modifiable center position by vector v.
     * @param v the translation vector
     */
    public void translate(final KVector v) {
        vertex.add(v);
        rect.x += v.x;
        rect.y += v.y;
    }
    
    /**
     * Helper to keep the center position and rectangle consistent.
     * @param p new center position
     */
    public void setCenterPosition(final KVector p) {
        translate(p.clone().sub(vertex));
    }
    
    /**
     * Returns the underlap between two nodes, i.e. how far the nodes can be pushed together on the line connecting
     * their centers without collision.
     * 
     * 
     * @param other the other node
     * @return the underlap
     */
    public double underlap(final Node other) {
        
        double horizontalCenterDistance = Math.abs(rect.getCenter().x - other.rect.getCenter().x);
        double verticalCenterDistance = Math.abs(rect.getCenter().y - other.rect.getCenter().y);
        double horizontalUnderlap = 0;
        double verticalUnderlap = 0;
        double hScale = 1;
        double vScale = 1;
        if (horizontalCenterDistance > rect.width / 2.0 + other.rect.width / 2.0) {
            horizontalUnderlap = Math.min(Math.abs(rect.x - (other.rect.x + other.rect.width)),
                    Math.abs(rect.x + rect.width - other.rect.x));
            hScale = 1 - horizontalUnderlap / horizontalCenterDistance;
        }
        if (verticalCenterDistance > rect.height / 2.0 + other.rect.height / 2.0) {
            verticalUnderlap = Math.min(Math.abs(rect.y - (other.rect.y + other.rect.height)),
                  Math.abs(rect.y + rect.height - other.rect.y));
            vScale = 1 - verticalUnderlap / verticalCenterDistance;
        }
        double scale = Math.min(hScale, vScale);
        return (1 - scale) * Math.sqrt(horizontalCenterDistance * horizontalCenterDistance
                                     + verticalCenterDistance * verticalCenterDistance);
    }
    
    /**
     * Distance between two nodes, i.e. how far the other node can be moved in the direction of v without colliding
     * with this node.
     * @param other the other node
     * @param v the direction
     * @return the direction dependent distance
     */
    public double distance(final Node other, final KVector v) {
        double result = Double.POSITIVE_INFINITY;
        for (TEdge e1 : Utils.getRectEdges(this.rect)) {
            for (TEdge e2 : Utils.getRectEdges(other.rect)) {
                double distance = ElkMath.distance(e1.u, e1.v, e2.u, e2.v, v);
                result = Math.min(result, distance);
            }
        }
        return result;
    }
    
    /**
     * Determines if two nodes touch.
     * @param other the other node
     * @return true if they touch
     */
    public boolean touches(final Node other) {
        return DoubleMath.fuzzyCompare(this.rect.x, other.rect.x + other.rect.width, InternalProperties.FUZZINESS) <= 0
            && DoubleMath.fuzzyCompare(other.rect.x, this.rect.x + this.rect.width, InternalProperties.FUZZINESS) <= 0
            && DoubleMath.fuzzyCompare(this.rect.y, other.rect.y + other.rect.height, InternalProperties.FUZZINESS) <= 0
            && DoubleMath.fuzzyCompare(other.rect.y, this.rect.y + this.rect.height, InternalProperties.FUZZINESS) <= 0;
    }
}
