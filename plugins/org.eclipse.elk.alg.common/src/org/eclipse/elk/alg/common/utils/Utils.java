/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.common.utils;

import java.util.List;

import org.eclipse.elk.alg.common.TEdge;
import org.eclipse.elk.core.math.ElkRectangle;

import com.google.common.collect.Lists;

/**
 * Utilities for geometric calculations.
 */
public final class Utils {
    
    /** Hidden default constructor. */
    private Utils() { }

    /**
     * Overlap of two axis aligned rectangles, i.e. the factor by that the line
     * connecting the center points of the rectangles has to be stretched so that the rectangles just touch.
     * The center points have to be distinct.
     * @param r1
     * @param r2
     * @return the overlap
     */
    public static double overlap(final ElkRectangle r1, final ElkRectangle r2) {
        double horizontalOverlap = Math.min(Math.abs(r1.x - (r2.x + r2.width)),
                                            Math.abs(r1.x + r1.width - r2.x));
        double verticalOverlap = Math.min(Math.abs(r1.y - (r2.y + r2.height)),
                                          Math.abs(r1.y + r1.height - r2.y));
        double horizontalCenterDistance = Math.abs((r1.x + r1.width / 2.0) - (r2.x + r2.width / 2.0));
        if (horizontalCenterDistance > r1.width / 2.0 + r2.width / 2.0) {
            return 1.0;
        }
        double verticalCenterDistance = Math.abs((r1.y + r1.height / 2.0) - (r2.y + r2.height / 2.0));
        if (verticalCenterDistance > r1.height / 2.0 + r2.height / 2.0) {
            return 1.0;
        }
        if (horizontalCenterDistance == 0.0 && verticalCenterDistance == 0.0) {
            return 0;
        }
        if (horizontalCenterDistance == 0.0) {
            return verticalOverlap / verticalCenterDistance + 1;
        }
        if (verticalCenterDistance == 0.0) {
            return horizontalOverlap / horizontalCenterDistance + 1;
        }
        return Math.min(horizontalOverlap / horizontalCenterDistance, verticalOverlap / verticalCenterDistance) + 1;
    }
    
    /**
     * returns a list of the four edges of a rectangle.
     * @param r the rectangle
     * @return the edges
     */
    public static List<TEdge> getRectEdges(final ElkRectangle r) {
        List<TEdge> rectEdegs = Lists.newArrayList();
        rectEdegs.add(new TEdge(r.getTopLeft(), r.getTopRight()));
        rectEdegs.add(new TEdge(r.getTopLeft(), r.getBottomLeft()));
        rectEdegs.add(new TEdge(r.getBottomRight(), r.getTopRight()));
        rectEdegs.add(new TEdge(r.getBottomRight(), r.getBottomLeft()));
        return rectEdegs;
    }
}
