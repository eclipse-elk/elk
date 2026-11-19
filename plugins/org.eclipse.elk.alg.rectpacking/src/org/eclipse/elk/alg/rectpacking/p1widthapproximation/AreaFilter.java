/*******************************************************************************
 * Copyright (c) 2018, 2020 Kiel University and others.
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
package org.eclipse.elk.alg.rectpacking.p1widthapproximation;

import java.util.ArrayList;
import java.util.List;

import org.eclipse.elk.alg.rectpacking.util.DrawingData;
import org.eclipse.elk.core.math.ElkPadding;

/**
 * This class implements a concrete strategy in the strategy pattern given by {@link BestCandidateFilter}.
 * <p>
 * This class offers a method that filters the given list regarding the area. The candidate or candidates with equal
 * area with the smallest area are returned in the filtered list.
 * </p>
 */
public class AreaFilter implements BestCandidateFilter {

    @Override
    public List<DrawingData> filterList(final List<DrawingData> candidates, final double aspectRatio,
            final ElkPadding padding) {
        List<DrawingData> remainingCandidates = new ArrayList<DrawingData>();
        double minArea = Double.POSITIVE_INFINITY;
        for (DrawingData opt : candidates) {
            minArea = Math.min(minArea, (opt.getDrawingWidth() + padding.getHorizontal())
                    * (opt.getDrawingHeight() + padding.getVertical()));
        }
        for (DrawingData candidate : candidates) {
            if ((candidate.getDrawingWidth() + padding.getHorizontal())
                    * (candidate.getDrawingHeight() + padding.getVertical()) == minArea) {
                remainingCandidates.add(candidate);
            }
        }
        return remainingCandidates;
    }
}
