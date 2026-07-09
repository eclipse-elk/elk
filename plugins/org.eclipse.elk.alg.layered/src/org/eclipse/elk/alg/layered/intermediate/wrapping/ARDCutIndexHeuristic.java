/*******************************************************************************
 * Copyright (c) 2016 Kiel University and others.
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
package org.eclipse.elk.alg.layered.intermediate.wrapping;

import java.util.List;

import org.eclipse.elk.alg.layered.graph.LGraph;

import com.google.common.collect.Lists;

/**
 * Calculates cut indexes based on a desired aspect ratio. Given a layering and estimations on the width and height of
 * each layer, first the number of desired 'chunks' is calculated, which is then used to select equidistant cut indexes.
 */
public class ARDCutIndexHeuristic implements ICutIndexCalculator {

    @Override
    public List<Integer> getCutIndexes(final LGraph graph, final GraphStats gs) {
        int rows = getChunkCount(gs);

        // the number of cuts is one less than the number of rows
        List<Integer> cuts = Lists.newArrayList();
        double step = gs.longestPath / (double) rows;
        for (int idx = 1; idx < rows; ++idx) {
            cuts.add((int) Math.round(idx * step));
        }

        return cuts;
    }

    /**
     * Given the estimations on width and height in {@code gs}, the method seeks a number of 'chunks' s.t.
     * 
     * <pre>
     * desiredAr = ((longestPath / chunks) * maxWidth) / (chunks * maxHeight)
     * </pre>
     * 
     * @param gs
     *            statistics about the graph.
     * @return the number of chunks the graph should be split into. The resulting number of cuts should therefore be one
     *         less.
     */
    public static int getChunkCount(final GraphStats gs) {
        double rowsd = Math.sqrt(gs.getSumWidth() / (gs.dar * gs.getMaxHeight()));
        int rows = (int) Math.round(rowsd);
        rows = Math.min(rows, gs.longestPath);
        return rows;
    }

    @Override
    public boolean guaranteeValid() {
        return false;
    }
    
}
