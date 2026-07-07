/*******************************************************************************
 * Copyright (c) 2016 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.layered.intermediate.wrapping;

import java.util.Collections;
import java.util.List;

import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.layered.options.WrappingStrategy;

/**
 * Calculate cut indexes that are used to cut a graph's layering into multiple chunks and 'wrap' the drawing to improve
 * the utilization of a given drawing area.
 * 
 * @see WrappingStrategy
 */
public interface ICutIndexCalculator {

    /**
     * @param graph
     *            a graph.
     * @param gs
     *            some further precalculated information about the graph.
     * @return a list of indexes at which the graph's layering should be split.
     */
    List<Integer> getCutIndexes(LGraph graph, GraphStats gs);

    /**
     * @return whether it is guaranteed that all cuts are valid.
     */
    boolean guaranteeValid();
    
    /**
     * Simple {@link ICutIndexCalculator} that reads the manually specified cut indexes 
     * from the layout option {@link LayeredOptions#WRAPPING_CUTTING_CUTS}.
     */
    public static class ManualCutIndexCalculator implements ICutIndexCalculator {
        @Override
        public List<Integer> getCutIndexes(final LGraph graph, final GraphStats gs) {
            List<Integer> cuts = graph.getProperty(LayeredOptions.WRAPPING_CUTTING_CUTS);
            if (cuts != null) {
                return cuts;
            } else {
                return Collections.emptyList();
            }
        }
        @Override
        public boolean guaranteeValid() {
            return false;
        }
    }
    
}
