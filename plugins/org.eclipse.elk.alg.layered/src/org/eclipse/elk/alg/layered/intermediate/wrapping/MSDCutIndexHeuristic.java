/*******************************************************************************
 * Copyright (c) 2016 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.alg.layered.intermediate.wrapping;

import java.util.List;
import java.util.ArrayList;

import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.options.LayeredOptions;


/**
 * As opposed to the simple {@link ARDCutIndexHeuristic}, this heuristic tries to directly optimize the max scale value
 * of the resulting drawing. The underlying idea is to find, for a fixed number of cut indexes, those indexes that
 * minimize the 'width part' of max scale.
 * 
 * Additionally, it allows to slightly alter the number of cut indexes via the
 * {@link LayeredOptions#WRAPPING_CUTTING_MSD_FREEDOM} layout option at the expense of more iterations.
 */
public class MSDCutIndexHeuristic implements ICutIndexCalculator {

    @Override
    public boolean guaranteeValid() {
        return false;
    }
    
    @Override
    public List<Integer> getCutIndexes(final LGraph graph, final GraphStats gs) {
        
        // minimize the max of the sums of the widths
        double[] widths = gs.getWidths();
        double[] heights = gs.getHeights();
        
        // record the accumulated width at each index of the original layering
        double[] widthAtIndex = new double[widths.length];
        widthAtIndex[0] = widths[0];
        double total = widths[0];
        for (int i = 1; i < widths.length; i++) {
            widthAtIndex[i] = widthAtIndex[i - 1] + widths[i];
            total += widths[i];
        }
        
        // initial guess on a good number of cuts
        int cutCnt = ARDCutIndexHeuristic.getChunkCount(gs) - 1;
        int freedom = graph.getProperty(LayeredOptions.WRAPPING_CUTTING_MSD_FREEDOM);

        double bestMaxScale = Double.NEGATIVE_INFINITY;
        List<Integer> bestCuts = new ArrayList<>();
        
        // now find the best set of cut indexes
        for (int m = Math.max(0, cutCnt - freedom); 
                 m <= Math.min(gs.longestPath - 1, cutCnt + freedom); 
                 m++) {
            
            // calculate cuts
            double rowSum = total / (m + 1);
            double sumSoFar = 0;
            int index = 1;
            List<Integer> cuts = new ArrayList<>();
            
            // maximum of any row width
            double width = Double.NEGATIVE_INFINITY;
            double lastCutWidth = 0;
            // sum of the row height maximums
            double height = 0;
            double rowHeightMax = heights[0];
            
            if (m == 0) {
                width = total;
                height = gs.getMaxHeight();
            } else {
                while (index < gs.longestPath) {
                    
                    if (widthAtIndex[index - 1] - sumSoFar >= rowSum) {
    
                        // cut _before_ index
                        cuts.add(index);  
                        
                        // update state
                        width = Math.max(width, widthAtIndex[index - 1] - lastCutWidth);
                        height += rowHeightMax;
                        
                        sumSoFar += (widthAtIndex[index - 1] - sumSoFar);
                        lastCutWidth = widthAtIndex[index - 1];
                        rowHeightMax = heights[index];
                    }
                    
                    rowHeightMax = Math.max(rowHeightMax, heights[index]);
                    
                    index++;
                }
                
                // add heights of last row
                height += rowHeightMax;
            }
            
            // are they better?
            double maxScale =  Math.min(1d / width, (1d / gs.dar) / height);
            
            if (maxScale > bestMaxScale) {
                bestMaxScale = maxScale;
                bestCuts = cuts;
            }
        }
        
        return bestCuts;
    }
    
}
