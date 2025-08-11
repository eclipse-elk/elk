/*******************************************************************************
 * Copyright (c) 2024 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 * 
 * SPDX-License-Identifier: EPL-2.0 
 *******************************************************************************/
package org.eclipse.elk.core.options;

import org.eclipse.elk.graph.ElkNode;

/**
 * Utility functions for reuse across different size approximators.
 *
 */
public class TopdownSizeApproximatorUtil {
    
    /**
     * Dynamically calculate the multiplier to be applied for the side length of the input node based on the number 
     * of children (with and without hierarchy) it and its siblings have. The distribution is mapped to a log scale,
     * which is divided into a number of categories that determine the multiplier.
     * 
     * @param originalGraph the graph to obtain the category multiplier for
     * @return the sidelength multiplier according to the category i.e. Category i => 2^i
     */
    public static double getSizeCategoryMultiplier(final ElkNode originalGraph) {
        ElkNode parent = originalGraph.getParent();
        int thisGraphsSize = getGraphSize(originalGraph);
        
        int CATEGORIES = originalGraph.getProperty(CoreOptions.TOPDOWN_SIZE_CATEGORIES);
        
        
        if (parent != null) {
            // 1. compute distribution of node sizes
            int sizeMinFound = Integer.MAX_VALUE;
            int sizeMaxFound = Integer.MIN_VALUE;
            
            for (ElkNode child : parent.getChildren()) {
                int size = getGraphSize(child);
                
                if (size > sizeMaxFound) {
                    sizeMaxFound = size;
                }
                if (size < sizeMinFound) {
                    sizeMinFound = size;
                }
            }
            
            
            double sizeMin = 1;
            double sizeMax = Math.pow(4, CATEGORIES);
            // shift the range to encompass the largest graph in the local neighbourhood
            if (sizeMaxFound > sizeMax) {
                sizeMax = sizeMaxFound;
            }
            
            // 2. set cutoffs at quarter percentiles on logarithmic scale 
            double x = (Math.log(sizeMax) - Math.log(sizeMin)) / CATEGORIES;
            double factor = Math.exp(x);
            
            // 3. assign node size according to dynamic cutoffs
            double cutoff = sizeMin * factor;
            for (int i = 0; i < CATEGORIES; i++) {
                if (thisGraphsSize < cutoff) {
                    return Math.pow(2, i);
                } else {
                    cutoff *= factor;
                }
            }
            // largest category
            return Math.pow(2, CATEGORIES-1);
            
        } else {
            return 1.0;
        }
        
    }
    
    /**
     * Returns the "size" of the graph defined as the sum of the children's weights.
     * Each simple node containing no children is counted with a weight of 1.
     * Each node with further children is counted with a weight defined in 
     * `CoreOptions.TOPDOWN_SIZE_CATEGORIES_HIERARCHICAL_NODE_WEIGHT`
     * @param originalGraph the graph
     * @return the size of the graph
     */
    public static int getGraphSize(final ElkNode originalGraph) {
        
        int sum = 0;
        
        final int HIERARCHICAL_NODE_WEIGHT = originalGraph.getProperty(CoreOptions.TOPDOWN_SIZE_CATEGORIES_HIERARCHICAL_NODE_WEIGHT);
        for (ElkNode child : originalGraph.getChildren()) {
            if (child.getChildren() != null && child.getChildren().size() > 0) {
                sum += HIERARCHICAL_NODE_WEIGHT;
            } else {
                sum += 1;
            }
        }
        return sum;
    }

}
