/*******************************************************************************
 * Copyright (c) 2010, 2015 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.layered.options;

import org.eclipse.elk.alg.layered.LayeredPhases;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.p2layers.BreadthFirstModelOrderLayerer;
import org.eclipse.elk.alg.layered.p2layers.CoffmanGrahamLayerer;
import org.eclipse.elk.alg.layered.p2layers.DepthFirstModelOrderLayerer;
import org.eclipse.elk.alg.layered.p2layers.InteractiveLayerer;
import org.eclipse.elk.alg.layered.p2layers.LongestPathLayerer;
import org.eclipse.elk.alg.layered.p2layers.LongestPathSourceLayerer;
import org.eclipse.elk.alg.layered.p2layers.MinWidthLayerer;
import org.eclipse.elk.alg.layered.p2layers.NetworkSimplexLayerer;
import org.eclipse.elk.alg.layered.p2layers.StretchWidthLayerer;
import org.eclipse.elk.core.alg.ILayoutPhase;
import org.eclipse.elk.core.alg.ILayoutPhaseFactory;
import org.eclipse.elk.graph.properties.AdvancedPropertyValue;
import org.eclipse.elk.graph.properties.ExperimentalPropertyValue;

/**
 * Enumeration of and factory for the different available layering strategies.
 */
public enum LayeringStrategy implements ILayoutPhaseFactory<LayeredPhases, LGraph> {

    /**
     * All nodes will be layered with minimal edge length by using the network-simplex-algorithm.
     */
    NETWORK_SIMPLEX,
    /**
     * All nodes will be layered according to the longest path to any sink.
     */
    LONGEST_PATH,
    /**
     * All nodes will be layered according to the longest path to any source.
     */
    LONGEST_PATH_SOURCE,
    /** 
     * Allows to restrict the number of original nodes in any layer. Essentially 
     * solves a precedence-constrained scheduling problem. 
     */
    @AdvancedPropertyValue
    COFFMAN_GRAHAM,
    /**
     * Nodes are put into layers according to their relative position. The actual positions
     * as given in the input diagram are considered here. This means that if the user moves
     * a node, that movement is reflected in the layering of the graph.
     */
    @AdvancedPropertyValue
    INTERACTIVE,
    /**
     * Nodes are put into Layers according to the average out-degree and their rank.
     * Similar to {@link #LONGEST_PATH}, however, tries to reduce the maximum number of 
     * nodes per layer.
     */
    @ExperimentalPropertyValue
    STRETCH_WIDTH,
    /**
     * Implementation of the 'MinWidth' heuristic for solving the NP-hard minimum-width layering
     * problem with consideration of dummy nodes.
     */
    @ExperimentalPropertyValue
    MIN_WIDTH,
    /**
     * Assumes an intended breadth first model order of nodes and creates layers accordingly.
     * No node with a higher model order is in a layer before a node with a lower model order.
     * Dummy labels are only placed in layers between real nodes.
     */
    @ExperimentalPropertyValue
    BF_MODEL_ORDER,
    /**
     * Assumes an intended depth first model order of nodes and creates layers accordingly.
     * Nodes are placed in the next layer until a node does not connect the previous layer.
     * Dummy labels are only placed in layers between real nodes.
     */
    @ExperimentalPropertyValue
    DF_MODEL_ORDER;
    
    @Override
    public ILayoutPhase<LayeredPhases, LGraph> create() {
        switch (this) {
        case NETWORK_SIMPLEX:
            return new NetworkSimplexLayerer();
            
        case LONGEST_PATH:
            return new LongestPathLayerer();
            
        case COFFMAN_GRAHAM:
            return new CoffmanGrahamLayerer();
            
        case INTERACTIVE:
            return new InteractiveLayerer();
            
        case STRETCH_WIDTH:
            return new StretchWidthLayerer();
          
        case MIN_WIDTH:
            return new MinWidthLayerer();
            
        case LONGEST_PATH_SOURCE:
            return new LongestPathSourceLayerer();
            
        case BF_MODEL_ORDER:
            return new BreadthFirstModelOrderLayerer();
            
        case DF_MODEL_ORDER:
            return new DepthFirstModelOrderLayerer();

        default:
            throw new IllegalArgumentException(
                    "No implementation is available for the layerer " + this.toString());
        }
    }

}
