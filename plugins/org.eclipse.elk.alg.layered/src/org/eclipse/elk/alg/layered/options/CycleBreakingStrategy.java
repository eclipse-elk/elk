/*******************************************************************************
 * Copyright (c) 2010, 2017 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.alg.layered.options;

import org.eclipse.elk.alg.layered.LayeredPhases;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.p1cycles.BFSNodeOrderCycleBreaker;
import org.eclipse.elk.alg.layered.p1cycles.DFSNodeOrderCycleBreaker;
import org.eclipse.elk.alg.layered.p1cycles.DepthFirstCycleBreaker;
import org.eclipse.elk.alg.layered.p1cycles.GreedyCycleBreaker;
import org.eclipse.elk.alg.layered.p1cycles.GreedyModelOrderCycleBreaker;
import org.eclipse.elk.alg.layered.p1cycles.InteractiveCycleBreaker;
import org.eclipse.elk.alg.layered.p1cycles.ModelOrderCycleBreaker;
import org.eclipse.elk.alg.layered.p1cycles.SCCNodeTypeCycleBreaker;
import org.eclipse.elk.alg.layered.p1cycles.SCConnectivity;
import org.eclipse.elk.alg.layered.p1cycles.StrictGroupOrderCycleBreaker;
import org.eclipse.elk.core.alg.ILayoutPhase;
import org.eclipse.elk.core.alg.ILayoutPhaseFactory;
import org.eclipse.elk.graph.properties.AdvancedPropertyValue;

/**
 * Enumeration of and factory for the different available cycle breaking strategies.
 * 
 * @author msp
 * @author cds
 * FIXME explain all the strategies again in detail
 */
public enum CycleBreakingStrategy implements ILayoutPhaseFactory<LayeredPhases, LGraph> {

    /**
     * Applies a greedy heuristic to minimize the number of reversed edges.
     */
    GREEDY,
    /**
     * Applies a depth-first traversal to find and reverse edges to make the graph acyclic.
     * This uses the edge order as the iteration order.
     */
    DEPTH_FIRST,
    /**
     * Reacts on user interaction by respecting initial node positions. The actual positions
     * as given in the input diagram are considered here. This means that if the user moves
     * a node, that movement is reflected in the decision which edges to reverse.
     */
    @AdvancedPropertyValue
    INTERACTIVE,
    
    /**
     * Reacts to the input model by respecting the initial ordering in the model file.
     * This ordering is used to identify backwards edges.
     */
    MODEL_ORDER,
    
    /**
     * Applies a greedy heuristic to minimize the number of reversed edges but uses the model order as a tie-breaker.
     */
    GREEDY_MODEL_ORDER,

    /**
     *  Uses Tarjan's algorithm to calculate the strongly connected components. Then it determines the order in these 
     *  components, i.e. the cycles using the in-/out-degree of the nodes. E.g. first Tarjan, then greedy.
     *  FIXME is this really true
     */
    SCC_CONNECTIVITY,

    /**
     *  Uses the Strongly Connected Component approach and selects the minima/maxima based on the model order.
     *  FIXME, what does this really do?
     */
    SCC_NODE_TYPE,

    /**
     *  Strictly enforces a order between different nodes, uses the model order for same types.
     *  FIXME I think this could be done by adding the groupmo strategy to another strategy.
     */
    STRICT_GROUP_ORDER,
    /**
     * Applies a depth-first traversal to find and reverse edges to make the graph acyclic.
     * This uses the node order as the iteration order.
     */
    DFS_NODE_ORDER,

    /**
     * Applies a breadth-first traversal to find and reverse edges to make the graph acyclic.
     * This uses the node order as the iteration order.
     */
    BFS_NODE_ORDER;
    

    @Override
    public ILayoutPhase<LayeredPhases, LGraph> create() {
        switch (this) {
        case GREEDY:
            return new GreedyCycleBreaker();
            
        case DEPTH_FIRST:
            return new DepthFirstCycleBreaker();
            
        case INTERACTIVE:
            return new InteractiveCycleBreaker();
            
        case MODEL_ORDER:
            return new ModelOrderCycleBreaker();
            
        case GREEDY_MODEL_ORDER:
            return new GreedyModelOrderCycleBreaker();

        case SCC_CONNECTIVITY:
            return new SCConnectivity();

        case SCC_NODE_TYPE:
            return new SCCNodeTypeCycleBreaker();

        case STRICT_GROUP_ORDER:
            return new StrictGroupOrderCycleBreaker();

        case DFS_NODE_ORDER:
            return new DFSNodeOrderCycleBreaker();

        case BFS_NODE_ORDER:
            return new BFSNodeOrderCycleBreaker();

        default:
            throw new IllegalArgumentException(
                    "No implementation is available for the cycle breaker " + this.toString());
        }
    }

}
