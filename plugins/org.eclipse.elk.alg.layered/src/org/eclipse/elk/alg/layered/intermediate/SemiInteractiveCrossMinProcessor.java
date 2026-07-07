/*******************************************************************************
 * Copyright (c) 2016, 2020 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.layered.intermediate;

import java.util.Optional;

import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LNode.NodeType;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.core.alg.ILayoutProcessor;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.util.IElkProgressMonitor;

/**
 * As opposed to the fully interactive {@link org.eclipse.elk.alg.layered.p3order.InteractiveCrossingMinimizer
 * InteractiveCrossingMinimizer}, this processor inserts ordering constraints between pairs of regular nodes 
 * ({@link NodeType#NORMAL}) of the same layer. Apart from this it relies on a crossing minimizer that supports 
 * such kind of ordering constraints, e.g. our {@link org.eclipse.elk.alg.layered.p3order.LayerSweepCrossingMinimizer
 * LayerSweepCrossingMinimizer}. The processor must be run after the layering phase.
 * 
 * As opposed to the interactive crossing minimizer this processor relies on the {@link LayeredOptions#POSITION}
 * property's y-coordinates to derive the desired order. A user may only want to specify the relative order of 
 * certain pairs of nodes. This is supported by only enforcing an order between nodes with the 
 * {@link LayeredOptions#POSITION} option set. 
 * 
 * Note that the {@link org.eclipse.elk.alg.layered.intermediate.greedyswitch.GreedySwitchHeuristic 
 * GreedySwitchHeuristic} is not compatible with this processor and may break the user-defined order.
 * 
 * <dl>
 *   <dt>Precondition:</dt><dd>a layered graph.</dd>
 *   <dt>Postcondition:</dt><dd>inserted in-layer constraints to preserve regular nodes' order</dd>
 *   <dt>Slots:</dt><dd>Before phase 3.</dd>
 *   <dt>Same-slot dependencies:</dt><dd>Any layer-altering processor, such as {@link NodePromotion},
 *                                       must be executed prior to this processor.</dd>
 * </dl>
 */
public class SemiInteractiveCrossMinProcessor implements ILayoutProcessor<LGraph> {

    @Override
    public void process(final LGraph layeredGraph, final IElkProgressMonitor progressMonitor) {
        progressMonitor.begin("Semi-Interactive Crossing Minimization Processor", 1);

        boolean addedConstraints = false;
        
        for (Layer l : layeredGraph) {
            // #1 extract relevant nodes
            // #2 sort them with ascending y coordinate
            // #3 introduce pair-wise in-layer constraints
            Optional<LNode> reduced = l.getNodes().stream()
                .filter(n -> n.getType() == NodeType.NORMAL)
                .filter(n -> n.getAllProperties().containsKey(LayeredOptions.POSITION))
                .sorted((n1, n2) -> {
                    KVector origPos1 = n1.getProperty(LayeredOptions.POSITION);
                    KVector origPos2 = n2.getProperty(LayeredOptions.POSITION);
                    return Double.compare(origPos1.y, origPos2.y);
                })
                .reduce((prev, cur) -> {
                    prev.getProperty(InternalProperties.IN_LAYER_SUCCESSOR_CONSTRAINTS).add(cur);
                    return cur;
                });
            addedConstraints |= reduced.isPresent();
        }
        
        // If we added in-layer successor constraints, make subsequent phases aware
        if (addedConstraints) {
            layeredGraph.setProperty(InternalProperties.IN_LAYER_SUCCESSOR_CONSTRAINTS_BETWEEN_NON_DUMMIES, true);
        }

        progressMonitor.done();
    }
}
