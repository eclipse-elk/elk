/*******************************************************************************
 * Copyright (c) 2021, 2025 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.alg.layered.p1cycles;

import java.util.List;

import org.eclipse.elk.alg.layered.LayeredPhases;
import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.intermediate.IntermediateProcessorStrategy;
import org.eclipse.elk.alg.layered.options.GroupOrderStrategy;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.layered.options.PortType;
import org.eclipse.elk.core.alg.ILayoutPhase;
import org.eclipse.elk.core.alg.LayoutProcessorConfiguration;
import org.eclipse.elk.core.util.IElkProgressMonitor;

import com.google.common.collect.Lists;

/**
 * A cycle breaker that reverses all edges that go against the model order or group model order,
 * i.e. edges from high model order to low model order.
 * 
 * <dl>
 *   <dt>Precondition:</dt>
 *      <dd>no self loops</dd>
 *   <dt>Postcondition:</dt>
 *      <dd>the graph has no cycles</dd>
 * </dl>
 */
public final class ModelOrderCycleBreaker implements ILayoutPhase<LayeredPhases, LGraph> {

    /** intermediate processing configuration. */
    private static final LayoutProcessorConfiguration<LayeredPhases, LGraph> INTERMEDIATE_PROCESSING_CONFIGURATION =
        LayoutProcessorConfiguration.<LayeredPhases, LGraph>create()
            .addAfter(LayeredPhases.P5_EDGE_ROUTING, IntermediateProcessorStrategy.REVERSED_EDGE_RESTORER);

    @Override
    public LayoutProcessorConfiguration<LayeredPhases, LGraph> getLayoutProcessorConfiguration(final LGraph graph) {
        return INTERMEDIATE_PROCESSING_CONFIGURATION;
    }

    @Override
    public void process(final LGraph layeredGraph, final IElkProgressMonitor monitor) {
        monitor.begin("Model order cycle breaking", 1);
        
        // gather edges that point to the wrong direction
        List<LEdge> revEdges = Lists.newArrayList();
        
        // One needs an offset to make sure that the model order of nodes with port constraints is
        // always lower/higher than that of other nodes.
        // E.g. A node with the LAST constraint needs to have a model order m = modelOrder + offset
        // such that m > m(n) with m(n) being the model order of a normal node n (without constraints).
        // Such that the highest model order has to be used as an offset
        int offset = Math.max(layeredGraph.getLayerlessNodes().size(), layeredGraph.getProperty(InternalProperties.MAX_MODEL_ORDER_NODES));
        int bigOffset = offset * layeredGraph.getProperty(InternalProperties.CB_NUM_MODEL_ORDER_GROUPS);
        boolean enforceGroupModelOrder = layeredGraph.getProperty(
                LayeredOptions.CONSIDER_MODEL_ORDER_GROUP_MODEL_ORDER_CB_GROUP_ORDER_STRATEGY) == GroupOrderStrategy.ENFORCED;
        for (LNode source : layeredGraph.getLayerlessNodes()) {
            GroupModelOrderCalculator calculator = new GroupModelOrderCalculator();
            int modelOrderSource = enforceGroupModelOrder
                    ? calculator.computeConstraintGroupModelOrder(source, bigOffset, offset)
                    : calculator.computeConstraintModelOrder(source, offset);
            
            for (LPort port : source.getPorts(PortType.OUTPUT)) {
                for (LEdge edge : port.getOutgoingEdges()) {
                LNode target = edge.getTarget().getNode();
                    int modelOrderTarget = enforceGroupModelOrder
                            ? calculator.computeConstraintGroupModelOrder(target, bigOffset, offset)
                            : calculator.computeConstraintModelOrder(target, offset);
                    if (modelOrderTarget < modelOrderSource) {
                        revEdges.add(edge);
                    }
                }
            }
        }
        
        // reverse the gathered edges
        for (LEdge edge : revEdges) {
            edge.reverse(layeredGraph, true);
            layeredGraph.setProperty(InternalProperties.CYCLIC, true);
        }
        revEdges.clear();
        monitor.done();
    }
}
