/*******************************************************************************
 * Copyright (c) 2023 mwr and others.
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
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.layered.options.PortType;
import org.eclipse.elk.core.alg.ILayoutPhase;
import org.eclipse.elk.core.alg.LayoutProcessorConfiguration;
import org.eclipse.elk.core.util.IElkProgressMonitor;

import com.google.common.collect.Lists;

/**
 * This simplistic approach for the Group Model Order problem is based on the Model Order Cycle Breaker. 
 * It uses the GroupID to decide edge reversals between different node-types. For connections within the same type of
 * nodes it uses the normal Model Order Cycle breaker approach. 
 * NOTICE: For languages that alternate frequently between different nodes (e.g. LinguaFranca) this approach degrades
 * into simply using the GroupID without the influence of Model Order. 
 * @author mwr
 *
 */
public class StrictGroupOrderCycleBreaker implements ILayoutPhase<LayeredPhases, LGraph> {

    private int firstSeparateModelOrder;
    private int lastSeparateModelOrder;

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
        
        // Reset FIRST_SEPARATE and LAST_SEPARATE counters.
        firstSeparateModelOrder = 0;
        lastSeparateModelOrder = 0;
        
        // gather edges that point to the wrong direction
        List<LEdge> revEdges = Lists.newArrayList();
        
        // One needs an offset to make sure that the model order of nodes with port constraints is
        // always lower/higher than that of other nodes.
        // E.g. A node with the LAST constraint needs to have a model order m = modelOrder + offset
        // such that m > m(n) with m(n) being the model order of a normal node n (without constraints).
        // Such that the highest model order has to be used as an offset
        int offset = layeredGraph.getLayerlessNodes().size();
        for (LNode node : layeredGraph.getLayerlessNodes()) {
            if (node.hasProperty(InternalProperties.MODEL_ORDER)) {
                offset = Math.max(offset, node.getProperty(InternalProperties.MODEL_ORDER) + 1);
            }
        }
        
        final int offset_final = offset;
        
        for (LNode source : layeredGraph.getLayerlessNodes()) {
            int modelOrderSource = computeConstraintModelOrder(source, offset);
            int groupIDSource = source.getProperty(LayeredOptions.CONSIDER_MODEL_ORDER_GROUP_I_D);
            
            source.getOutgoingEdges().forEach(edge -> {
                LNode target = edge.getTarget().getNode();
                int groupIDTarget = target.getProperty(LayeredOptions.CONSIDER_MODEL_ORDER_GROUP_I_D);
                int modelOrderTarget = computeConstraintModelOrder(target, offset_final);
                // Nodes with -1 as GroupID are external ports, which are never part of a cycle.
                if (groupIDTarget == -1 || groupIDSource == -1) {
                    return;
                }
                // If groups do not match, use groupID to reverse
                if (groupIDTarget < groupIDSource) {
                    revEdges.add(edge);
                } 
                // If the groups match use the Model Order
                else if (modelOrderTarget < modelOrderSource && groupIDTarget == groupIDSource) {
                    revEdges.add(edge);
                }
            });
        }
        
        // reverse the gathered edges
        for (LEdge edge : revEdges) {
            edge.reverse(layeredGraph, true);
            layeredGraph.setProperty(InternalProperties.CYCLIC, true);
        }
        revEdges.clear();
        monitor.done();
    }
    
    /**
     * Set model order to a value such that the constraint is respected and the ordering between nodes with
     * the same constraint is preserved.
     * The order should be FIRST_SEPARATE < FIRST < NORMAL < LAST < LAST_SEPARATE. The offset is used to make sure the 
     * all nodes have unique model orders.
     * @param node The LNode
     * @param offset The offset between FIRST, FIRST_SEPARATE, NORMAL, LAST_SEPARATE, and LAST nodes for unique order
     * @return A unique model order
     */
    private int computeConstraintModelOrder(final LNode node, final int offset) {
        int modelOrder = 0;
        switch (node.getProperty(LayeredOptions.LAYERING_LAYER_CONSTRAINT)) {
        case FIRST_SEPARATE:
            modelOrder = 2 * -offset + firstSeparateModelOrder;
            firstSeparateModelOrder++;
            break;
        case FIRST:
            modelOrder = -offset;
            break;
        case LAST:
            modelOrder = offset;
            break;
        case LAST_SEPARATE:
            modelOrder = 2 * offset + lastSeparateModelOrder;
            lastSeparateModelOrder++;
            break;
        default:
            break;
        }
        if (node.hasProperty(InternalProperties.MODEL_ORDER)) {
            modelOrder += node.getProperty(InternalProperties.MODEL_ORDER);
        }
        return modelOrder;
    }

}
