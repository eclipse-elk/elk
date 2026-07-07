/*******************************************************************************
 * Copyright (c) 2011, 2015 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.layered.p1cycles;

import java.util.List;

import org.eclipse.elk.alg.layered.LayeredPhases;
import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.intermediate.IntermediateProcessorStrategy;
import org.eclipse.elk.alg.layered.options.PortType;
import org.eclipse.elk.core.alg.ILayoutPhase;
import org.eclipse.elk.core.alg.LayoutProcessorConfiguration;
import org.eclipse.elk.core.util.IElkProgressMonitor;

import com.google.common.collect.Lists;

/**
 * A cycle breaker that responds to user interaction by respecting the direction of
 * edges as given in the original drawing.
 * 
 * <dl>
 *   <dt>Precondition:</dt><dd>none</dd>
 *   <dt>Postcondition:</dt><dd>the graph has no cycles</dd>
 * </dl>
 * 
 * @author msp
 */
public final class InteractiveCycleBreaker implements ILayoutPhase<LayeredPhases, LGraph> {

    /** intermediate processing configuration. */
    private static final LayoutProcessorConfiguration<LayeredPhases, LGraph> INTERMEDIATE_PROCESSING_CONFIGURATION =
        LayoutProcessorConfiguration.<LayeredPhases, LGraph>create()
            .addBefore(LayeredPhases.P1_CYCLE_BREAKING,
                    IntermediateProcessorStrategy.INTERACTIVE_EXTERNAL_PORT_POSITIONER)
            .addAfter(LayeredPhases.P5_EDGE_ROUTING, IntermediateProcessorStrategy.REVERSED_EDGE_RESTORER);

    @Override
    public LayoutProcessorConfiguration<LayeredPhases, LGraph> getLayoutProcessorConfiguration(final LGraph graph) {
        return INTERMEDIATE_PROCESSING_CONFIGURATION;
    }

    @Override
    public void process(final LGraph layeredGraph, final IElkProgressMonitor monitor) {
        monitor.begin("Interactive cycle breaking", 1);
        
        // gather edges that point to the wrong direction
        List<LEdge> revEdges = Lists.newArrayList();
        for (LNode source : layeredGraph.getLayerlessNodes()) {
            source.id = 1;
            double sourcex = source.getInteractiveReferencePoint().x;
            for (LPort port : source.getPorts(PortType.OUTPUT)) {
                for (LEdge edge : port.getOutgoingEdges()) {
                    LNode target = edge.getTarget().getNode();
                    if (target != source) {
                        double targetx = target.getInteractiveReferencePoint().x;
                        if (targetx < sourcex) {
                            revEdges.add(edge);
                        }
                    }
                }
            }
        }
        // reverse the gathered edges
        for (LEdge edge : revEdges) {
            edge.reverse(layeredGraph, true);
        }
        
        // perform an additional check for cycles - maybe we missed something
        // (could happen if some nodes have the same horizontal position)
        revEdges.clear();
        for (LNode node : layeredGraph.getLayerlessNodes()) {
            // unvisited nodes have id = 1
            if (node.id > 0) {
                findCycles(node, revEdges);
            }
        }
        // again, reverse the edges that were marked
        for (LEdge edge : revEdges) {
            edge.reverse(layeredGraph, true);
        }
        
        revEdges.clear();
        monitor.done();
    }
    
    /**
     * Perform a DFS starting on the given node and mark back edges in order to break cycles.
     * 
     * @param node1 a node
     * @param revEdges list of edges that will be reversed
     */
    private void findCycles(final LNode node1, final List<LEdge> revEdges) {
        // nodes with negative id are part of the currently inspected path
        node1.id = -1;
        for (LPort port : node1.getPorts(PortType.OUTPUT)) {
            for (LEdge edge : port.getOutgoingEdges()) {
                LNode node2 = edge.getTarget().getNode();
                if (node1 != node2) {
                    if (node2.id < 0) {
                        // a node of the current path is found --> cycle
                        revEdges.add(edge);
                    } else if (node2.id > 0) {
                        // the node has not been visited yet --> expand the current path
                        findCycles(node2, revEdges);
                    }
                }
            }
        }
        // nodes with id = 0 have been already visited and are ignored if encountered again
        node1.id = 0;
    }

}
