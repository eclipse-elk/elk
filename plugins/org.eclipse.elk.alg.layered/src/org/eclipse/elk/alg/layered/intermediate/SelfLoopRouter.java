/*******************************************************************************
 * Copyright (c) 2019 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.layered.intermediate;

import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode.NodeType;
import org.eclipse.elk.alg.layered.intermediate.loops.SelfLoopHolder;
import org.eclipse.elk.alg.layered.intermediate.loops.routing.AbstractSelfLoopRouter;
import org.eclipse.elk.alg.layered.intermediate.loops.routing.LabelPlacer;
import org.eclipse.elk.alg.layered.intermediate.loops.routing.OrthogonalSelfLoopRouter;
import org.eclipse.elk.alg.layered.intermediate.loops.routing.PolylineSelfLoopRouter;
import org.eclipse.elk.alg.layered.intermediate.loops.routing.RoutingDirector;
import org.eclipse.elk.alg.layered.intermediate.loops.routing.RoutingSlotAssigner;
import org.eclipse.elk.alg.layered.intermediate.loops.routing.SplineSelfLoopRouter;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.core.alg.ILayoutProcessor;
import org.eclipse.elk.core.labels.ILabelManager;
import org.eclipse.elk.core.labels.LabelManagementOptions;
import org.eclipse.elk.core.util.IElkProgressMonitor;

/**
 * Computes bend points for self loops and places self loop labels. This is done in several steps:
 * 
 * <dl>
 * <dt>Preconditions:</dt>
 *   <dd>Self loop ports have been restored.</dd>
 *   <dd>The node's inner margins for ports and port labels have been computed.</dd>
 * <dt>Postconditions:</dt>
 *   <dd>All self loops are properly routed.</dd>
 *   <dd>All self loop labels are placed.</dd>
 * <dt>Slots:</dt>
 *   <dd>Before phase 4.</dd>
 * <dt>Same-slot dependencies:</dt>
 *   <dd>{@link SelfLoopPortRestorer}</dd>
 *   <dd>{@link InnermostNodeMarginCalculator}</dd>
 * </dl>
 */
public class SelfLoopRouter implements ILayoutProcessor<LGraph> {
    
    /** Routing director we'll use to compute how the loops are routed around the node. */
    private final RoutingDirector routingDirector = new RoutingDirector();
    /** Label placer we'll use to determine the placement of labels. */
    private final LabelPlacer labelPlacer = new LabelPlacer();
    /** The crossing minimizer we use to assign self loops to routing slots. */
    private final RoutingSlotAssigner routingSlotAssigner = new RoutingSlotAssigner();
    
    
    @Override
    public void process(final LGraph graph, final IElkProgressMonitor progressMonitor) {
        progressMonitor.begin("Self-Loop routing", 1);
        
        AbstractSelfLoopRouter router = routerForGraph(graph);
        ILabelManager labelManager = graph.getProperty(LabelManagementOptions.LABEL_MANAGER);
        
        // Process every node that actually has self loops
        graph.getLayers().stream()
            .flatMap(layer -> layer.getNodes().stream())
            .filter(lNode -> lNode.getType() == NodeType.NORMAL)
            .filter(lNode -> lNode.hasProperty(InternalProperties.SELF_LOOP_HOLDER))
            .map(lNode -> lNode.getProperty(InternalProperties.SELF_LOOP_HOLDER))
            .forEach(slHolder -> processNode(slHolder, labelManager, router, progressMonitor));
        
        progressMonitor.done();
    }

    private AbstractSelfLoopRouter routerForGraph(final LGraph graph) {
        switch (graph.getProperty(LayeredOptions.EDGE_ROUTING)) {
        case POLYLINE:
            return new PolylineSelfLoopRouter();
            
        case SPLINES:
            return new SplineSelfLoopRouter();
            
        default:
            return new OrthogonalSelfLoopRouter();
        }
    }
    
    private void processNode(final SelfLoopHolder slHolder, final ILabelManager labelManager,
            final AbstractSelfLoopRouter slRouter, final IElkProgressMonitor monitor) {
        
        // Compute how each hyper loop is routed around the node (from a "leftmost" to a "rightmost" port involved in
        // the loop, but those must be computed and will influence the number of crossings)
        routingDirector.determineLoopRoutes(slHolder);
        
        // Place self loop labels. This will allow the routing slot assigner to make sure that no two overlapping labels
        // end up in the same slot.
        labelPlacer.placeLabels(slHolder, labelManager, monitor);
        
        // Find out which port side each hyper loop appears on and assign routing slots such that the self loop "trunks"
        // (the pieces of the self loop that seem to run around the node, parallel to the respective node boundary) do
        // not intersect
        routingSlotAssigner.assignRoutingSlots(slHolder);
        
        // Finally route the self loops
        slRouter.routeSelfLoops(slHolder);
    }

}
