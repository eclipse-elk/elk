/*******************************************************************************
 * Copyright (c) 2018, 2019 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.layered.intermediate;

import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.intermediate.loops.SelfLoopHolder;
import org.eclipse.elk.alg.layered.intermediate.loops.SelfLoopPort;
import org.eclipse.elk.alg.layered.options.GraphProperties;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.core.alg.ILayoutProcessor;
import org.eclipse.elk.core.options.PortConstraints;
import org.eclipse.elk.core.util.IElkProgressMonitor;

/**
 * Finds regular nodes with self loops and preprocesses those loops.
 *
 * <dl>
 *   <dt>Precondition:</dt>
 *     <dd>A layered graph.</dd>
 *   <dt>Postcondition:</dt>
 *     <dd>Each node with self loops has a {@link SelfLoopHolder} created for it and stored under the
 *         {@link InternalProperties#SELF_LOOP_HOLDER} property.
 *     <dd>All self loop edges are removed from their ports to not confuse subsequent phases.</dd>
 *     <dd>Unless port orders are fixed, all ports with only self loop edges connected to them are removed from their
 *         nodes in order to not confuse subsequent phases.</dd>
 *   <dt>Slots:</dt>
 *     <dd>Before phase 2.</dd>
 *   <dt>Same-slot dependencies:</dt>
 *     <dd>None.</dd>
 * </dl>
 */
public class SelfLoopPreProcessor implements ILayoutProcessor<LGraph> {

    @Override
    public void process(final LGraph graph, final IElkProgressMonitor progressMonitor) {
        progressMonitor.begin("Self-Loop pre-processing", 1);

        for (LNode lnode : graph.getLayerlessNodes()) {
            if (SelfLoopHolder.needsSelfLoopProcessing(lnode)) {
                SelfLoopHolder slHolder = SelfLoopHolder.install(lnode);
                hideSelfLoops(slHolder);
                hidePorts(slHolder);
            }
        }

        progressMonitor.done();
    }

    /**
     * Hides all self loop edges. These would confuse subsequent phases and are thus simply removed from the graph, to
     * be restored once edge routing has finished.
     */
    private void hideSelfLoops(final SelfLoopHolder slHolder) {
        slHolder.getSLHyperLoops().stream()
            .flatMap(slLoop -> slLoop.getSLEdges().stream())
            .map(slEdge -> slEdge.getLEdge())
            .forEach(lEdge -> hideSelfLoop(lEdge));
    }

    /**
     * Hides the given edge by removing it from its ports.
     */
    private void hideSelfLoop(final LEdge lEdge) {
        lEdge.setSource(null);
        lEdge.setTarget(null);
    }

    /**
     * Possibly hides all ports whose only incident edges are self loops. This is only done if port constraints are not
     * at least {@link PortConstraints#FIXED_ORDER}.
     */
    private void hidePorts(final SelfLoopHolder slHolder) {
        LNode lNode = slHolder.getLNode();
        LGraph nestedGraph = lNode.getNestedGraph();
        
        /* There are two cases in which we want to refrain from hiding ports:
         * 1. The port order is already fixed. Then the code that would compute a proper port order won't be confused by
         *    the self-loop ports, so there's no need to hide them.
         * 2. The self loop holder has another graph inside of it which contains external ports. In that case, the
         *    crossing minimization will compute proper port orders. The self loop holder will have the port order
         *    applied and its port constraints set to FIXED_POS, so this is basically the first case with extra steps.
         */
        boolean orderFixed = lNode.getProperty(LayeredOptions.PORT_CONSTRAINTS).isOrderFixed();
        boolean hierarchyMode = nestedGraph != null && nestedGraph.getProperty(InternalProperties.GRAPH_PROPERTIES)
                .contains(GraphProperties.EXTERNAL_PORTS);
        
        if (orderFixed || hierarchyMode) {
            // No need to hide any ports
            return;
        }
        
        for (SelfLoopPort slPort : slHolder.getSLPortMap().values()) {
            if (slPort.hadOnlySelfLoops()) {
                // Hide the port
                LPort lPort = slPort.getLPort();
                lPort.setNode(null);
                
                // Remember that we actually did so
                slPort.setHidden(true);
                slHolder.setPortsHidden(true);
                
                // Originally, we removed external port dummy this port belongs to, if any (#352). Now, we should never
                // run into this case.
                assert lPort.getProperty(InternalProperties.PORT_DUMMY) == null;
            }
        }
    }

}
