/*******************************************************************************
 * Copyright (c) 2018, 2019 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.

 * This Source Code may also be made available under the following Secondary
 * Licenses when the conditions for such availability set forth in the Eclipse
 * Public License v. 2.0 are satisfied: GPL-3.0 which is available at
 * https://www.gnu.org/licenses/gpl-3.0-standalone.html.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.layered.intermediate;

import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LNode.NodeType;
import org.eclipse.elk.alg.layered.intermediate.loops.SelfLoopEdge;
import org.eclipse.elk.alg.layered.intermediate.loops.SelfLoopHolder;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.core.alg.ILayoutProcessor;
import org.eclipse.elk.core.util.IElkProgressMonitor;

/**
 * Finds regular nodes with self loops and postprocesses those loops.
 *
 * <dl>
 *   <dt>Precondition:</dt>
 *     <dd>A layered graph.</dd>
 *     <dd>Self loops are routed.</dd>
 *     <dd>Node coordinates are set.</dd>
 *   <dt>Postcondition:</dt>
 *     <dd>All self loops are restored and routed.</dd>
 *   <dt>Slots:</dt>
 *     <dd>After phase 5.</dd>
 *   <dt>Same-slot dependencies:</dt>
 *     <dd>None.</dd>
 * </dl>
 */
public class SelfLoopPostProcessor implements ILayoutProcessor<LGraph> {

    @Override
    public void process(final LGraph graph, final IElkProgressMonitor progressMonitor) {
        progressMonitor.begin("Self-Loop post-processing", 1);

        graph.getLayers().stream()
            .flatMap(layer -> layer.getNodes().stream())
            .filter(lNode -> lNode.getType() == NodeType.NORMAL)
            .filter(lNode -> lNode.hasProperty(InternalProperties.SELF_LOOP_HOLDER))
            .forEach(lNode -> processNode(lNode));

        progressMonitor.done();
    }

    private void processNode(final LNode lNode) {
        SelfLoopHolder slHolder = lNode.getProperty(InternalProperties.SELF_LOOP_HOLDER);
        
        slHolder.getSLHyperLoops().stream()
            .flatMap(slLoop -> slLoop.getSLEdges().stream())
            .forEach(slEdge -> restoreEdge(lNode, slEdge));
        
        slHolder.getSLHyperLoops().stream()
            .filter(slLoop -> slLoop.getSLLabels() != null)
            .forEach(slLoop -> slLoop.getSLLabels().applyPlacement(lNode.getPosition()));
    }

    private void restoreEdge(final LNode lNode, final SelfLoopEdge slEdge) {
        LEdge lEdge = slEdge.getLEdge();
        lEdge.setSource(slEdge.getSLSource().getLPort());
        lEdge.setTarget(slEdge.getSLTarget().getLPort());
        
        lEdge.getBendPoints().offset(lNode.getPosition());
    }

}
