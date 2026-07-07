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
import org.eclipse.elk.alg.layered.intermediate.loops.ordering.PortRestorer;
import org.eclipse.elk.alg.layered.intermediate.loops.ordering.PortSideAssigner;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.core.alg.ILayoutProcessor;
import org.eclipse.elk.core.util.IElkProgressMonitor;

/**
 * Restores self loop ports and computes self loop types. Does not restore the self loops themselves.
 * 
 * <dl>
 * <dt>Preconditions:</dt>
 *   <dd>A layered graph whose self-loop nodes have been preprocessed by {@link SelfLoopPreprocessor}.</dd>
 * <dt>Postconditions:</dt>
 *   <dd>All ports previously removed are restored.</dd>
 * <dt>Slots:</dt>
 *   <dd>Before phase 4.</dd>
 * <dt>Same-slot dependencies:</dt>
 *   <dd>None</dd>
 * </dl>
 */
public class SelfLoopPortRestorer implements ILayoutProcessor<LGraph> {
    
    /** This thing will assign hidden ports to port sides if port constraints were initially free. */
    private final PortSideAssigner portSideAssigner = new PortSideAssigner();
    /** Port restorer we'll use for fixed side port constraints. */
    private final PortRestorer portRestorer = new PortRestorer();

    @Override
    public void process(final LGraph graph, final IElkProgressMonitor progressMonitor) {
        progressMonitor.begin("Self-Loop ordering", 1);
        
        // Process every node that actually has self loops
        graph.getLayers().stream()
            .flatMap(layer -> layer.getNodes().stream())
            .filter(lNode -> lNode.getType() == NodeType.NORMAL)
            .filter(lNode -> lNode.hasProperty(InternalProperties.SELF_LOOP_HOLDER))
            .map(lNode -> lNode.getProperty(InternalProperties.SELF_LOOP_HOLDER))
            .forEach(slHolder -> processNode(slHolder, progressMonitor));
        
        progressMonitor.done();
    }

    private void processNode(final SelfLoopHolder slHolder, final IElkProgressMonitor monitor) {
        // Restore and order pure self loop ports if they were previously hidden
        if (slHolder.arePortsHidden()) {
            switch (slHolder.getLNode().getProperty(InternalProperties.ORIGINAL_PORT_CONSTRAINTS)) {
            case UNDEFINED:
            case FREE:
                // We need to assign port sides first and then fall through to restore ports
                portSideAssigner.assignPortSides(slHolder);
                
            case FIXED_SIDE:
                // Restore ports (which by now have port sides assigned to them). After this call, arePortsHidden() will
                // report false
                computeSelfLoopTypes(slHolder);
                portRestorer.restorePorts(slHolder, monitor);
                break;
                
            default:
                // This should not happen. If ports were hidden this must have been because their order was not fixed
                assert false;
            }
            
        } else {
            // Ensure that self loops types are computed
            computeSelfLoopTypes(slHolder);
        }
    }
    
    private void computeSelfLoopTypes(final SelfLoopHolder slHolder) {
        slHolder.getSLHyperLoops().stream().forEach(slLoop -> slLoop.computePortsPerSide());
    }

}
