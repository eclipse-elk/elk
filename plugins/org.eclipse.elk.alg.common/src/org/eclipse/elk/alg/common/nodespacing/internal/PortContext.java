/*******************************************************************************
 * Copyright (c) 2017, 2020 Kiel University and others.
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
package org.eclipse.elk.alg.common.nodespacing.internal;

import org.eclipse.elk.alg.common.nodespacing.cellsystem.LabelCell;
import org.eclipse.elk.core.math.ElkMargin;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.options.PortLabelPlacement;
import org.eclipse.elk.core.util.adapters.GraphAdapters.PortAdapter;

/**
 * Data holder class to be passed around to avoid having too much state in the size calculation classes. Port contexts
 * are part of {@link NodeContext node contexts}. The position of a port calculated as part of the algorithm should 
 * first be stored in {@link #portPosition} and only be applied at the end of the algorithm, if required.
 */
public final class PortContext {
    
    // CHECKSTYLEOFF Visibility
    // This is a purely internal data holder class, so we ditch setters for public fields.
    
    /////////////////////////////////////////////////////////////////////////////////
    // Convenience Access to Things
    
    /** The node the port belongs to. */
    public final NodeContext parentNodeContext;
    /** The port we calculate stuff for. */
    public final PortAdapter<?> port;
    /** The port's position, to be modified by the algorithm and possibly applied later. */
    public final KVector portPosition;
    /** Whether the port's labels need to be placed next to the port. */
    public final boolean labelsNextToPort;


    /////////////////////////////////////////////////////////////////////////////////
    // Calculated Things
    
    /**
     * Margin aroung the port to assume when placing the port. If node labels are taken into consideration, this will
     * for example include the label cell. When placing the ports, this is the size the port will be assumed to have.
     */
    public ElkMargin portMargin = new ElkMargin();
    /** The cell we place our port labels in. */
    public LabelCell portLabelCell;
    
    
    
    /////////////////////////////////////////////////////////////////////////////////
    // Creation
    
    /**
     * Creates a new context object for the given port, fully initialized with the port's settings.
     * 
     * @param parentNodeContext the port's parent node context.
     * @param port the port to create the context for.
     */
    public PortContext(final NodeContext parentNodeContext, final PortAdapter<?> port) {
        this.parentNodeContext = parentNodeContext;
        this.port = port;
        this.portPosition = new KVector(port.getPosition());
        
        final boolean portLabelsNextToPort =
                parentNodeContext.portLabelsPlacement.contains(PortLabelPlacement.NEXT_TO_PORT_IF_POSSIBLE);
        
        // Whether labels are supposed to be placed next to their port is determined differently depending on whether
        // they are to be placed inside or outside
        if (parentNodeContext.portLabelsPlacement.contains(PortLabelPlacement.INSIDE)) {
            if (parentNodeContext.treatAsCompoundNode) {
                // There might be connections to the inside. That means that we may want to place port labels next to
                // their port, if possible
                labelsNextToPort = portLabelsNextToPort && !port.hasCompoundConnections();
                
            } else {
                labelsNextToPort = true;
            }
        } else if (parentNodeContext.portLabelsPlacement.contains(PortLabelPlacement.OUTSIDE)) {
            if (portLabelsNextToPort) {
                // We can place a label next to an outside port if that port doesn't have incident edges
                labelsNextToPort =
                        !(port.getIncomingEdges().iterator().hasNext() || port.getOutgoingEdges().iterator().hasNext());
            } else {
                // We are not to place labels next to the port anyway
                labelsNextToPort = false;
            }
        } else {
            labelsNextToPort = false;
        }
    }
    
    
    /////////////////////////////////////////////////////////////////////////////////
    // Application
    
    /**
     * Applies the port position stored in this context to the actual port.
     */
    public void applyPortPosition() {
        port.setPosition(portPosition);
    }
    
}
