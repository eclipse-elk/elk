/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
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
package org.eclipse.elk.alg.common.nodespacing.internal.algorithm;

import org.eclipse.elk.alg.common.nodespacing.cellsystem.LabelCell;
import org.eclipse.elk.alg.common.nodespacing.internal.NodeContext;
import org.eclipse.elk.alg.common.nodespacing.internal.PortContext;
import org.eclipse.elk.core.options.PortLabelPlacement;
import org.eclipse.elk.core.options.PortSide;
import org.eclipse.elk.core.util.adapters.GraphAdapters.PortAdapter;

/**
 * Creates port context objects and assigns volatile IDs to all ports. Also, unless port labels are fixed, the labels
 * are added to the port context label cells.
 */
public final class PortContextCreator {
    
    /**
     * No instance required.
     */
    private PortContextCreator() {
        
    }
    
    
    /**
     * Creates and initializes port context objects for each of the node's ports. This also involves setting up the
     * label cell associated with each port, although the node label placement therein
     * @param nodeContext
     */
    public static void createPortContexts(final NodeContext nodeContext, final boolean ignoreInsidePortLabels) {
        // Hue hue hue... gimme color?
        boolean imPortLabels =
                !ignoreInsidePortLabels || !nodeContext.portLabelsPlacement.contains(PortLabelPlacement.INSIDE);

        int volatileId = 0;
        for (PortAdapter<?> port : nodeContext.node.getPorts()) {
            if (port.getSide() == PortSide.UNDEFINED) {
                throw new IllegalArgumentException("Label and node size calculator can only be used with ports that "
                        + "have port sides assigned.");
            }
            
            port.setVolatileId(volatileId++);
            createPortContext(nodeContext, port, imPortLabels);
        }
    }
    
    /**
     * Creates a port context for the given adapter and initializes it properly.
     */
    private static void createPortContext(final NodeContext nodeContext, final PortAdapter<?> port,
            final boolean imPortLabels) {
        
        PortContext portContext = new PortContext(nodeContext, port);
        nodeContext.portContexts.put(port.getSide(), portContext);
        
        // If the port has labels and if port labels are to be placed, we need to remember them
        if (imPortLabels && !PortLabelPlacement.isFixed(nodeContext.portLabelsPlacement)) {
            portContext.portLabelCell = new LabelCell(nodeContext.labelLabelSpacing);
            port.getLabels().forEach(label -> portContext.portLabelCell.addLabel(label));
        }
    }
    
}
