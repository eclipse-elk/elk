/*******************************************************************************
 * Copyright (c) 2022 Kiel University and others.
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
package org.eclipse.elk.alg.layered.components;

import java.util.List;

import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.options.PortSide;

/**
 * Places components such that rows are formed while still maintaining the model order.
 * Since components cannot be expanded to create visual cues the placement does not utilize
 * subrows or stacks inside a row since the resulting placement can be ambiguous.
 */
public class ModelOrderRowGraphPlacer extends SimpleRowGraphPlacer {
    
    @Override
    public void sortComponents(final List<LGraph> components, final LGraph target) {
        // Does nothing since components are already sorted by model order
    }
    
    @Override
    public void placeComponents(final List<LGraph> components, final LGraph target, final double maxRowWidth, final
            double componentSpacing) {
        // Place nodes iteratively into rows while considering the model order.
        double xpos = 0, ypos = 0, highestBox = 0, broadestRow = componentSpacing;
        LGraph lastComponent = null;
        double startXOfRow = 0;
        for (LGraph graph : components) {
            KVector size = graph.getSize();
            if ((xpos + size.x > maxRowWidth
                    && !graph.getProperty(InternalProperties.EXT_PORT_CONNECTIONS).contains(PortSide.NORTH))
                    || (lastComponent != null
                        && lastComponent.getProperty(InternalProperties.EXT_PORT_CONNECTIONS).contains(PortSide.EAST))
                    || graph.getProperty(InternalProperties.EXT_PORT_CONNECTIONS).contains(PortSide.WEST)) {
                // Components with NORTH connection are allowed to violate the width constrain.
                // Place the graph into the next row
                // Previous EAST ports require a new row.
                // A WEST port requires a new row.
                xpos = startXOfRow;
                ypos += highestBox + componentSpacing;
                highestBox = 0;
            }
            KVector offset = graph.getOffset();
            // North ports should be placed such that they don't intersect with prior components.
            if (graph.getProperty(InternalProperties.EXT_PORT_CONNECTIONS).contains(PortSide.NORTH)) {
                xpos = broadestRow + componentSpacing;
            }
            offsetGraph(graph, xpos + offset.x, ypos + offset.y);
            broadestRow = Math.max(broadestRow, xpos + size.x);
            // South ports block of everything below them.
            if (graph.getProperty(InternalProperties.EXT_PORT_CONNECTIONS).contains(PortSide.SOUTH)) {
                startXOfRow = Math.max(startXOfRow, xpos + size.x + componentSpacing);
            }
            offset.reset();
            highestBox = Math.max(highestBox, size.y);
            xpos += size.x + componentSpacing;
            lastComponent = graph;
        }

        target.getSize().x = broadestRow;
        target.getSize().y = ypos + highestBox;
    }
}
