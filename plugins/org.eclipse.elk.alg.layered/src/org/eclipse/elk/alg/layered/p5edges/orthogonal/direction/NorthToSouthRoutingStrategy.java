/*******************************************************************************
 * Copyright (c) 2010, 2020 Kiel University and others.
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
package org.eclipse.elk.alg.layered.p5edges.orthogonal.direction;

import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.p5edges.orthogonal.HyperEdgeSegment;
import org.eclipse.elk.alg.layered.p5edges.orthogonal.OrthogonalRoutingGenerator;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.options.PortSide;

/**
 * Routing strategy for routing layers from north to south.
 */
class NorthToSouthRoutingStrategy extends BaseRoutingDirectionStrategy {

    @Override
    public double getPortPositionOnHyperNode(final LPort port) {
        return port.getNode().getPosition().x + port.getPosition().x + port.getAnchor().x;
    }

    @Override
    public PortSide getSourcePortSide() {
        return PortSide.SOUTH;
    }

    @Override
    public PortSide getTargetPortSide() {
        return PortSide.NORTH;
    }

    @Override
    public void calculateBendPoints(final HyperEdgeSegment segment, final double startPos, final double edgeSpacing) {
        // We don't do anything with dummy segments; they are dealt with when their partner is processed
        if (segment.isDummy()) {
            return;
        }
        
        // Calculate coordinates for each port's bend points
        double segmentY = startPos + segment.getRoutingSlot() * edgeSpacing;

        for (LPort port : segment.getPorts()) {
            double sourceX = port.getAbsoluteAnchor().x;

            for (LEdge edge : port.getOutgoingEdges()) {
                if (!edge.isSelfLoop()) {
                    LPort target = edge.getTarget();
                    double targetX = target.getAbsoluteAnchor().x;
                    
                    if (Math.abs(sourceX - targetX) > OrthogonalRoutingGenerator.TOLERANCE) {
                        // We'll update these if we find that the segment was split
                        double currentY = segmentY;
                        HyperEdgeSegment currentSegment = segment;
                        
                        KVector bend = new KVector(sourceX, currentY);
                        edge.getBendPoints().add(bend);
                        addJunctionPointIfNecessary(edge, currentSegment, bend, false);
                        
                        // If this segment was split, we need two additional bend points
                        HyperEdgeSegment splitPartner = segment.getSplitPartner();
                        if (splitPartner != null) {
                            double splitX = splitPartner.getIncomingConnectionCoordinates().get(0);
                            
                            bend = new KVector(splitX, currentY);
                            edge.getBendPoints().add(bend);
                            addJunctionPointIfNecessary(edge, currentSegment, bend, false);
                            
                            // Advance to the split partner's routing slot
                            currentY = startPos + splitPartner.getRoutingSlot() * edgeSpacing;
                            currentSegment = splitPartner;
                            
                            bend = new KVector(splitX, currentY);
                            edge.getBendPoints().add(bend);
                            addJunctionPointIfNecessary(edge, currentSegment, bend, false);
                        }

                        bend = new KVector(targetX, currentY);
                        edge.getBendPoints().add(bend);
                        addJunctionPointIfNecessary(edge, currentSegment, bend, false);
                    }
                }
            }
        }
    }
}