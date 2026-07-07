/*******************************************************************************
 * Copyright (c) 2010, 2020 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
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
 * Routing strategy for routing layers from west to east.
 */
class WestToEastRoutingStrategy extends BaseRoutingDirectionStrategy {
    
    @Override
    public double getPortPositionOnHyperNode(final LPort port) {
        return port.getNode().getPosition().y + port.getPosition().y + port.getAnchor().y;
    }

    @Override
    public PortSide getSourcePortSide() {
        return PortSide.EAST;
    }

    @Override
    public PortSide getTargetPortSide() {
        return PortSide.WEST;
    }

    @Override
    public void calculateBendPoints(final HyperEdgeSegment segment, final double startPos, final double edgeSpacing) {
        // We don't do anything with dummy segments; they are dealt with when their partner is processed
        if (segment.isDummy()) {
            return;
        }
        
        // Calculate coordinates for each port's bend points
        double segmentX = startPos + segment.getRoutingSlot() * edgeSpacing;

        for (LPort port : segment.getPorts()) {
            double sourceY = port.getAbsoluteAnchor().y;

            for (LEdge edge : port.getOutgoingEdges()) {
                if (!edge.isSelfLoop()) {
                    LPort target = edge.getTarget();
                    double targetY = target.getAbsoluteAnchor().y;
                    
                    if (Math.abs(sourceY - targetY) > OrthogonalRoutingGenerator.TOLERANCE) {
                        // We'll update these if we find that the segment was split
                        double currentX = segmentX;
                        HyperEdgeSegment currentSegment = segment;
                        
                        KVector bend = new KVector(currentX, sourceY);
                        edge.getBendPoints().add(bend);
                        addJunctionPointIfNecessary(edge, currentSegment, bend, true);
                        
                        // If this segment was split, we need two additional bend points
                        HyperEdgeSegment splitPartner = segment.getSplitPartner();
                        if (splitPartner != null) {
                            double splitY = splitPartner.getIncomingConnectionCoordinates().get(0);
                            
                            bend = new KVector(currentX, splitY);
                            edge.getBendPoints().add(bend);
                            addJunctionPointIfNecessary(edge, currentSegment, bend, true);
                            
                            // Advance to the split partner's routing slot
                            currentX = startPos + splitPartner.getRoutingSlot() * edgeSpacing;
                            currentSegment = splitPartner;
                            
                            bend = new KVector(currentX, splitY);
                            edge.getBendPoints().add(bend);
                            addJunctionPointIfNecessary(edge, currentSegment, bend, true);
                        }

                        bend = new KVector(currentX, targetY);
                        edge.getBendPoints().add(bend);
                        addJunctionPointIfNecessary(edge, currentSegment, bend, true);
                    }
                }
            }
        }
    }
    
}