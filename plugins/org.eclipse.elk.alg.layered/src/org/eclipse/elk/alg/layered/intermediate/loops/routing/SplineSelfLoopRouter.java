/*******************************************************************************
 * Copyright (c) 2019 Kiel University and others.
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
package org.eclipse.elk.alg.layered.intermediate.loops.routing;

import org.eclipse.elk.alg.layered.graph.LGraphUtil;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.intermediate.loops.SelfLoopEdge;
import org.eclipse.elk.alg.layered.intermediate.loops.SelfLoopPort;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.layered.p5edges.splines.NubSpline;
import org.eclipse.elk.alg.layered.p5edges.splines.SplinesMath;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.math.KVectorChain;
import org.eclipse.elk.core.options.PortSide;

/**
 * Routes self loops with splines.
 */
public class SplineSelfLoopRouter extends OrthogonalSelfLoopRouter {
    
    /** Dimension of our self loop splines. */
    private static final int DIM = 3;
    /** To avoid magic number warning. */
    private static final double HALF = 0.5;
    
    
    /**
     * Turns a vector chain of orthogonal bend points into polyline bend points by cutting the corners.
     */
    @Override
    protected KVectorChain modifyBendPoints(final SelfLoopEdge slEdge, final EdgeRoutingDirection routingDirection,
            final KVectorChain bendPoints) {
        
        double edgeLabelDistance = LGraphUtil.getIndividualOrInherited(
                slEdge.getSLHyperLoop().getSLHolder().getLNode(), LayeredOptions.SPACING_EDGE_LABEL);
        
        // For the splines to be routed correctly, we also have to include the source and target positions, so we build
        // up the new list of bend points bit by bit
        KVectorChain splineBendPoints = new KVectorChain(relativePortAnchor(slEdge.getSLSource()));
        addSplineControlPoints(slEdge, routingDirection, bendPoints, splineBendPoints, edgeLabelDistance);
        splineBendPoints.add(relativePortAnchor(slEdge.getSLTarget()));
        
        return new NubSpline(true, DIM, splineBendPoints).getBezierCP();
    }
    
    private KVector relativePortAnchor(final SelfLoopPort slPort) {
        LPort lPort = slPort.getLPort();
        return new KVector(lPort.getPosition()).add(lPort.getAnchor());
    }
    
    /**
     * Inserts spline control points between each consecutive pair of bend points as computed by the orthogonal self
     * loop router. The spline control points are slightly offset away from the node.
     */
    private void addSplineControlPoints(final SelfLoopEdge slEdge, final EdgeRoutingDirection routingDirection,
            final KVectorChain orthoBendPoints, final KVectorChain newBendPoints, final double edgeLabelDistance) {
        
        assert orthoBendPoints.size() >= 2;
        
        // We want to insert a new bend point between each pair of consecutive bend points in the old list
        PortSide currPortSide = slEdge.getSLSource().getLPort().getSide();
        KVector firstBP = orthoBendPoints.get(0);
        
        for (int secondBPIndex = 1; secondBPIndex < orthoBendPoints.size(); secondBPIndex++) {
            // Next pair
            KVector secondBP = orthoBendPoints.get(secondBPIndex);
            
            // The first bend point will go straight into our list before we compute a new bend point
            newBendPoints.add(firstBP);
            
            // Compute a middle bend point and move it away from the node a little
            KVector midBP = new KVector(firstBP).add(secondBP).scale(HALF);
            KVector offset = new KVector(SplinesMath.portSideToDirection(currPortSide))
                    .scale(edgeLabelDistance);
            midBP.add(offset);
            
            newBendPoints.add(midBP);
            
            // Advance to the next pair of bend points on the next port side
            firstBP = secondBP;
            currPortSide = routingDirection == EdgeRoutingDirection.CLOCKWISE
                    ? currPortSide.right()
                    : currPortSide.left();
        }
        
        // Add the last of the original bend points
        newBendPoints.add(orthoBendPoints.getLast());
    }
    
}
