/*******************************************************************************
 * Copyright (c) 2019 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.alg.layered.intermediate.loops.routing;

import java.util.Iterator;

import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.intermediate.loops.SelfLoopEdge;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.math.KVectorChain;
import org.eclipse.elk.alg.layered.compaction.oned.CompareFuzzy;


/**
 * Routes self loops with polylines. This is basically the same routing as computed by the
 * {@link OrthogonalSelfLoopRouter}, but with the corners cut a little.
 */
public class PolylineSelfLoopRouter extends OrthogonalSelfLoopRouter {
    
    /** Distance by which corner bend points are offset to cut their corner. */
    private static final double CORNER_DISTANCE = 10;

    /**
     * Turns a vector chain of orthogonal bend points into polyline bend points by cutting the corners.
     */
    @Override
    protected KVectorChain modifyBendPoints(final SelfLoopEdge slEdge, final EdgeRoutingDirection routingDirection,
            final KVectorChain bendPoints) {
        
        // Add the source and target points
        LPort lSourcePort = slEdge.getSLSource().getLPort();
        bendPoints.add(0, lSourcePort.getPosition().clone().add(lSourcePort.getAnchor()));

        LPort lTargetPort = slEdge.getSLTarget().getLPort();
        bendPoints.add(lTargetPort.getPosition().clone().add(lTargetPort.getAnchor()));
        
        return cutCorners(bendPoints, CORNER_DISTANCE);
    }
    
    /** Tolerance for double comparisons. */
    private static final double TOLERANCE = 0.01;
    
    /**
     * Turns the given list of bend points into a new list of bend points by cutting the corners. The first and last
     * bend point always remain untouched (and are not included in the returned list of bend points anyway). The other
     * bend points are replaced by two each which are ideally the given distance away from the original bend point.
     * 
     * This method is public for it to be accessible by our unit tests.
     */
    public KVectorChain cutCorners(final KVectorChain bendPoints, final double distance) {
        // The incoming list should consist of more than just the two end points
        assert bendPoints.size() > 2;
        
        KVectorChain result = new KVectorChain();
        
        // We always look at three bend points at a time, where the one in the middle is the corner we want to cut (this
        // could be improved if KVectorChain wasn't a LinkedList and thus had faster direct access to elements)
        Iterator<KVector> bpIterator = bendPoints.iterator();
        KVector previous = null;
        KVector corner = bpIterator.next();
        KVector next = bpIterator.next();
        
        while (bpIterator.hasNext()) {
            // Move to the next corner
            previous = corner;
            corner = next;
            next = bpIterator.next();
            
            // Orthogonal routing
            assert areOrthogonallyRouted(previous, corner, next);
            
            // Compute the new bend points step by step. Start by finding out how much we need to offset the corner to
            // get to the previous and next bend points. Then, limit how far we travel along that vector. Through all of
            // this, remember that the offsets always have one coordinate at 0.
            KVector offset1 = nearZeroToZero(previous.clone().sub(corner));
            KVector offset2 = nearZeroToZero(next.clone().sub(corner));
            
            // We usually use the standard distance, but that might be too much
            double effectiveDistance = distance;
            effectiveDistance = Math.min(effectiveDistance, Math.abs(offset1.x + offset1.y) / 2);
            effectiveDistance = Math.min(effectiveDistance, Math.abs(offset2.x + offset2.y) / 2);
            
            // Limit the offset vectors to our effective distance, either in the positive or negative direction
            offset1.x = Math.signum(offset1.x) * effectiveDistance;
            offset1.y = Math.signum(offset1.y) * effectiveDistance;
            offset2.x = Math.signum(offset2.x) * effectiveDistance;
            offset2.y = Math.signum(offset2.y) * effectiveDistance;
            
            // Compute the effective bend points and add them to our result
            result.add(offset1.add(corner));
            result.add(offset2.add(corner));
        }
        
        return result;
    }

    /**
     * Determines whether the three vectors constitute an orthogonal routing, minus double tolerances.
     */
    private boolean areOrthogonallyRouted(final KVector previous, final KVector corner, final KVector next) {
        boolean verticalHorizontal = CompareFuzzy.fuzzyEquals(previous.x, corner.x, TOLERANCE)
                && CompareFuzzy.fuzzyEquals(corner.y, next.y, TOLERANCE);
        boolean horizontalVertical = CompareFuzzy.fuzzyEquals(previous.y, corner.y, TOLERANCE)
                && CompareFuzzy.fuzzyEquals(corner.x, next.x, TOLERANCE);
        
        return verticalHorizontal || horizontalVertical;
    }

    /**
     * Sets coordinates close to zero to exactly zero.
     */
    private KVector nearZeroToZero(final KVector vector) {
        if (vector.x >= -TOLERANCE && vector.x <= TOLERANCE) {
            vector.x = 0;
        }
        
        if (vector.y >= -TOLERANCE && vector.y <= TOLERANCE) {
            vector.y = 0;
        }
        
        return vector;
    }
    
}
