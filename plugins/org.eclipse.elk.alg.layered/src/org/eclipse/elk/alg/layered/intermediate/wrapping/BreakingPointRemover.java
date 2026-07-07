/*******************************************************************************
 * Copyright (c) 2016 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.layered.intermediate.wrapping;

import java.util.Collection;
import java.util.List;

import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LNode.NodeType;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.intermediate.LongEdgeJoiner;
import org.eclipse.elk.alg.layered.intermediate.ReversedEdgeRestorer;
import org.eclipse.elk.alg.layered.intermediate.wrapping.BreakingPointInserter.BPInfo;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.layered.p5edges.splines.SplineSegment;
import org.eclipse.elk.core.alg.ILayoutProcessor;
import org.eclipse.elk.core.math.KVectorChain;
import org.eclipse.elk.core.options.EdgeRouting;
import org.eclipse.elk.core.util.IElkProgressMonitor;

import com.google.common.collect.Lists;

/**
 * The processor locates any {@link NodeType#BREAKING_POINT} end dummies 
 * and removes them together with the start dummy and the edge dummies from the graph.
 * Before removing though it transfers the split edge routes to the original (split) edge.
 * Care has to be taken when assembling the new edge route. Different {@link EdgeRouting} styles 
 * require different bend points. 
 * 
 * The {@link LongEdgeJoiner} took care of joining long edges an in-layer dummies,
 * therefore the scenario we are given is:
 * <pre>node --nodeStartEdge--> startMarker --startEndEdge--> endMarker --e--> node</pre>
 * 
 * <dl>
 *   <dt>Preconditions:</dt>
 *     <dd>the graph is fully laid out, in particular, edges have been routed</dd>
 *   <dt>Postconditions:</dt>
 *     <dd>{@link NodeType#BREAKING_POINT}s have been removed.</dd>
 *     <dd>the 'backward' wrapping edges already point into the correct direction, 
 *         the {@link ReversedEdgeRestorer} does <em>not</em> have to revert them.</dd>
 *   <dt>Slots:</dt>
 *     <dd>After phase 5.</dd>
 *   <dt>Same-slot dependencies:</dt>
 *     <dd>Must run after {@link org.eclipse.elk.alg.layered.intermediate.LongEdgeJoiner LongEdgeJoiner}</dd>
 *     <dd>Must run before {@link org.eclipse.elk.alg.layered.intermediate.compaction.HorizontalGraphCompactor
 *         HorizontalGraphCompactor}</dd>
 * </dl>
 * 
 * @see BreakingPointInserter
 * @see BreakingPointProcessor
 */
public class BreakingPointRemover implements ILayoutProcessor<LGraph> {

    private EdgeRouting edgeRouting;
    
    @Override
    public void process(final LGraph graph, final IElkProgressMonitor progressMonitor) {
        progressMonitor.begin("Breaking Point Removing", 1);
        
        edgeRouting = graph.getProperty(LayeredOptions.EDGE_ROUTING);
        for (Layer l : graph.getLayers()) {
            for (LNode node : Lists.newArrayList(l.getNodes())) {
                
                if (BPInfo.isEnd(node)) {
                    BPInfo bpi = node.getProperty(InternalProperties.BREAKING_POINT_INFO);
                    if (bpi.next == null) {
                        remove(graph, bpi);
                    }
                }
            }
        }
        
        progressMonitor.done();
    }
    
    private void remove(final LGraph graph, final BPInfo bpi) {
      
        // assemble the new edge route
        KVectorChain newBends = new KVectorChain();
        
//        System.out.println("NS: " + bpi.nodeStartEdge.getBendPoints());
//        System.out.println("EN: " + bpi.originalEdge.getBendPoints());
        
        // depending on the edge routing style we have to act differently
        switch (edgeRouting) {
        
            // for splines, the final bendpoints have not been computed yet
            //  (which will be done by the 'FinalSplineBendpointsCalculator').
            // so far, only the intermediate nubsplines have been attached via properties
            case SPLINES:  
                // gather the computed spline information 
                List<SplineSegment> s1 = bpi.nodeStartEdge.getProperty(InternalProperties.SPLINE_ROUTE_START);
                List<SplineSegment> s2 = bpi.startEndEdge.getProperty(InternalProperties.SPLINE_ROUTE_START);
                List<SplineSegment> s3 = bpi.originalEdge.getProperty(InternalProperties.SPLINE_ROUTE_START);
                List<LEdge> e1 = bpi.nodeStartEdge.getProperty(InternalProperties.SPLINE_EDGE_CHAIN);
                List<LEdge> e2 = bpi.startEndEdge.getProperty(InternalProperties.SPLINE_EDGE_CHAIN);
                List<LEdge> e3 = bpi.originalEdge.getProperty(InternalProperties.SPLINE_EDGE_CHAIN);

                // join them (... and remember to reverse some of the segments)
                List<SplineSegment> joinedSegments = Lists.newArrayList();
                joinedSegments.addAll(s1);
                s2.forEach(s -> s.inverseOrder = true);
                joinedSegments.addAll(Lists.reverse(s2));
                joinedSegments.addAll(s3);
                
                List<LEdge> joinedEdges = Lists.newArrayList();
                joinedEdges.addAll(e1);
                joinedEdges.addAll(Lists.reverse(e2));
                joinedEdges.addAll(e3);

                // transfer the information to the original edge
                bpi.originalEdge.setProperty(InternalProperties.SPLINE_ROUTE_START, joinedSegments);
                bpi.originalEdge.setProperty(InternalProperties.SPLINE_EDGE_CHAIN, joinedEdges);
                bpi.originalEdge.setProperty(InternalProperties.SPLINE_SURVIVING_EDGE, bpi.originalEdge);
                
                // cleanup
                bpi.nodeStartEdge.setProperty(InternalProperties.SPLINE_ROUTE_START, null);
                bpi.nodeStartEdge.setProperty(InternalProperties.SPLINE_EDGE_CHAIN, null);
                bpi.startEndEdge.setProperty(InternalProperties.SPLINE_ROUTE_START, null);
                bpi.startEndEdge.setProperty(InternalProperties.SPLINE_EDGE_CHAIN, null);
            break;
        
            // after executing the polyline router the edge segments are simply joined,
            // note that the positions of bpi.start and bpi.end must be added as bend points 
            case POLYLINE:
                newBends.addAll(bpi.nodeStartEdge.getBendPoints());
                newBends.add(bpi.start.getPosition());
                newBends.addAll(Lists.reverse(bpi.startEndEdge.getBendPoints()));
                newBends.add(bpi.end.getPosition());
                newBends.addAll(bpi.originalEdge.getBendPoints());
                break;

            // for orthogonal edge routing it is guaranteed, that incident edges of a node  
            // enter the node parallel to the x axis. Thus, the positions 
            // of bpi.start and bpi.end can be dropped since they lie on a straight line
            default: // ORTHOGONAL
                newBends.addAll(bpi.nodeStartEdge.getBendPoints());
                newBends.addAll(Lists.reverse(bpi.startEndEdge.getBendPoints()));
                newBends.addAll(bpi.originalEdge.getBendPoints());
        }

        // restore original edge
        bpi.originalEdge.getBendPoints().clear();
        bpi.originalEdge.getBendPoints().addAll(newBends);
        bpi.originalEdge.setSource(bpi.nodeStartEdge.getSource());
        
        // collect any created junction points (order can be arbitrary)
        KVectorChain junctionPointsOne = bpi.nodeStartEdge.getProperty(LayeredOptions.JUNCTION_POINTS);
        KVectorChain junctionPointsTwo = bpi.startEndEdge.getProperty(LayeredOptions.JUNCTION_POINTS);
        KVectorChain junctionPointsThree = bpi.originalEdge.getProperty(LayeredOptions.JUNCTION_POINTS);
        if (junctionPointsOne != null || junctionPointsTwo != null || junctionPointsThree != null) {
            KVectorChain newJunctionPoints = new KVectorChain();
            addNullSafe(newJunctionPoints, junctionPointsThree);
            addNullSafe(newJunctionPoints, junctionPointsTwo);
            addNullSafe(newJunctionPoints, junctionPointsOne);
            bpi.originalEdge.setProperty(LayeredOptions.JUNCTION_POINTS, newJunctionPoints);
        }
        
        // remove all the dummy stuff
        bpi.startEndEdge.setSource(null);
        bpi.startEndEdge.setTarget(null);
        bpi.nodeStartEdge.setSource(null);
        bpi.nodeStartEdge.setTarget(null);
        bpi.end.setLayer(null);
        bpi.start.setLayer(null);
        
        if (bpi.prev != null) {
            remove(graph, bpi.prev);
        }
    }
    
    private <T> boolean addNullSafe(final Collection<T> container, final Collection<T> toAdd) {
        if (container == null || toAdd == null) {
            return false;
        }
        return container.addAll(toAdd);
    }
}
