/*******************************************************************************
 * Copyright (c) 2013, 2015 Kiel University and others.
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
package org.eclipse.elk.alg.layered.compound;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LGraphUtil;
import org.eclipse.elk.alg.layered.graph.LLabel;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.layered.p5edges.orthogonal.OrthogonalRoutingGenerator;
import org.eclipse.elk.core.alg.ILayoutProcessor;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.math.KVectorChain;
import org.eclipse.elk.core.util.IElkProgressMonitor;

import com.google.common.base.Predicate;
import com.google.common.collect.Iterables;
import com.google.common.collect.Multimap;
import com.google.common.collect.Sets;

/**
 * Postprocess a compound graph by restoring cross-hierarchy edges that have previously been split
 * by the {@link CompoundGraphPreprocessor}.
 * 
 * <dl>
 *   <dt>Precondition:</dt>
 *     <dd>a compound graph with no layers and no cross-hierarchy edges, but with external ports
 *       and with fully specified layout.</dd>
 *   <dt>Postcondition:</dt>
 *     <dd>a compound graph with no layers and with the original cross-hierarchy edges;
 *       the layout applied to these edges conforms to the rules of the ElkGraph meta model.</dd>
 * </dl>
 *
 * @author msp
 * @author cds
 */
public class CompoundGraphPostprocessor implements ILayoutProcessor<LGraph> {
    
    /**
     * A predicate that checks if a given cross hierarchy edge has junction points.
     */
    private static final Predicate<CrossHierarchyEdge> HAS_JUNCTION_POINTS_PREDICATE =
            new Predicate<CrossHierarchyEdge>() {
        
        public boolean apply(final CrossHierarchyEdge chEdge) {
            KVectorChain jps = chEdge.getEdge().getProperty(LayeredOptions.JUNCTION_POINTS);
            return jps != null && !jps.isEmpty();
        }
        
    };
    
    
    @Override
    public void process(final LGraph graph, final IElkProgressMonitor monitor) {
        monitor.begin("Compound graph postprocessor", 1);
        
        // whether bend points should be added whenever crossing a hierarchy boundary
        boolean addUnnecessaryBendpoints = graph.getProperty(LayeredOptions.UNNECESSARY_BENDPOINTS);
        
        // restore the cross-hierarchy map that was built by the preprocessor
        Multimap<LEdge, CrossHierarchyEdge> crossHierarchyMap = graph.getProperty(
                InternalProperties.CROSS_HIERARCHY_MAP);
        
        // remember all dummy edges we encounter; these need to be removed at the end
        Set<LEdge> dummyEdges = Sets.newHashSet();
        
        // iterate over all original edges
        for (LEdge origEdge : crossHierarchyMap.keySet()) {
            // find all cross-hierarchy edges the original edge was split into, and sort them from source to target
            List<CrossHierarchyEdge> crossHierarchyEdges = new ArrayList<CrossHierarchyEdge>(
                    crossHierarchyMap.get(origEdge));
            Collections.sort(crossHierarchyEdges, new CrossHierarchyEdgeComparator(graph));
           
            // find the original source and target ports for the original edge
            LPort sourcePort = crossHierarchyEdges.get(0).getActualSource();
            LPort targetPort = crossHierarchyEdges.get(crossHierarchyEdges.size() - 1).getActualTarget();
            
            // determine the reference graph for all bend points
            LNode referenceNode = sourcePort.getNode();
            LGraph referenceGraph;
            if (LGraphUtil.isDescendant(targetPort.getNode(), referenceNode)) {
                referenceGraph = referenceNode.getNestedGraph();
            } else {
                referenceGraph = referenceNode.getGraph();
            }

            // check whether there are any junction points
            KVectorChain junctionPoints = clearJunctionPoints(origEdge, crossHierarchyEdges);
            
            // reset bend points (we have computed new ones anyway)
            origEdge.getBendPoints().clear();
            
            // apply the computed layouts to the cross-hierarchy edge
            KVector lastPoint = null;
            for (CrossHierarchyEdge chEdge : crossHierarchyEdges) {
                // transform all coordinates from the graph of the dummy edge to the reference graph
                KVector offset = new KVector();
                LGraphUtil.changeCoordSystem(offset, chEdge.getGraph(), referenceGraph);
                
                LEdge ledge = chEdge.getEdge();
                KVectorChain bendPoints = new KVectorChain();
                bendPoints.addAllAsCopies(0, ledge.getBendPoints());
                bendPoints.offset(offset);
                
                // Note: if an NPE occurs here, that means ELK Layered has replaced the original edge
                KVector sourcePoint = new KVector(ledge.getSource().getAbsoluteAnchor());
                KVector targetPoint = new KVector(ledge.getTarget().getAbsoluteAnchor());
                sourcePoint.add(offset);
                targetPoint.add(offset);

                if (lastPoint != null) {
                    KVector nextPoint;
                    if (bendPoints.isEmpty()) {
                        nextPoint = targetPoint;
                    } else {
                        nextPoint = bendPoints.getFirst();
                    }
                    
                    // we add the source point as a bend point to properly connect the hierarchy levels
                    // either if the last point of the previous hierarchy edge segment is a certain
                    // level of tolerance away or if we are required to add unnecessary bend points
                    boolean xDiffEnough =
                            Math.abs(lastPoint.x - nextPoint.x) > OrthogonalRoutingGenerator.TOLERANCE;
                    boolean yDiffEnough =
                            Math.abs(lastPoint.y - nextPoint.y) > OrthogonalRoutingGenerator.TOLERANCE;
                    
                    if ((!addUnnecessaryBendpoints && xDiffEnough && yDiffEnough)
                            || (addUnnecessaryBendpoints && (xDiffEnough || yDiffEnough))) {
                        
                        origEdge.getBendPoints().add(sourcePoint);
                    }
                }

                origEdge.getBendPoints().addAll(bendPoints);
                
                if (bendPoints.isEmpty()) {
                    lastPoint = sourcePoint;
                } else {
                    lastPoint = bendPoints.getLast();
                }
                
                // copy junction points
                copyJunctionPoints(ledge, junctionPoints, offset);
                
                // add offset to target port with a special property
                if (chEdge.getActualTarget() == targetPort) {
                    if (targetPort.getNode().getGraph() != chEdge.getGraph()) {
                        // the target port is in a different coordinate system -- recompute the offset
                        offset = new KVector();
                        LGraphUtil.changeCoordSystem(offset, targetPort.getNode().getGraph(),
                                referenceGraph);
                    }
                    origEdge.setProperty(InternalProperties.TARGET_OFFSET, offset);
                }
                
                // copy labels back to the original edge
                copyLabelsBack(ledge, origEdge, referenceGraph);
                
                // remember the dummy edge for later removal (dummy edges may be in use by several
                // different original edges, which is why we cannot just go and remove it now)
                dummyEdges.add(ledge);
            }
            
            // restore the original source port and target port
            origEdge.setSource(sourcePort);
            origEdge.setTarget(targetPort);
        }
        
        // remove the dummy edges from the graph (dummy ports and dummy nodes are retained)
        for (LEdge dummyEdge : dummyEdges) {
            dummyEdge.setSource(null);
            dummyEdge.setTarget(null);
        }

        monitor.done();
    }
    
    /**
     * Clears an original edge's list of junction points and returns them. The returned vector chain is already set
     * on the original edge, so all the calling code needs to do to update the edge's junction points is to update
     * the returned vector chain.
     */
    private static KVectorChain clearJunctionPoints(final LEdge origEdge,
            final List<CrossHierarchyEdge> crossHierarchyEdges) {
        
        KVectorChain junctionPoints = origEdge.getProperty(LayeredOptions.JUNCTION_POINTS);
        if (Iterables.any(crossHierarchyEdges, HAS_JUNCTION_POINTS_PREDICATE)) {
            // if so, make sure the original edge has an empty non-null junction point list
            if (junctionPoints == null) {
                junctionPoints = new KVectorChain();
                origEdge.setProperty(LayeredOptions.JUNCTION_POINTS, junctionPoints);
            } else {
                junctionPoints.clear();
            }
        } else if (junctionPoints != null) {
            origEdge.setProperty(LayeredOptions.JUNCTION_POINTS, null);
        }
        
        return junctionPoints;
    }

    /**
     * Copies the junction points of the source to the target, adding the given offset.
     */
    private static void copyJunctionPoints(final LEdge source, final KVectorChain target, final KVector offset) {
        KVectorChain ledgeJPs = source.getProperty(LayeredOptions.JUNCTION_POINTS);
        if (ledgeJPs != null) {
            KVectorChain jpCopies = new KVectorChain();
            jpCopies.addAllAsCopies(0, ledgeJPs);
            jpCopies.offset(offset);
            
            target.addAll(jpCopies);
        }
    }

    /**
     * Copies the labels from the given hierarchy segment back to the original hierarchical edge.
     */
    private static void copyLabelsBack(final LEdge hierarchySegment, final LEdge origEdge,
            final LGraph referenceGraph) {
        
        Iterator<LLabel> labelIterator = hierarchySegment.getLabels().listIterator();
        while (labelIterator.hasNext()) {
            LLabel currLabel = labelIterator.next();
            if (currLabel.getProperty(InternalProperties.ORIGINAL_LABEL_EDGE) != origEdge) {
                continue;
            }
            
            LGraphUtil.changeCoordSystem(
                    currLabel.getPosition(),
                    hierarchySegment.getSource().getNode().getGraph(),
                    referenceGraph);
            labelIterator.remove();
            origEdge.getLabels().add(currLabel);
        }
    }

}
