/*******************************************************************************
 * Copyright (c) 2010, 2015 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.layered.intermediate;

import java.util.List;
import java.util.ListIterator;

import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LLabel;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LNode.NodeType;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.core.alg.ILayoutProcessor;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.math.KVectorChain;
import org.eclipse.elk.core.options.PortSide;
import org.eclipse.elk.core.util.IElkProgressMonitor;

/**
 * Removes dummy nodes due to edge splitting (dummy nodes that have the node type
 * {@link NodeType#LONG_EDGE}). If an edge is split into a
 * chain of edges <i>e1, e2, ..., ek</i>, the first edge <i>e1</i> is retained, while the other
 * edges <i>e2, ..., ek</i> are discarded. This fact should be respected by all processors that
 * create dummy nodes: they should always put the original edge as first edge in the chain of edges,
 * so the original edge is restored.
 * 
 * <p>
 * The actual implementation that joins long edges is provided by this class as a public utility method
 * to be used by other processors.
 * </p>
 * 
 * <dl>
 *   <dt>Preconditions:</dt>
 *     <dd>a layered graph</dd>
 *     <dd>nodes are placed</dd>
 *     <dd>edges are routed.</dd>
 *   <dt>Postconditions:</dt>
 *     <dd>there are no dummy nodes of type {@link NodeType#LONG_EDGE} in the graph's layers.</dd>
 *     <dd>the dummy nodes' {@link LNode#getLayer() layer} fields 
 *         have <strong>not</strong> been set to {@code null} though.</dd>
 *   <dt>Slots:</dt>
 *     <dd>After phase 5.</dd>
 *   <dt>Same-slot dependencies:</dt>
 *     <dd>{@link HierarchicalPortOrthogonalEdgeRouter}</dd>
 * </dl>
 *
 * @author cds
 */
public final class LongEdgeJoiner implements ILayoutProcessor<LGraph> {

    @Override
    public void process(final LGraph layeredGraph, final IElkProgressMonitor monitor) {
        monitor.begin("Edge joining", 1);
        
        final boolean addUnnecessaryBendpoints =
                layeredGraph.getProperty(LayeredOptions.UNNECESSARY_BENDPOINTS);
        
        // Iterate through the layers
        for (Layer layer : layeredGraph) {
            // Get a list iterator for the layer's nodes (since we might be
            // removing dummy nodes from it)
            ListIterator<LNode> nodeIterator = layer.getNodes().listIterator();
            
            while (nodeIterator.hasNext()) {
                LNode node = nodeIterator.next();
                
                // Check if it's a dummy edge we're looking for
                if (node.getType() == NodeType.LONG_EDGE) {
                    joinAt(node, addUnnecessaryBendpoints);
                    
                    // Remove the node
                    nodeIterator.remove();
                }
            }
        }
        
        monitor.done();
    }

    /**
     * Joins the edges connected to the given dummy node. The dummy node is then ready to be removed
     * from the graph.
     * 
     * @param longEdgeDummy
     *            the dummy node whose incident edges to join.
     * @param addUnnecessaryBendpoints
     *            {@code true} if a bend point should be added to the edges at the position of the
     *            dummy node.
     */
    public static void joinAt(final LNode longEdgeDummy, final boolean addUnnecessaryBendpoints) {
        // Get the input and output port (of which we assume to have only one, on the western side and
        // on the eastern side, respectively); the incoming edges are retained, and the outgoing edges
        // are discarded
        List<LEdge> inputPortEdges =
            longEdgeDummy.getPorts(PortSide.WEST).iterator().next().getIncomingEdges();
        List<LEdge> outputPortEdges =
            longEdgeDummy.getPorts(PortSide.EAST).iterator().next().getOutgoingEdges();
        int edgeCount = inputPortEdges.size();
        
        // If we are to add unnecessary bend points, we need to know where. We take the position of the
        // first port we find. (It doesn't really matter which port we're using, so we opt to keep it
        // surprisingly simple.)
        KVector unnecessaryBendpoint = longEdgeDummy.getPorts().get(0).getAbsoluteAnchor();
        
        // The following code assumes that edges with the same indices in the two lists originate from
        // the same long edge, which is true for the current implementation of LongEdgeSplitter and
        // HyperedgeDummyMerger
        while (edgeCount-- > 0) {
            // Get the two edges
            LEdge survivingEdge = inputPortEdges.get(0);
            LEdge droppedEdge = outputPortEdges.get(0);
            
            // The surviving edge's target needs to be set to the old target of the dropped edge.
            // However, this doesn't replace the dropped edge with the surviving edge in the list of
            // incoming edges of the (new) target port, but instead appends the surviving edge. That in
            // turn messes with the implicit assumption that edges with the same index on input and
            // output ports of long edge dummies belong to each other. Thus, we need to ensure that the
            // surviving edge is at the correct index in the list of incoming edges. Hence the
            // complicated code below. (KIPRA-1670)
            List<LEdge> targetIncomingEdges = droppedEdge.getTarget().getIncomingEdges();
            int droppedEdgeListIndex = targetIncomingEdges.indexOf(droppedEdge);
            survivingEdge.setTargetAndInsertAtIndex(droppedEdge.getTarget(), droppedEdgeListIndex);
            
            // Remove the dropped edge from the graph
            droppedEdge.setSource(null);
            droppedEdge.setTarget(null);
            
            // Join their bend points and add possibly an unnecessary one
            KVectorChain survivingBendPoints = survivingEdge.getBendPoints();
            
            if (addUnnecessaryBendpoints) {
                survivingBendPoints.add(new KVector(unnecessaryBendpoint));
            }
            
            for (KVector bendPoint : droppedEdge.getBendPoints()) {
                survivingBendPoints.add(new KVector(bendPoint));
            }
            
            // Join their labels
            List<LLabel> survivingLabels = survivingEdge.getLabels();
            for (LLabel label: droppedEdge.getLabels()) {
                survivingLabels.add(label);
            }
            
            // Join their junction points
            KVectorChain survivingJunctionPoints = survivingEdge.getProperty(
                    LayeredOptions.JUNCTION_POINTS);
            KVectorChain droppedJunctionsPoints = droppedEdge.getProperty(
                    LayeredOptions.JUNCTION_POINTS);
            if (droppedJunctionsPoints != null) {
                if (survivingJunctionPoints == null) {
                    survivingJunctionPoints = new KVectorChain();
                    survivingEdge.setProperty(LayeredOptions.JUNCTION_POINTS,
                            survivingJunctionPoints);
                }
                for (KVector jp : droppedJunctionsPoints) {
                    survivingJunctionPoints.add(new KVector(jp));
                }
            }
        }
    }

}
