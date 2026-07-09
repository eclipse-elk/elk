/*******************************************************************************
 * Copyright (c) 2011, 2019 Kiel University and others.
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
package org.eclipse.elk.alg.layered.p3order;

import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;

import org.eclipse.elk.alg.layered.LayeredPhases;
import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LNode.NodeType;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.intermediate.IntermediateProcessorStrategy;
import org.eclipse.elk.alg.layered.options.GraphProperties;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.p3order.counting.IInitializable;
import org.eclipse.elk.core.alg.ILayoutPhase;
import org.eclipse.elk.core.alg.LayoutProcessorConfiguration;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.math.KVectorChain;
import org.eclipse.elk.core.util.IElkProgressMonitor;

/**
 * A crossing minimizer that allows user interaction by respecting previous node positions.
 * 
 * <dl>
 *   <dt>Preconditions:</dt>
 *     <dd>The graph has a proper layering</dd>
 *     <dd>All nodes have at least fixed port sides.</dd>
 *   <dt>Postconditions:</dt>
 *     <dd>The order of nodes in each layer is rearranged according to previous positions give
 *       by the input graph.</dd>
 *     <dd>Long edge dummy nodes have their calculated original position set as property
 *       {@link InternalProperties#ORIGINAL_DUMMY_NODE_POSITION}.</dd>
 * </dl>
 */
public final class InteractiveCrossingMinimizer implements ILayoutPhase<LayeredPhases, LGraph> {

    /** intermediate processing configuration. */
    private static final LayoutProcessorConfiguration<LayeredPhases, LGraph> INTERMEDIATE_PROCESSING_CONFIGURATION =
        LayoutProcessorConfiguration.<LayeredPhases, LGraph>create()
            .addBefore(LayeredPhases.P3_NODE_ORDERING, IntermediateProcessorStrategy.LONG_EDGE_SPLITTER)
            .addBefore(LayeredPhases.P4_NODE_PLACEMENT, IntermediateProcessorStrategy.IN_LAYER_CONSTRAINT_PROCESSOR)
            .addAfter(LayeredPhases.P5_EDGE_ROUTING, IntermediateProcessorStrategy.LONG_EDGE_JOINER);
    
    @Override
    public LayoutProcessorConfiguration<LayeredPhases, LGraph> getLayoutProcessorConfiguration(final LGraph graph) {
        LayoutProcessorConfiguration<LayeredPhases, LGraph> configuration =
                LayoutProcessorConfiguration.createFrom(INTERMEDIATE_PROCESSING_CONFIGURATION);
        
        if (graph.getProperty(InternalProperties.GRAPH_PROPERTIES).contains(GraphProperties.NON_FREE_PORTS)) {
            configuration.addBefore(LayeredPhases.P3_NODE_ORDERING, IntermediateProcessorStrategy.PORT_LIST_SORTER);
        }
        
        return configuration;
    }

    @Override
    public void process(final LGraph layeredGraph, final IElkProgressMonitor monitor) {
        monitor.begin("Interactive crossing minimization", 1);
        
        // Set ID's for each layer since they will be used by the port distribution code to index into arrays
        int layerIndex = 0;
        for (Layer layer : layeredGraph.getLayers()) {
            layer.id = layerIndex++;
        }
        
        LNode[][] nodeOrder = layeredGraph.toNodeArray();
        AbstractBarycenterPortDistributor portDistributor = new NodeRelativePortDistributor(nodeOrder.length);
        IInitializable.init(Arrays.asList(portDistributor), nodeOrder);
        int portCount = 0;
        layerIndex = 0;
        for (Layer layer : layeredGraph) {
            // determine a horizontal position for edge bend points comparison
            double horizPos = 0;
            int nodeCount = 0;
            for (LNode node : layer.getNodes()) {
                if (node.getPosition().x > 0) {
                    horizPos += node.getPosition().x + node.getSize().x / 2;
                    nodeCount++;
                }
                for (LPort port : node.getPorts()) {
                    port.id = portCount++;
                }
            }
            
            if (nodeCount > 0) {
                horizPos /= nodeCount;
            }
            
            // create an array of vertical node positions
            final double[] pos = new double[layer.getNodes().size()];
            int nextIndex = 0;
            for (LNode node : layer) {
                node.id = nextIndex++;
                pos[node.id] = getPos(node, horizPos);
                
                // if we have a long edge dummy node, save the calculated position in a property
                // to be used by the interactive node placer (for dummy nodes other than long edge
                // dummies, we haven't calculated meaningful positions)
                if (node.getType() == NodeType.LONG_EDGE) {
                    node.setProperty(InternalProperties.ORIGINAL_DUMMY_NODE_POSITION, pos[node.id]);
                }
            }
            
            // sort the nodes using the position array
            Collections.sort(layer.getNodes(), new Comparator<LNode>() {
                public int compare(final LNode node1, final LNode node2) {
                    int compare = Double.compare(pos[node1.id], pos[node2.id]);
                    
                    if (compare == 0) {
                        // The two nodes have the same y coordinate. Check for node successor
                        // constraints
                        List<LNode> node1Successors =
                                node1.getProperty(InternalProperties.IN_LAYER_SUCCESSOR_CONSTRAINTS);
                        List<LNode> node2Successors =
                                node2.getProperty(InternalProperties.IN_LAYER_SUCCESSOR_CONSTRAINTS);
                        
                        if (node1Successors.contains(node2)) {
                            return -1;
                        } else if (node2Successors.contains(node1)) {
                            return 1;
                        }
                    }
                    
                    return compare;
                }
            });
            portDistributor.distributePortsWhileSweeping(nodeOrder, layerIndex, true);
            layerIndex++;
        }
        
        monitor.done();
    }
    
    /**
     * Determine a vertical position for the given node.
     * 
     * @param node a node
     * @param horizPos the horizontal position at which to measure (relevant for edges)
     * @return the vertical position used for sorting
     */
    private double getPos(final LNode node, final double horizPos) {
        switch (node.getType()) {
        case LONG_EDGE:
            LEdge edge = (LEdge) node.getProperty(InternalProperties.ORIGIN);
            
            // reconstruct the original bend points from the node annotations
            KVectorChain bendpoints = edge.getProperty(InternalProperties.ORIGINAL_BENDPOINTS);
            if (bendpoints == null) {
                bendpoints = new KVectorChain();
            } else if (edge.getProperty(InternalProperties.REVERSED)) {
                bendpoints = KVectorChain.reverse(bendpoints);
            }
            
            // Check if we can determine the position just by using the source point, if we can determine it
            LPort source = node.getProperty(InternalProperties.LONG_EDGE_SOURCE);
            if (source != null) {
                KVector sourcePoint = source.getAbsoluteAnchor();
                if (horizPos <= sourcePoint.x) {
                    return sourcePoint.y;
                }
                
                bendpoints.addFirst(sourcePoint);
            }
            
            // Check if we can determine the position just by using the target point
            LPort target = node.getProperty(InternalProperties.LONG_EDGE_TARGET);
            if (target != null) {
                KVector targetPoint = target.getAbsoluteAnchor();
                if (targetPoint.x <= horizPos) {
                    return targetPoint.y;
                }
                
                bendpoints.addLast(targetPoint);
            }
            
            // Find the two points along the edge that the horizontal point lies between
            if (bendpoints.size() >= 2) {
                Iterator<KVector> pointIter = bendpoints.iterator();
                KVector point1 = pointIter.next();
                KVector point2 = pointIter.next();
                while (point2.x < horizPos && pointIter.hasNext()) {
                    point1 = point2;
                    point2 = pointIter.next();
                }
                return point1.y + (horizPos - point1.x) / (point2.x - point1.x) * (point2.y - point1.y);
            }
            
            break;
            
        case NORTH_SOUTH_PORT:
            // Get one of the ports the dummy node was created for, and its original node
            LPort originPort = (LPort) node.getPorts().get(0).getProperty(InternalProperties.ORIGIN);
            LNode originNode = originPort.getNode();
            
            switch (originPort.getSide()) {
            case NORTH:
                // Use the position of the node's northern side. This causes northern dummies to
                // have the same y coordinate as the node they were created from if TOP_LEFT is
                // used as the anchor point. We solve this when sorting the nodes by y coordinate
                // by respecting node successor constraints.
                return originNode.getPosition().y;
            
            case SOUTH:
                // Use the position of the node's southern side
                return originNode.getPosition().y + originNode.getSize().y;
            }
            
            break;

        // FIXME What about the other node types?
        }
        
        // the fallback solution is to take the previous position of the node's anchor point
        return node.getInteractiveReferencePoint().y;
    }

}
