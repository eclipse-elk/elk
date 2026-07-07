/*******************************************************************************
 * Copyright (c) 2011, 2015 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.layered.intermediate;

import java.util.List;

import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.core.alg.ILayoutProcessor;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.math.KVectorChain;
import org.eclipse.elk.core.options.PortSide;
import org.eclipse.elk.core.util.IElkProgressMonitor;

import com.google.common.collect.Lists;

/**
 * Improves the placement of hypernodes by moving them such that they replace the join
 * points of connected edges.
 * 
 * <dl>
 *   <dt>Precondition:</dt><dd>a layered graph with all five phases done</dd>
 *   <dt>Postcondition:</dt><dd>the position of some hypernodes as well as some bend
 *     points of connected edges may be changed</dd>
 *   <dt>Slots:</dt><dd>after phase 5</dd>
 * </dl>
 *
 * @author msp
 */
public final class HypernodesProcessor implements ILayoutProcessor<LGraph> {

    @Override
    public void process(final LGraph layeredGraph, final IElkProgressMonitor monitor) {
        monitor.begin("Hypernodes processing", 1);
        
        for (Layer layer : layeredGraph) {
            for (LNode node : layer) {
                if (node.getProperty(LayeredOptions.HYPERNODE) && node.getPorts().size() <= 2) {
                    int topEdges = 0, rightEdges = 0, bottomEdges = 0, leftEdges = 0;
                    for (LPort port : node.getPorts()) {
                        switch (port.getSide()) {
                        case NORTH:
                            topEdges++;
                            break;
                        case EAST:
                            rightEdges++;
                            break;
                        case SOUTH:
                            bottomEdges++;
                            break;
                        case WEST:
                            leftEdges++;
                            break;
                        }
                    }
                    // don't move the node if there are any edges to the top or bottom
                    if (topEdges == 0 && bottomEdges == 0) {
                        moveHypernode(layeredGraph, node, leftEdges <= rightEdges);
                    }
                }
            }
        }
        
        monitor.done();
    }
    
    /**
     * Move the given hypernode either towards the previous layer or towards the next layer.
     * 
     * @param layeredGraph a layered graph
     * @param hypernode a node that is marked with {@link LayeredOptions.HYPERNODE}
     * @param right if true, the node is moved right (to the next layer), else it
     *     is moved left (to the previous layer)
     */
    private void moveHypernode(final LGraph layeredGraph, final LNode hypernode,
            final boolean right) {
        // find edges that constitute the first join point of the hyperedge
        List<LEdge> bendEdges = Lists.newArrayList();
        double bendx = Integer.MAX_VALUE, diffx = Integer.MAX_VALUE, diffy = Integer.MAX_VALUE;
        if (right) {
            bendx = layeredGraph.getSize().x;
            for (LPort port : hypernode.getPorts()) {
                for (LEdge edge : port.getOutgoingEdges()) {
                    if (!edge.getBendPoints().isEmpty()) {
                        KVector firstPoint = edge.getBendPoints().getFirst();
                        if (firstPoint.x < bendx) {
                            diffx = bendx - firstPoint.x;
                            diffy = Integer.MAX_VALUE;
                            bendEdges.clear();
                            bendx = firstPoint.x;
                        }
                        if (firstPoint.x <= bendx) {
                            bendEdges.add(edge);
                            if (edge.getBendPoints().size() > 1) {
                                diffy = Math.min(diffy, Math.abs(edge.getBendPoints().get(1).y
                                        - firstPoint.y));
                            }
                        }
                    }
                }
            }
        } else {
            for (LPort port : hypernode.getPorts()) {
                for (LEdge edge : port.getIncomingEdges()) {
                    if (!edge.getBendPoints().isEmpty()) {
                        KVector lastPoint = edge.getBendPoints().getLast();
                        if (lastPoint.x > bendx) {
                            diffx = lastPoint.x - bendx;
                            diffy = Integer.MAX_VALUE;
                            bendEdges.clear();
                            bendx = lastPoint.x;
                        }
                        if (lastPoint.x >= bendx) {
                            bendEdges.add(edge);
                            if (edge.getBendPoints().size() > 1) {
                                diffy = Math.min(diffy, Math.abs(edge.getBendPoints().get(
                                        edge.getBendPoints().size() - 2).y - lastPoint.y));
                            }
                        }
                    }
                }
            }
        }
        
        if (!bendEdges.isEmpty() && diffx > hypernode.getSize().x / 2
                && diffy > hypernode.getSize().y / 2) {
            // create new ports for the edges
            LPort northPort = new LPort();
            northPort.setNode(hypernode);
            northPort.setSide(PortSide.NORTH);
            northPort.getPosition().x = hypernode.getSize().x / 2;
            LPort southPort = new LPort();
            southPort.setNode(hypernode);
            southPort.setSide(PortSide.SOUTH);
            southPort.getPosition().x = hypernode.getSize().x / 2;
            southPort.getPosition().y = hypernode.getSize().y;
            // replace the first bend point by the new ports
            for (LEdge edge : bendEdges) {
                KVector first, second;
                if (right) {
                    first = edge.getBendPoints().removeFirst();
                    second = edge.getBendPoints().isEmpty() ? edge.getTarget().getAbsoluteAnchor()
                            : edge.getBendPoints().getFirst();
                    if (second.y >= first.y) {
                        edge.setSource(southPort);
                    } else {
                        edge.setSource(northPort);
                    }
                } else {
                    first = edge.getBendPoints().removeLast();
                    second = edge.getBendPoints().isEmpty() ? edge.getSource().getAbsoluteAnchor()
                            : edge.getBendPoints().getLast();
                    if (second.y >= first.y) {
                        edge.setTarget(southPort);
                    } else {
                        edge.setTarget(northPort);
                    }
                }
                // remove junction points that collide with the eliminated bend point
                KVectorChain junctionPoints = edge.getProperty(LayeredOptions.JUNCTION_POINTS);
                if (junctionPoints != null) {
                    junctionPoints.remove(first);
                }
            }
            // move the node to new position
            hypernode.getPosition().x = bendx - hypernode.getSize().x / 2;
        }
    }

}
