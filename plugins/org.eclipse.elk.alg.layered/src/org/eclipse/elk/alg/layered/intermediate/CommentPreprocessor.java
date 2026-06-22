/*******************************************************************************
 * Copyright (c) 2011, 2015 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.alg.layered.intermediate;

import java.util.Iterator;
import java.util.List;
import java.util.ArrayList;

import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LLabel;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.core.alg.ILayoutProcessor;
import org.eclipse.elk.core.options.PortSide;
import org.eclipse.elk.core.util.IElkProgressMonitor;


/**
 * A pre-processor for comment boxes. Looks for comments that have exactly one connection
 * to a normal node and removes them from the graph. Such comments are put either into
 * the {@link InternalProperties#TOP_COMMENTS} or the {@link InternalProperties#BOTTOM_COMMENTS} list
 * of the connected node and processed later by the {@link CommentPostprocessor}.
 * Other comments are processed normally, i.e. they are treated as regular nodes, but
 * their incident edges may be reversed.
 *
 * <dl>
 *   <dt>Precondition:</dt>
 *      <dd>none</dd>
 *   <dt>Postcondition:</dt>
 *      <dd>Comments with only one connection to a port of degree 1 are removed and stored for later
 *      processing.</dd>
 *   <dt>Slots:</dt>
 *      <dd>Before phase 1.</dd>
 * </dl>
 */
public final class CommentPreprocessor implements ILayoutProcessor<LGraph> {

    @Override
    public void process(final LGraph layeredGraph, final IElkProgressMonitor monitor) {
        monitor.begin("Comment pre-processing", 1);
        
        int commentBoxCount = 0;
        
        Iterator<LNode> nodeIter = layeredGraph.getLayerlessNodes().iterator();
        while (nodeIter.hasNext()) {
            LNode node = nodeIter.next();
            if (node.getProperty(LayeredOptions.COMMENT_BOX)) {
                commentBoxCount++;
                
                int edgeCount = 0;
                LEdge edge = null;
                LPort oppositePort = null;
                for (LPort port : node.getPorts()) {
                    edgeCount += port.getDegree();
                    if (port.getIncomingEdges().size() == 1) {
                        edge = port.getIncomingEdges().get(0);
                        oppositePort = edge.getSource();
                    }
                    if (port.getOutgoingEdges().size() == 1) {
                        edge = port.getOutgoingEdges().get(0);
                        oppositePort = edge.getTarget();
                    }
                }
                
                if (edgeCount == 1 && oppositePort.getDegree() == 1
                        && !oppositePort.getNode().getProperty(LayeredOptions.COMMENT_BOX)) {
                    // found a comment that has exactly one connection
                    processBox(node, edge, oppositePort, oppositePort.getNode());
                    nodeIter.remove();
                } else {
                    // reverse edges that are oddly connected
                    List<LEdge> revEdges = new ArrayList<>();
                    for (LPort port : node.getPorts()) {
                        for (LEdge outedge : port.getOutgoingEdges()) {
                            if (!outedge.getTarget().getOutgoingEdges().isEmpty()) {
                                revEdges.add(outedge);
                            }
                        }
                        
                        for (LEdge inedge : port.getIncomingEdges()) {
                            if (!inedge.getSource().getIncomingEdges().isEmpty()) {
                                revEdges.add(inedge);
                            }
                        }
                    }
                    
                    for (LEdge re : revEdges) {
                        re.reverse(layeredGraph, true);
                    }
                }
            }
        }
        
        if (monitor.isLoggingEnabled()) {
            monitor.log("Found " + commentBoxCount + " comment boxes");
        }
        
        monitor.done();
    }
    
    /**
     * Process a comment box by putting it into a property of the corresponding node.
     * 
     * @param box a comment box
     * @param edge the edge that connects the box with the real node
     * @param oppositePort the port of the real node to which the edge is incident
     * @param realNode the normal node that is connected with the comment
     */
    private void processBox(final LNode box, final LEdge edge, final LPort oppositePort,
            final LNode realNode) {
        boolean topFirst, onlyTop = false, onlyBottom = false;
        if (realNode.getProperty(LayeredOptions.PORT_CONSTRAINTS).isSideFixed()) {
            boolean hasNorth = false, hasSouth = false;
            portLoop: for (LPort port1 : realNode.getPorts()) {
                for (LPort port2 : port1.getConnectedPorts()) {
                    if (!port2.getNode().getProperty(LayeredOptions.COMMENT_BOX)) {
                        if (port1.getSide() == PortSide.NORTH) {
                            hasNorth = true;
                            break portLoop;
                        }
                        if (port1.getSide() == PortSide.SOUTH) {
                            hasSouth = true;
                            break portLoop;
                        }
                    }
                }
            }
            onlyTop = hasSouth && !hasNorth;
            onlyBottom = hasNorth && !hasSouth;
        }
        if (!onlyTop && !onlyBottom && !realNode.getLabels().isEmpty()) {
            double labelPos = 0;
            for (LLabel label : realNode.getLabels()) {
                labelPos += label.getPosition().y + label.getSize().y / 2;
            }
            labelPos /= realNode.getLabels().size();
            topFirst = labelPos >= realNode.getSize().y / 2;
        } else {
            topFirst = !onlyBottom;
        }
        
        List<LNode> boxList;
        if (topFirst) {
            // determine the position to use, favoring the top position
            List<LNode> topBoxes = realNode.getProperty(InternalProperties.TOP_COMMENTS);
            if (topBoxes == null) {
                boxList = new ArrayList<>();
                realNode.setProperty(InternalProperties.TOP_COMMENTS, boxList);
            } else if (onlyTop) {
                boxList = topBoxes;
            } else {
                List<LNode> bottomBoxes = realNode.getProperty(InternalProperties.BOTTOM_COMMENTS);
                if (bottomBoxes == null) {
                    boxList = new ArrayList<>();
                    realNode.setProperty(InternalProperties.BOTTOM_COMMENTS, boxList);
                } else {
                    if (topBoxes.size() <= bottomBoxes.size()) {
                        boxList = topBoxes;
                    } else {
                        boxList = bottomBoxes;
                    }
                }
            }
        } else {
            // determine the position to use, favoring the bottom position
            List<LNode> bottomBoxes = realNode.getProperty(InternalProperties.BOTTOM_COMMENTS);
            if (bottomBoxes == null) {
                boxList = new ArrayList<>();
                realNode.setProperty(InternalProperties.BOTTOM_COMMENTS, boxList);
            } else if (onlyBottom) {
                boxList = bottomBoxes;
            } else {
                List<LNode> topBoxes = realNode.getProperty(InternalProperties.TOP_COMMENTS);
                if (topBoxes == null) {
                    boxList = new ArrayList<>();
                    realNode.setProperty(InternalProperties.TOP_COMMENTS, boxList);
                } else {
                    if (bottomBoxes.size() <= topBoxes.size()) {
                        boxList = bottomBoxes;
                    } else {
                        boxList = topBoxes;
                    }
                }
            }
        }
        
        // add the comment box to one of the two possible lists
        boxList.add(box);
        
        // set the opposite port as property for the comment box
        box.setProperty(InternalProperties.COMMENT_CONN_PORT, oppositePort);
        // detach the edge and the opposite port
        if (edge.getTarget() == oppositePort) {
            edge.setTarget(null);
            if (oppositePort.getDegree() == 0) {
                oppositePort.setNode(null);
            }
            removeHierarchicalPortDummyNode(oppositePort);
        } else {
            edge.setSource(null);
            if (oppositePort.getDegree() == 0) {
                oppositePort.setNode(null);
            }
        }
        edge.getBendPoints().clear();
    }

    private void removeHierarchicalPortDummyNode(final LPort oppositePort) {
        LNode dummy = oppositePort.getProperty(InternalProperties.PORT_DUMMY);
        if (dummy != null) {
            Layer layer = dummy.getLayer();
            layer.getNodes().remove(dummy);
            if (layer.getNodes().isEmpty()) {
                dummy.getGraph().getLayers().remove(layer);
            }
        }
    }
}
