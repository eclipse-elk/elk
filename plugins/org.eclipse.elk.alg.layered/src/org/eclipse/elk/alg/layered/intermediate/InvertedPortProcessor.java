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
package org.eclipse.elk.alg.layered.intermediate;

import java.util.List;
import java.util.ListIterator;
import java.util.ArrayList;

import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LGraphUtil;
import org.eclipse.elk.alg.layered.graph.LLabel;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LNode.NodeType;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.layered.options.PortType;
import org.eclipse.elk.core.alg.ILayoutProcessor;
import org.eclipse.elk.core.options.EdgeLabelPlacement;
import org.eclipse.elk.core.options.PortConstraints;
import org.eclipse.elk.core.options.PortSide;
import org.eclipse.elk.core.util.IElkProgressMonitor;


/**
 * Inserts dummy nodes to cope with inverted ports.
 * 
 * <p>
 * The problem is with edges coming from the left of a node being connected to a port that's on its
 * right side, or the other way around. Let a node of that kind be in layer {@code i}. This
 * processor now takes the offending edge and connects it to a new dummy node, also in layer
 * {@code i}. Finally, the dummy is connected with the offending port. This means that once one of
 * these cases occurs in the graph, the layering is not proper anymore.
 * </p>
 * 
 * <p>
 * The dummy nodes are decorated with a
 * {@link org.eclipse.elk.alg.layered.options.LayeredOptions#NODE_TYPE} property. They are treated
 * just like ordinary {@link org.eclipse.elk.alg.layered.options.NodeType#LONG_EDGE} dummy
 * nodes
 * </p>
 * 
 * <p>
 * This processor supports self-loops by not doing anything about them. That is, no dummy nodes are
 * created for edges whose source and target node are identical.
 * </p>
 * 
 * <p>
 * Note: the following phases must support in-layer connections for this to work.
 * </p>
 * 
 * <dl>
 *   <dt>Precondition:</dt>
 *     <dd>a layered graph</dd>
 *     <dd>nodes have fixed port sides.</dd>
 *   <dt>Postcondition:</dt>
 *     <dd>dummy nodes have been inserted for edges connected to ports on odd sides
 *     <dd>the graph may contain new in-layer connections.</dd>
 *   <dt>Slots:</dt>
 *     <dd>Before phase 3.</dd>
 *   <dt>Same-slot dependencies:</dt>
 *     <dd>{@link PortSideProcessor}</dd>
 * </dl>
 * 
 * @see PortSideProcessor
 */
public final class InvertedPortProcessor implements ILayoutProcessor<LGraph> {

    @Override
    public void process(final LGraph layeredGraph, final IElkProgressMonitor monitor) {
        monitor.begin("Inverted port preprocessing", 1);
        
        // Retrieve the layers in the graph
        List<Layer> layers = layeredGraph.getLayers();
        
        // Iterate through the layers and for each layer create a list of dummy nodes
        // that were created, but not yet assigned to the layer (to avoid concurrent
        // modification exceptions)
        ListIterator<Layer> layerIterator = layers.listIterator();
        Layer currentLayer = null;
        List<LNode> unassignedNodes = new ArrayList<>();
        
        while (layerIterator.hasNext()) {
            // Update previous and current layers
            Layer previousLayer = currentLayer;
            currentLayer = layerIterator.next();
            
            // If the previous layer had unassigned nodes, assign them now and clear the list
            for (LNode node : unassignedNodes) {
                node.setLayer(previousLayer);
            }
            unassignedNodes.clear();
            
            // Iterate through the layer's nodes
            for (LNode node : currentLayer) {
                // Skip dummy nodes
                if (node.getType() != NodeType.NORMAL) {
                    continue;
                }
                
                // Skip nodes whose port sides are not fixed (because in that case, the odd
                // port side problem won't appear)
                if (!node.getProperty(LayeredOptions.PORT_CONSTRAINTS).isSideFixed()) {
                    continue;
                }
                
                // Look for input ports on the right side
                for (LPort port : node.getPorts(PortType.INPUT, PortSide.EAST)) {
                    // For every edge going into this port, insert dummy nodes (do this using
                    // a copy of the current list of edges, since the edges are modified when
                    // dummy nodes are created)
                    List<LEdge> edges = port.getIncomingEdges();
                    LEdge[] edgeArray = LGraphUtil.toEdgeArray(edges);
                    
                    for (LEdge edge : edgeArray) {
                        createEastPortSideDummies(layeredGraph, port, edge, unassignedNodes);
                    }
                }
                
                // Look for ports on the left side connected to edges going to higher layers
                for (LPort port : node.getPorts(PortType.OUTPUT, PortSide.WEST)) {
                    // For every edge going out of this port, insert dummy nodes (do this using
                    // a copy of the current list of edges, since the edges are modified when
                    // dummy nodes are created)
                    List<LEdge> edges = port.getOutgoingEdges();
                    LEdge[] edgeArray = LGraphUtil.toEdgeArray(edges);
                    
                    for (LEdge edge : edgeArray) {
                        createWestPortSideDummies(layeredGraph, port, edge, unassignedNodes);
                    }
                }
            }
        }
        
        // There may be unassigned nodes left
        for (LNode node : unassignedNodes) {
            node.setLayer(currentLayer);
        }
        
        monitor.done();
    }

    /**
     * Creates the necessary dummy nodes for an input port on the east side of a node, provided that
     * the edge connects two different nodes.
     * 
     * @param layeredGraph
     *            the layered graph.
     * @param eastwardPort
     *            the offending port.
     * @param edge
     *            the edge connected to the port.
     * @param layerNodeList
     *            list of unassigned nodes belonging to the layer of the node the port belongs to.
     *            The new dummy node is added to this list and must be assigned to the layer later.
     */
    private void createEastPortSideDummies(final LGraph layeredGraph, final LPort eastwardPort,
            final LEdge edge, final List<LNode> layerNodeList) {
        
        assert edge.getTarget() == eastwardPort;
        
        // Ignore self loops
        if (edge.getSource().getNode() == eastwardPort.getNode()) {
            return;
        }
        
        // Dummy node in the same layer
        LNode dummy = new LNode(layeredGraph);
        dummy.setType(NodeType.LONG_EDGE);
        dummy.setProperty(InternalProperties.ORIGIN, edge);
        dummy.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_POS);
        layerNodeList.add(dummy);
        
        LPort dummyInput = new LPort();
        dummyInput.setNode(dummy);
        dummyInput.setSide(PortSide.WEST);
        
        LPort dummyOutput = new LPort();
        dummyOutput.setNode(dummy);
        dummyOutput.setSide(PortSide.EAST);
        
        // Reroute the original edge
        edge.setTarget(dummyInput);
        
        // Connect the dummy with the original port
        LEdge dummyEdge = new LEdge();
        dummyEdge.copyProperties(edge);
        dummyEdge.setProperty(LayeredOptions.JUNCTION_POINTS, null);
        dummyEdge.setSource(dummyOutput);
        dummyEdge.setTarget(eastwardPort);
        
        // Set LONG_EDGE_SOURCE and LONG_EDGE_TARGET properties on the LONG_EDGE dummy
        setLongEdgeSourceAndTarget(dummy, dummyInput, dummyOutput, eastwardPort);
        
        // Move head labels from the old edge over to the new one
        ListIterator<LLabel> labelIterator = edge.getLabels().listIterator();
        while (labelIterator.hasNext()) {
            LLabel label = labelIterator.next();
            EdgeLabelPlacement labelPlacement = label.getProperty(LayeredOptions.EDGE_LABELS_PLACEMENT);
            
            if (labelPlacement == EdgeLabelPlacement.HEAD) {
                // Remember which edge the label originally belonged to, unless it already knows
                if (!label.hasProperty(InternalProperties.END_LABEL_EDGE)) {
                    label.setProperty(InternalProperties.END_LABEL_EDGE, edge);
                }
                
                labelIterator.remove();
                dummyEdge.getLabels().add(label);
            }
        }
    }

    /**
     * Creates the necessary dummy nodes for an output port on the west side of a node, provided
     * that the edge connects two different nodes.
     * 
     * @param layeredGraph
     *            the layered graph
     * @param westwardPort
     *            the offending port.
     * @param edge
     *            the edge connected to the port.
     * @param layerNodeList
     *            list of unassigned nodes belonging to the layer of the node the port belongs to.
     *            The new dummy node is added to this list and must be assigned to the layer later.
     */
    private void createWestPortSideDummies(final LGraph layeredGraph, final LPort westwardPort,
            final LEdge edge, final List<LNode> layerNodeList) {
        
        assert edge.getSource() == westwardPort;
        
        // Ignore self loops
        if (edge.getTarget().getNode() == westwardPort.getNode()) {
            return;
        }
        
        // Dummy node in the same layer
        LNode dummy = new LNode(layeredGraph);
        dummy.setType(NodeType.LONG_EDGE);
        dummy.setProperty(InternalProperties.ORIGIN, edge);
        dummy.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_POS);
        layerNodeList.add(dummy);
        
        LPort dummyInput = new LPort();
        dummyInput.setNode(dummy);
        dummyInput.setSide(PortSide.WEST);
        
        LPort dummyOutput = new LPort();
        dummyOutput.setNode(dummy);
        dummyOutput.setSide(PortSide.EAST);
        
        // Reroute the original edge
        LPort originalTarget = edge.getTarget();
        edge.setTarget(dummyInput);
        
        // Connect the dummy with the original port
        LEdge dummyEdge = new LEdge();
        dummyEdge.copyProperties(edge);
        dummyEdge.setProperty(LayeredOptions.JUNCTION_POINTS, null);
        dummyEdge.setSource(dummyOutput);
        dummyEdge.setTarget(originalTarget);
        
        // Move any head labels over to the new dummy edge
        ListIterator<LLabel> labelIterator = edge.getLabels().listIterator();
        while (labelIterator.hasNext()) {
            LLabel label = labelIterator.next();
            
            if (label.getProperty(LayeredOptions.EDGE_LABELS_PLACEMENT) == EdgeLabelPlacement.HEAD) {
                // Remember which edge the label belonged to originally
                assert !label.hasProperty(InternalProperties.END_LABEL_EDGE);
                label.setProperty(InternalProperties.END_LABEL_EDGE, edge);
                
                labelIterator.remove();
                dummyEdge.getLabels().add(label);
            }
        }
        
        // Set LONG_EDGE_SOURCE and LONG_EDGE_TARGET properties on the LONG_EDGE dummy
        setLongEdgeSourceAndTarget(dummy, dummyInput, dummyOutput, westwardPort);
    }
    
    /**
     * Properly sets the
     * {@link org.eclipse.elk.alg.layered.options.LayeredOptions#LONG_EDGE_SOURCE} and
     * {@link org.eclipse.elk.alg.layered.options.LayeredOptions#LONG_EDGE_TARGET} properties for
     * the given long edge dummy. This is required for the
     * {@link org.eclipse.elk.alg.layered.intermediate.HyperedgeDummyMerger} to work correctly.
     * 
     * @param longEdgeDummy
     *            the long edge dummy whose properties to set.
     * @param dummyInputPort
     *            the dummy node's input port.
     * @param dummyOutputPort
     *            the dummy node's output port.
     * @param oddPort
     *            the odd port that prompted the dummy to be created.
     */
    private void setLongEdgeSourceAndTarget(final LNode longEdgeDummy, final LPort dummyInputPort,
            final LPort dummyOutputPort, final LPort oddPort) {
        
        // There's exactly one edge connected to the input and output port
        LPort sourcePort = dummyInputPort.getIncomingEdges().get(0).getSource();
        LNode sourceNode = sourcePort.getNode();
        NodeType sourceNodeType = sourceNode.getType();
        
        LPort targetPort = dummyOutputPort.getOutgoingEdges().get(0).getTarget();
        LNode targetNode = targetPort.getNode();
        NodeType targetNodeType = targetNode.getType();
        
        // Set the LONG_EDGE_SOURCE property
        if (sourceNodeType == NodeType.LONG_EDGE) {
            // The source is a LONG_EDGE node; use its LONG_EDGE_SOURCE
            longEdgeDummy.setProperty(InternalProperties.LONG_EDGE_SOURCE,
                    sourceNode.getProperty(InternalProperties.LONG_EDGE_SOURCE));
        } else {
            // The target is the original node; use it
            longEdgeDummy.setProperty(InternalProperties.LONG_EDGE_SOURCE, sourcePort);
        }

        // Set the LONG_EDGE_TARGET property
        if (targetNodeType == NodeType.LONG_EDGE) {
            // The target is a LONG_EDGE node; use its LONG_EDGE_TARGET
            longEdgeDummy.setProperty(InternalProperties.LONG_EDGE_TARGET,
                    targetNode.getProperty(InternalProperties.LONG_EDGE_TARGET));
        } else {
            // The target is the original node; use it
            longEdgeDummy.setProperty(InternalProperties.LONG_EDGE_TARGET, targetPort);
        }
    }

}
