/*******************************************************************************
 * Copyright (c) 2010, 2015 Kiel University and others.
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

import java.util.ListIterator;

import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LLabel;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LNode.NodeType;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.core.alg.ILayoutProcessor;
import org.eclipse.elk.core.options.EdgeLabelPlacement;
import org.eclipse.elk.core.options.PortConstraints;
import org.eclipse.elk.core.options.PortSide;
import org.eclipse.elk.core.util.IElkProgressMonitor;

/**
 * Splits the long edges of the layered graph to obtain a proper layering. For each edge that
 * connects two nodes that are more than one layer apart from each other, create a dummy node to
 * split the edge. The resulting layering is <i>proper</i>, i.e. all edges connect only nodes from
 * subsequent layers.
 * 
 * <p>
 * The dummy nodes retain a reference to the ports the original long edge's source and target ports.
 * </p>
 * 
 * <p>
 * The actual method that splits edges is declared as a public utility method to be called from other
 * processors that need to split edges.
 * </p>
 * 
 * <dl>
 *   <dt>Precondition:</dt>
 *     <dd>a layered graph.</dd>
 *   <dt>Postcondition:</dt>
 *     <dd>the graph is properly layered.</dd>
 *   <dt>Slots:</dt>
 *     <dd>Before phase 3.</dd>
 *   <dt>Same-slot dependencies:</dt>
 *     <dd>{@link LayerConstraintProcessor}</dd>
 * </dl>
 */
public final class LongEdgeSplitter implements ILayoutProcessor<LGraph> {

    @Override
    public void process(final LGraph layeredGraph, final IElkProgressMonitor monitor) {
        monitor.begin("Edge splitting", 1);
        
        if (layeredGraph.getLayers().size() <= 2) {
            monitor.done();
            return;
        }
        
        // Iterate through the layers
        ListIterator<Layer> layerIter = layeredGraph.getLayers().listIterator();
        Layer nextLayer = layerIter.next();
        while (layerIter.hasNext()) {
            Layer layer = nextLayer;
            nextLayer = layerIter.next();
            
            // Iterate through the nodes
            for (LNode node : layer) {
                // Iterate through the outgoing edges
                for (LPort port : node.getPorts()) {
                    // Iterate through the edges
                    for (LEdge edge : port.getOutgoingEdges()) {
                        LPort targetPort = edge.getTarget();
                        Layer targetLayer = targetPort.getNode().getLayer();
                        
                        // If the edge doesn't go to the current or next layer, split it
                        if (targetLayer != layer && targetLayer != nextLayer) {
                            // If there is no next layer, something is wrong
                            assert layerIter.hasNext();
                            
                            // Split the edge
                            splitEdge(edge, createDummyNode(layeredGraph, nextLayer, edge));
                        }
                    }
                }
            }
        }
        
        monitor.done();
    }

    /**
     * Creates a long edge dummy node to split the given edge at the given layer and adds it to the
     * layer.
     * 
     * @param layeredGraph the graph.
     * @param targetLayer the layer the dummy node will be inserted into.
     * @param edgeToSplit the edge that will later be split by the dummy node. The edge will be set as
     *                    the dummy node's {@link InternalProperties#ORIGIN ORIGIN}.
     * @return the created dummy node.
     */
    private LNode createDummyNode(final LGraph layeredGraph, final Layer targetLayer, final LEdge edgeToSplit) {
        LNode dummyNode = new LNode(layeredGraph);
        
        dummyNode.setType(NodeType.LONG_EDGE);
        dummyNode.setProperty(InternalProperties.ORIGIN, edgeToSplit);
        dummyNode.setProperty(LayeredOptions.PORT_CONSTRAINTS,
                PortConstraints.FIXED_POS);
        dummyNode.setLayer(targetLayer);
        
        return dummyNode;
    }

    /**
     * Does the actual work of splitting a given edge by rerouting it to the given dummy node and
     * introducing a new edge. Two ports are added to the dummy node and long edge properties are
     * configured for it. Also, any head labels the old edge has are moved to the new edge.
     * 
     * @param edge
     *            the edge to be split.
     * @param dummyNode
     *            the dummy node to split the edge with.
     * @return the new edge.
     */
    public static LEdge splitEdge(final LEdge edge, final LNode dummyNode) {
        LPort oldEdgeTarget = edge.getTarget();
        
        // Set thickness of the edge
        double thickness = edge.getProperty(LayeredOptions.EDGE_THICKNESS);
        if (thickness < 0) {
            thickness = 0;
            edge.setProperty(LayeredOptions.EDGE_THICKNESS, thickness);
        }
        dummyNode.getSize().y = thickness;
        double portPos = Math.floor(thickness / 2);
        
        // Create dummy input and output ports
        LPort dummyInput = new LPort();
        dummyInput.setSide(PortSide.WEST);
        dummyInput.setNode(dummyNode);
        dummyInput.getPosition().y = portPos;
        
        LPort dummyOutput = new LPort();
        dummyOutput.setSide(PortSide.EAST);
        dummyOutput.setNode(dummyNode);
        dummyOutput.getPosition().y = portPos;
        
        edge.setTarget(dummyInput);
        
        // Create a dummy edge
        LEdge dummyEdge = new LEdge();
        dummyEdge.copyProperties(edge);
        dummyEdge.setProperty(LayeredOptions.JUNCTION_POINTS, null);
        dummyEdge.setSource(dummyOutput);
        dummyEdge.setTarget(oldEdgeTarget);
        
        setDummyNodeProperties(dummyNode, edge, dummyEdge);
        moveHeadLabels(edge, dummyEdge);
        
        return dummyEdge;
    }
    
    /**
     * Sets the {@link InternalProperties#LONG_EDGE_SOURCE LONG_EDGE_SOURCE} and
     * {@link InternalProperties#LONG_EDGE_TARGET LONG_EDGE_TARGET} properties on the given dummy node.
     * 
     * @param dummyNode
     *            the dummy node.
     * @param inEdge
     *            the edge going into the dummy node.
     * @param outEdge
     *            the edge going out of the dummy node.
     */
    private static void setDummyNodeProperties(final LNode dummyNode, final LEdge inEdge,
            final LEdge outEdge) {
        
        LNode inEdgeSourceNode = inEdge.getSource().getNode();
        LNode outEdgeTargetNode = outEdge.getTarget().getNode();
        
        if (inEdgeSourceNode.getType() == NodeType.LONG_EDGE) {
            // The incoming edge originates from a long edge dummy node, so we can just copy its properties
            dummyNode.setProperty(InternalProperties.LONG_EDGE_SOURCE,
                    inEdgeSourceNode.getProperty(InternalProperties.LONG_EDGE_SOURCE));
            dummyNode.setProperty(InternalProperties.LONG_EDGE_TARGET,
                    inEdgeSourceNode.getProperty(InternalProperties.LONG_EDGE_TARGET));
            dummyNode.setProperty(InternalProperties.LONG_EDGE_HAS_LABEL_DUMMIES,
                    inEdgeSourceNode.getProperty(InternalProperties.LONG_EDGE_HAS_LABEL_DUMMIES));
            
        } else if (inEdgeSourceNode.getType() == NodeType.LABEL) {
            // The incoming edge originates from a label dummy node, so we can just copy its properties
            dummyNode.setProperty(InternalProperties.LONG_EDGE_SOURCE,
                    inEdgeSourceNode.getProperty(InternalProperties.LONG_EDGE_SOURCE));
            dummyNode.setProperty(InternalProperties.LONG_EDGE_TARGET,
                    inEdgeSourceNode.getProperty(InternalProperties.LONG_EDGE_TARGET));
            dummyNode.setProperty(InternalProperties.LONG_EDGE_HAS_LABEL_DUMMIES, true);
            
        } else if (outEdgeTargetNode.getType() == NodeType.LABEL) {
            // The outgoing edge points to a label dummy node, so we can just copy its properties
            dummyNode.setProperty(InternalProperties.LONG_EDGE_SOURCE,
                    outEdgeTargetNode.getProperty(InternalProperties.LONG_EDGE_SOURCE));
            dummyNode.setProperty(InternalProperties.LONG_EDGE_TARGET,
                    outEdgeTargetNode.getProperty(InternalProperties.LONG_EDGE_TARGET));
            dummyNode.setProperty(InternalProperties.LONG_EDGE_HAS_LABEL_DUMMIES, true);
            
        } else {
            // The source is the input edge's source port, the target is the output edge's target port. Also,
            // the long edge doesn't seem to have label dummy nodes
            dummyNode.setProperty(InternalProperties.LONG_EDGE_SOURCE, inEdge.getSource());
            dummyNode.setProperty(InternalProperties.LONG_EDGE_TARGET, outEdge.getTarget());
        }
    }
    
    /**
     * Moves all head labels from a given split edge to the new edge created to split it.
     * 
     * @param oldEdge
     *            the old edge whose head labels are to be moved.
     * @param newEdge
     *            the new edge whose head labels are to be moved.
     */
    private static void moveHeadLabels(final LEdge oldEdge, final LEdge newEdge) {
        ListIterator<LLabel> labelIterator = oldEdge.getLabels().listIterator();
        while (labelIterator.hasNext()) {
            LLabel label = labelIterator.next();
            EdgeLabelPlacement labelPlacement = label.getProperty(LayeredOptions.EDGE_LABELS_PLACEMENT);
            
            if (labelPlacement == EdgeLabelPlacement.HEAD) {
                labelIterator.remove();
                newEdge.getLabels().add(label);
                
                if (!label.hasProperty(InternalProperties.END_LABEL_EDGE)) {
                    label.setProperty(InternalProperties.END_LABEL_EDGE, oldEdge);
                }
            }
        }
    }

}
