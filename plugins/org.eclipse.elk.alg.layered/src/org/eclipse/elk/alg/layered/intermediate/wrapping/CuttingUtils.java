/*******************************************************************************
 * Copyright (c) 2016 Kiel University and others.
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
package org.eclipse.elk.alg.layered.intermediate.wrapping;

import java.util.List;
import java.util.ArrayList;

import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LNode.NodeType;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.intermediate.LongEdgeSplitter;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.core.options.PortConstraints;
import org.eclipse.elk.core.options.PortSide;
import org.eclipse.elk.core.util.IndividualSpacings;


/**
 * Utility methods for cut index calculations used by more than one heuristic.
 */
public final class CuttingUtils {

    private CuttingUtils() { }
    
    /**
     * Inserts two in-layer dummies and a varying number of long edge dummies for the passed edge.
     * 
     * Let {@code originalEdge = e = (u,v)}. In-layer dummies {@code il_1} and {@code il_2} are inserted into layers
     * {@code L(u)} and {@code L(v)}. Long edge dummies {@code d_1, ..., d_m} are inserted into the layers
     * {@code L(u)+1, ..., L(v)-1}.
     * 
     * The resulting chain of edges looks as follows (if {@code offsetFirstInLayerDummy} is 1).
     * 
     * <pre>
     *                               v
     *                               î
     * il_1 -> d_1 -> ... -> d_m -> il_2 
     *  î 
     *  u
     * </pre>
     * 
     * @param layeredGraph
     *            the underlying graph with a given layering.
     * @param originalEdge
     *            the edge to be split
     * @param offsetFirstInLayerDummy
     *            the dummy nodes are inserted at the end of each layer. For the first in-layer dummy this is not always
     *            desired and can be modified using this offset.
     */
    public static List<LEdge> insertDummies(final LGraph layeredGraph,
            final LEdge originalEdge, final int offsetFirstInLayerDummy) {
        
        // to visually separate the backward wrapping edges, add some additional spacing
        double edgeNodeSpacing = layeredGraph.getProperty(LayeredOptions.SPACING_EDGE_NODE);
        double additionalSpacing = layeredGraph.getProperty(LayeredOptions.WRAPPING_ADDITIONAL_EDGE_SPACING);
        IndividualSpacings is = new IndividualSpacings();
        is.setProperty(LayeredOptions.SPACING_EDGE_NODE, edgeNodeSpacing + additionalSpacing);
        
        LEdge edge = originalEdge;
        LPort targetPort = edge.getTarget();
        LNode src = edge.getSource().getNode();
        LNode tgt = edge.getTarget().getNode();
        
        int srcIndex = src.getLayer().getIndex();
        int tgtIndex = tgt.getLayer().getIndex();
        
        List<LEdge> createdEdges = new ArrayList<>();
        
        for (int i = srcIndex; i <= tgtIndex; i++) {
            
            // Create dummy node
            LNode dummyNode = new LNode(layeredGraph);
            dummyNode.setType(NodeType.LONG_EDGE);
            dummyNode.setProperty(InternalProperties.ORIGIN, edge);
            dummyNode.setProperty(LayeredOptions.PORT_CONSTRAINTS,
                    PortConstraints.FIXED_POS);
            dummyNode.setProperty(LayeredOptions.SPACING_INDIVIDUAL, is);
            
            Layer nextLayer = layeredGraph.getLayers().get(i);
            if (i == srcIndex) {
                dummyNode.setLayer(nextLayer.getNodes().size() - offsetFirstInLayerDummy, nextLayer);
            } else {
                dummyNode.setLayer(nextLayer);
            }

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
            dummyEdge.setTarget(targetPort);

            setDummyProperties(dummyNode, edge, dummyEdge);
            createdEdges.add(dummyEdge);
            edge = dummyEdge;
        }
        
        return createdEdges;
    }
    
    /**
     * Copied from {@link LongEdgeSplitter}.
     * 
     * Sets the source and target properties on the given dummy node.
     * 
     * @param dummy the dummy node.
     * @param inEdge the edge going into the dummy node.
     * @param outEdge the edge going out of the dummy node.
     */
    private static void setDummyProperties(final LNode dummy, final LEdge inEdge, final LEdge outEdge) {
        LNode inEdgeSourceNode = inEdge.getSource().getNode();
        
        if (inEdgeSourceNode.getType() == NodeType.LONG_EDGE) {
            // The incoming edge originates from a long edge dummy node, so we can
            // just copy its properties
            dummy.setProperty(InternalProperties.LONG_EDGE_SOURCE,
                    inEdgeSourceNode.getProperty(InternalProperties.LONG_EDGE_SOURCE));
            dummy.setProperty(InternalProperties.LONG_EDGE_TARGET,
                    inEdgeSourceNode.getProperty(InternalProperties.LONG_EDGE_TARGET));
        } else {
            // The source is the input edge's source port, the target is the output
            // edge's target port
            dummy.setProperty(InternalProperties.LONG_EDGE_SOURCE, inEdge.getSource());
            dummy.setProperty(InternalProperties.LONG_EDGE_TARGET, outEdge.getTarget());
        }
    }
}
