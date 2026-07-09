/*******************************************************************************
 * Copyright (c) 2014, 2019 Kiel University and others.
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
package org.eclipse.elk.alg.layered;

import java.io.StringWriter;
import java.util.Iterator;
import java.util.List;

import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LNode.NodeType;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.p4nodes.LinearSegmentsNodePlacer.LinearSegment;
import org.eclipse.elk.alg.layered.p5edges.orthogonal.HyperEdgeSegmentDependency;
import org.eclipse.elk.alg.layered.p5edges.orthogonal.HyperEdgeSegment;

/**
 * A utility class for debugging of ELK Layered.
 */
public final class DotDebugUtil {
    
    
    /**
     * Hidden constructor to avoid instantiation.
     */
    private DotDebugUtil() { }


    /**
     * Output a representation of the given graph in dot format. The following conventions are used:
     * <ul>
     * <li>Standard nodes are drawn as rectangles.</li>
     * <li>Dummy nodes are drawn as ellipses.</li>
     * <li>Nodes have a color that depends on their node type. (yellow for {@code LONG_EDGE},
     * turquoise for {@code ODD_PORT_SIDE}, dark blue for {@code NORTH_SOUTH_PORT})</li>
     * </ul>
     * 
     * @param lgraph
     *            the layered graph
     * @return debug graph in DOT format.
     */
    public static String createDebugGraph(final LGraph lgraph) {
        StringWriter writer = new StringWriter();
                    
        // Begin the digraph
        writer.write("digraph {\n");
        
        // Digraph options
        writer.write("    rankdir=LR;\n");
        
        // Write layerless nodes and edges
        writeLayer(writer, -1, lgraph.getLayerlessNodes());
        
        // Go through the layers
        int layerNumber = -1;
        for (Layer layer : lgraph) {
            layerNumber++;
            
            // Write the nodes and edges
            writeLayer(writer, layerNumber, layer.getNodes());
        }
        
        // Close the digraph. And the writer.
        writer.write("}\n");
        
        return writer.toString();
    }
    
    /**
     * Writes a debug graph for the given linear segments and their dependencies.
     * 
     * @param layeredGraph
     *            the layered graph.
     * @param segmentList
     *            the list of linear segments.
     * @param outgoingList
     *            the list of successors for each linear segment.
     * @return segments graph in DOT format.
     */
    public static String createDebugGraph(final LGraph layeredGraph, final List<LinearSegment> segmentList,
            final List<List<LinearSegment>> outgoingList) {

        StringWriter writer = new StringWriter();
        
        writer.write("digraph {\n");

        Iterator<LinearSegment> segmentIterator = segmentList.iterator();
        Iterator<List<LinearSegment>> successorsIterator = outgoingList.iterator();

        while (segmentIterator.hasNext()) {
            LinearSegment segment = segmentIterator.next();
            List<LinearSegment> successors = successorsIterator.next();

            writer.write("  " + segment.hashCode() + "[label=\"" + segment + "\"]\n");

            for (LinearSegment successor : successors) {
                writer.write("  " + segment.hashCode() + "->" + successor.hashCode() + "\n");
            }
        }

        writer.write("}\n");
        
        return writer.toString();
    }
    
    /**
     * Writes a debug graph for the given list of hypernodes.
     * 
     * @param layeredGraph the layered graph
     * @param hypernodes a list of hypernodes
     * @return hypernodes graph in DOT format.
     */
    public static String createDebugGraph(final LGraph layeredGraph, final List<HyperEdgeSegment> hypernodes) {
        StringWriter writer = new StringWriter();
        
        writer.write("digraph {\n");
        
        // Write hypernode information
        for (HyperEdgeSegment hypernode : hypernodes) {
            writer.write("  " + hypernode.hashCode() + "[label=\""
                    + hypernode.toString() + "\"]\n");
        }
        
        // Write dependency information
        for (HyperEdgeSegment hypernode : hypernodes) {
            for (HyperEdgeSegmentDependency dependency : hypernode.getOutgoingSegmentDependencies()) {
                writer.write("  " + hypernode.hashCode() + "->" + dependency.getTarget().hashCode()
                        + "[label=\"" + dependency.getType().name() + " (" + dependency.getWeight() + ")\"]\n");
            }
        }
        
        writer.write("}\n");

        return writer.toString();
    }
    
    /**
     * Writes the given list of nodes and their edges.
     * 
     * @param writer writer to write to.
     * @param layerNumber the layer number. {@code -1} for layerless nodes.
     * @param nodes the nodes in the layer.
     */
    private static void writeLayer(final StringWriter writer, final int layerNumber, final List<LNode> nodes) {
        if (nodes.isEmpty()) {
            return;
        }
        
        // Go through the layer's nodes
        int nodeNumber = -1;
        for (LNode node : nodes) {
            nodeNumber++;
            
            // The node's name in the output is its hash code (unique!)
            writer.write("        " + node.hashCode());
            
            // Options time!
            StringBuffer options = new StringBuffer();
            
            // Label
            options.append("label=\"");
            if (node.getType() == NodeType.NORMAL) {
                // Normal nodes display their name, if any
                if (node.getDesignation() != null) {
                    options.append(node.getDesignation().replace("\"", "\\\"") + " ");
                }
            } else {
                // Dummy nodes show their name (if set), or their node ID
                if (node.getDesignation() != null) {
                    options.append(node.getDesignation().replace("\"", "\\\"") + " ");
                } else {
                    options.append("n_" + node.id + " ");
                }
                if (node.getType() == NodeType.NORTH_SOUTH_PORT) {
                    Object origin = node.getProperty(InternalProperties.ORIGIN);
                    if (origin instanceof LNode) {
                        options.append("(" + ((LNode) origin).toString() + ")");
                    }
                }
            }
            options.append("(" + layerNumber + "," + nodeNumber + ")\",");
            
            // Node type
            if (node.getType() == NodeType.NORMAL) {
                options.append("shape=box,");
            } else {
                options.append("style=\"rounded,filled\",");
                
                String color = node.getType().getColor();
                if (color != null) {
                    options.append("color=\"" + color + "\",");
                }
            }
            
            // Print options, if any
            options.deleteCharAt(options.length() - 1);
            if (options.length() > 0) {
                writer.write("[" + options + "]");
            }
            
            // End the node line
            writer.write(";\n");
        }
        
        // Write the edges
        for (LNode node : nodes) {
            // Go through all edges and output those that have this node as their source
            for (LPort port : node.getPorts()) {
                for (LEdge edge : port.getOutgoingEdges()) {
                    writer.write("    " + node.hashCode() + " -> "
                            + edge.getTarget().getNode().hashCode());
                    writer.write(";\n");
                }
            }
        }
    }

}
