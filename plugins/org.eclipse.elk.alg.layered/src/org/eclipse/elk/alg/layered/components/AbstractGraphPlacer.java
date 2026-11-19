/*******************************************************************************
 * Copyright (c) 2012, 2015 Kiel University and others.
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
package org.eclipse.elk.alg.layered.components;

import java.util.Collection;
import java.util.List;

import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LLabel;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.math.KVectorChain;

/**
 * Takes a list of layered graphs and combines them into a single graph, placing them according to some
 * algorithm.
 * 
 * @author cds
 * @author msp
 */
abstract class AbstractGraphPlacer {
    
    /**
     * Computes a proper placement for the given graphs and combines them into a single graph.
     * 
     * @param components the graphs to be combined.
     * @param target the target graph into which the others shall be combined
     */
    public abstract void combine(List<LGraph> components, LGraph target);
    

    /**
     * Move the source graphs into the destination graph using a specified offset. The order of the
     * source graphs in the collection must not depend on their hash values. Otherwise, subsequent
     * layout calls will most likely produce different results. 
     * 
     * @param destGraph the destination graph.
     * @param sourceGraph the source graphs.
     * @param offsetx x coordinate offset.
     * @param offsety y coordinate offset.
     */
    protected void moveGraphs(final LGraph destGraph, final Collection<LGraph> sourceGraphs,
            final double offsetx, final double offsety) {
        
        for (LGraph sourceGraph : sourceGraphs) {
            moveGraph(destGraph, sourceGraph, offsetx, offsety);
        }
    }

    /**
     * Move the source graph into the destination graph using a specified offset. This method should
     * only be called once per graph, since it also applies the graph's set offset in addition to the
     * one given in the methods of this argument.
     * 
     * @param destGraph the destination graph.
     * @param sourceGraph the source graph.
     * @param offsetx x coordinate offset.
     * @param offsety y coordinate offset.
     */
    protected void moveGraph(final LGraph destGraph, final LGraph sourceGraph,
            final double offsetx, final double offsety) {
        KVector graphOffset = sourceGraph.getOffset().add(offsetx, offsety);
        
        for (LNode node : sourceGraph.getLayerlessNodes()) {
            node.getPosition().add(graphOffset);
            for (LPort port : node.getPorts()) {
                for (LEdge edge : port.getOutgoingEdges()) {
                    edge.getBendPoints().offset(graphOffset);
                    KVectorChain junctionPoints = edge.getProperty(LayeredOptions.JUNCTION_POINTS);
                    if (junctionPoints != null) {
                        junctionPoints.offset(graphOffset);
                    }
                    for (LLabel label : edge.getLabels()) {
                        label.getPosition().add(graphOffset);
                    }
                }
            }
            destGraph.getLayerlessNodes().add(node);
            node.setGraph(destGraph);
        }
    }
    
    /**
     * Offsets the given graphs by a given offset without moving their nodes to another graph. The order
     * of the graphs in the collection must not depend on their hash values. Otherwise, subsequent
     * layout calls will most likely produce different results. 
     * 
     * @param graph the graph to offset.
     * @param offsetx x coordinate offset.
     * @param offsety y coordinate offset.
     */
    protected void offsetGraphs(final Collection<LGraph> graphs, final double offsetx,
            final double offsety) {
        
        for (LGraph graph : graphs) {
            offsetGraph(graph, offsetx, offsety);
        }
    }
    
    /**
     * Offsets the given graph by a given offset without moving its nodes to another graph. This
     * method can be called as many times as required on a given graph: it does not take the graph's
     * offset into account.
     * 
     * @param graph the graph to offset.
     * @param offsetx x coordinate offset.
     * @param offsety y coordinate offset.
     */
    protected void offsetGraph(final LGraph graph, final double offsetx, final double offsety) {
        KVector graphOffset = new KVector(offsetx, offsety);
        
        for (LNode node : graph.getLayerlessNodes()) {
            node.getPosition().add(graphOffset);
            for (LPort port : node.getPorts()) {
                for (LEdge edge : port.getOutgoingEdges()) {
                    edge.getBendPoints().offset(graphOffset);
                    KVectorChain junctionPoints = edge.getProperty(LayeredOptions.JUNCTION_POINTS);
                    if (junctionPoints != null) {
                        junctionPoints.offset(graphOffset);
                    }
                    for (LLabel label : edge.getLabels()) {
                        label.getPosition().add(graphOffset);
                    }
                }
            }
        }
    }

}
