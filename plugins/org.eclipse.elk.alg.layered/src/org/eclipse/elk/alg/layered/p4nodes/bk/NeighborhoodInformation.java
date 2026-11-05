/*******************************************************************************
 * Copyright (c) 2015 Kiel University and others.
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
package org.eclipse.elk.alg.layered.p4nodes.bk;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.ArrayList;

import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.core.util.Pair;


/**
 * Class holds neighborhood information for a layered graph that is used during bk node placing.
 * Since this information is required multiple times but does not change during processing we
 * precalculate it in this class.
 * 
 * @author uru
 */
public final class NeighborhoodInformation {
    
    // Allow the fields of this container to be accessed from package siblings.
    // SUPPRESS CHECKSTYLE NEXT 24 VisibilityModifier
    
    /** Number of nodes in the graph. */
    public int nodeCount;
    /** For a layer {@code l} the entry at {@code layerIndex[l.id]} holds the index of layer {@code l}. */
    public int[] layerIndex;
    /** For a node {@code n} the entry at {@code nodeIndex[n.id]} holds the index of {@code n} within its layer. */
    public int[] nodeIndex;
    /**
     * For a node n holds leftNeighbors.get(n.id) holds a list with all left neighbors along with
     * any edge that connects n to its neighbor.
     */
    public List<List<Pair<LNode, LEdge>>> leftNeighbors;
    /** See javadoc of {@link #leftNeighbors}. */
    public List<List<Pair<LNode, LEdge>>> rightNeighbors;
    /** A comparator to sort a list of neighbor nodes. */
    private NeighborComparator neighborComparator;
    
    /**
     * Use {@link #buildFor(LGraph)} to construct an instance.
     */
    private NeighborhoodInformation() {
    }
    
    /**
     * Release allocated resources.
     */
    public void cleanup() {
        layerIndex = null;
        nodeIndex = null;
        leftNeighbors.clear();
        rightNeighbors.clear();
        neighborComparator = null;
    }
    
    /**
     * Creates and properly initializes the neighborhood information required by the
     * {@link BKNodePlacer}. This includes:
     * <ul>
     * <li>calculating the number of nodes in the graph</li>
     * <li>assigning a unique id to every layer and node</li>
     * <li>recording the index of every node in its layer</li>
     * <li>calculating left and right neighbors for every node</li>
     * </ul>
     * 
     * @param graph
     *            the underlying graph
     * @return a properly initialized instance
     */
    public static NeighborhoodInformation buildFor(final LGraph graph) {
        NeighborhoodInformation ni = new NeighborhoodInformation();
        
        ni.nodeCount = 0;
        for (Layer layer : graph) {
            ni.nodeCount += layer.getNodes().size();
        }
        
        // cache indexes of layers and of nodes
        int lId = 0;
        int lIndex = 0;
        ni.layerIndex = new int[graph.getLayers().size()];
        int nId = 0;
        int nIndex = 0;
        ni.nodeIndex = new int[ni.nodeCount];
        for (Layer l : graph.getLayers()) {
            l.id = lId++;
            ni.layerIndex[l.id] = lIndex++;
            nIndex = 0;
            for (LNode n : l.getNodes()) {
                n.id = nId++;
                ni.nodeIndex[n.id] = nIndex++;
            }
        }
        
        // we will need a comparator for the next step
        ni.neighborComparator = ni.new NeighborComparator();
        
        // determine all left and right neighbors of the graph's nodes
        ni.leftNeighbors = new ArrayList<>(ni.nodeCount);
        determineAllLeftNeighbors(ni, graph);
        ni.rightNeighbors = new ArrayList<>(ni.nodeCount);
        determineAllRightNeighbors(ni, graph);
        
        return ni;
    }
    
    
    /**
     * Give all right neighbors (originally known as lower neighbors) of a given node. A lower
     * neighbor is a node in a following layer that has an edge coming from the given node.
     */
    private static void determineAllRightNeighbors(final NeighborhoodInformation ni, final LGraph graph) {
        for (Layer l : graph) {
            for (LNode n : l) {

                List<Pair<LNode, LEdge>> result = new ArrayList<>();
                int maxPriority = 0;
                
                for (LEdge edge : n.getOutgoingEdges()) {
                    if (edge.isSelfLoop() || edge.isInLayerEdge()) {
                        continue;
                    } 
                    int edgePrio = edge.getProperty(LayeredOptions.PRIORITY_STRAIGHTNESS);
                    if (edgePrio > maxPriority) {
                        maxPriority = edgePrio;
                        result.clear();
                    }
                    if (edgePrio == maxPriority) {
                        result.add(Pair.of(edge.getTarget().getNode(), edge));
                    }
                }

                Collections.sort(result, ni.neighborComparator);

                ni.rightNeighbors.add(n.id, result);
            }
        }
    }
    

    /**
     * Gives all left neighbors (originally known as upper neighbors) of a given node. An upper
     * neighbor is a node in a previous layer that has an edge pointing to the given node.
     */
    private static void determineAllLeftNeighbors(final NeighborhoodInformation ni, final LGraph graph) {
        for (Layer l : graph) {
            for (LNode n : l) {
                
                List<Pair<LNode, LEdge>> result = new ArrayList<>();
                int maxPriority = 0;
                
                for (LEdge edge : n.getIncomingEdges()) {
                    if (edge.isSelfLoop() || edge.isInLayerEdge()) {
                        continue;
                    } 
                    int edgePrio = edge.getProperty(LayeredOptions.PRIORITY_STRAIGHTNESS);
                    if (edgePrio > maxPriority) {
                        maxPriority = edgePrio;
                        result.clear();
                    }
                    if (edgePrio == maxPriority) {
                        result.add(Pair.of(edge.getSource().getNode(), edge));
                    }
                }
                
                Collections.sort(result, ni.neighborComparator);
                
                ni.leftNeighbors.add(n.id, result);
            }
        }
    }
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    // Class NeighborComparator
    
    // SUPPRESS CHECKSTYLE NEXT 1 Javadoc
    private final class NeighborComparator implements Comparator<Pair<LNode, LEdge>> {
        @Override
        public int compare(final Pair<LNode, LEdge> o1, final Pair<LNode, LEdge> o2) {
            int cmp = nodeIndex[o1.getFirst().id] - nodeIndex[o2.getFirst().id];
            return (int) Math.signum(cmp);
        }
    }
}
