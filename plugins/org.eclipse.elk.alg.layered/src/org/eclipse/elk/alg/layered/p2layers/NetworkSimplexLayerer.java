/*******************************************************************************
 * Copyright (c) 2010, 2020 Kiel University and others.
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
package org.eclipse.elk.alg.layered.p2layers;

import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.ArrayList;
import java.util.HashMap;

import org.eclipse.elk.alg.common.networksimplex.NEdge;
import org.eclipse.elk.alg.common.networksimplex.NGraph;
import org.eclipse.elk.alg.common.networksimplex.NNode;
import org.eclipse.elk.alg.common.networksimplex.NetworkSimplex;
import org.eclipse.elk.alg.layered.LayeredPhases;
import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.intermediate.IntermediateProcessorStrategy;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.core.alg.ILayoutPhase;
import org.eclipse.elk.core.alg.LayoutProcessorConfiguration;
import org.eclipse.elk.core.util.IElkProgressMonitor;


/**
 * The main class of the network simplex layerer component. It offers an algorithm to determine an
 * optimal layering of all nodes in the graph concerning a minimal length of all edges using the
 * network simplex algorithm described in
 * <ul>
 * <li>Emden R. Gansner, Eleftherios Koutsofios, Stephen C. North, Kiem-Phong Vo, A technique for
 * drawing directed graphs. <i>Software Engineering</i> 19(3), pp. 214-230, 1993.</li>
 * </ul>
 * 
 * <dl>
 * <dt>Precondition:</dt>
 * <dd>the graph has no cycles</dd>
 * <dt>Postcondition:</dt>
 * <dd>all nodes have been assigned a layer such that edges connect only nodes from layers with
 * increasing indices</dd>
 * </dl>
 */
public final class NetworkSimplexLayerer implements ILayoutPhase<LayeredPhases, LGraph> {

    /** intermediate processing configuration. */
    private static final LayoutProcessorConfiguration<LayeredPhases, LGraph> BASELINE_PROCESSING_CONFIGURATION =
        LayoutProcessorConfiguration.<LayeredPhases, LGraph>create()
            .addBefore(LayeredPhases.P1_CYCLE_BREAKING,
                    IntermediateProcessorStrategy.EDGE_AND_LAYER_CONSTRAINT_EDGE_REVERSER)
            .addBefore(LayeredPhases.P2_LAYERING, IntermediateProcessorStrategy.LAYER_CONSTRAINT_PREPROCESSOR)
            .addBefore(LayeredPhases.P3_NODE_ORDERING, IntermediateProcessorStrategy.LAYER_CONSTRAINT_POSTPROCESSOR);


    // ================================== Attributes ==============================================

    /** The layered graph all methods in this class operate on. */
    private LGraph layeredGraph;

    /**
     * A {@code LinkedList} containing all nodes of the currently identified connected component by
     * {@code connectedComponents()}.
     * 
     * @see #connectedComponents(List)
     */
    private List<LNode> componentNodes;

    /**
     * A flag indicating whether a specified node has been visited during DFS-traversal. This array
     * has to be filled with {@code false} each time, before a DFS-based method is invoked.
     */
    private boolean[] nodeVisited;
    
    
    // =============================== Initialization Methods =====================================

    @Override
    public LayoutProcessorConfiguration<LayeredPhases, LGraph> getLayoutProcessorConfiguration(final LGraph graph) {
        return BASELINE_PROCESSING_CONFIGURATION;
    }

    /**
     * Helper method for the network simplex layerer. It determines all connected components of the
     * graph given by a {@code Collection} containing all nodes of the graph.
     * 
     * @param theNodes
     *            a {@code Collection} containing all nodes of the graph to determine the connected
     *            components
     * @return a {@code LinkedList} of {@code LinkedLists} containing all nodes of every connected
     *         component
     * 
     * @see NetworkSimplexLayerer#connectedComponentsDFS(LNode) connectedComponentsDFS()
     */
    private List<List<LNode>> connectedComponents(final List<LNode> theNodes) {
        // initialize required attributes
        if (nodeVisited == null || nodeVisited.length < theNodes.size()) {
            nodeVisited = new boolean[theNodes.size()];
        } else {
            Arrays.fill(nodeVisited, false);
        }
        componentNodes = new ArrayList<>();

        // re-index nodes
        int counter = 0;
        for (LNode node : theNodes) {
            node.id = counter++;
        }
        // determine connected components
        LinkedList<List<LNode>> components = new LinkedList<>();
        for (LNode node : theNodes) {
            if (!nodeVisited[node.id]) {
                connectedComponentsDFS(node);
                // connected component with the most nodes should be layered first to guarantee
                // reusability of attribute instances
                if (components.isEmpty() || components.getFirst().size() < componentNodes.size()) {
                    components.addFirst(componentNodes);
                } else {
                    components.addLast(componentNodes);
                }
                componentNodes = new ArrayList<>();
            }
        }
        return components;
    }

    /**
     * Helper method for the connected components determination. It determines all nodes, that are
     * connected with the input node (i.e. all nodes of that connected component the input node is
     * part of) and adds them to {@code componentNodes}.
     * 
     * @param node
     *            the root of the DFS-subtree
     * @return a {@code LinkedList} containing all nodes reachable through a path beginning at the
     *         input node (i.e. all nodes connected to the input node)
     * 
     * @see NetworkSimplexLayerer#connectedComponents(List) connectedComponents()
     * @see NetworkSimplexLayerer#componentNodes componentNodes
     */
    private void connectedComponentsDFS(final LNode node) {
        nodeVisited[node.id] = true;
        
        // node is part of the current connected component
        componentNodes.add(node);
        LNode opposite;
        
        // continue with next nodes, if not already visited
        for (LPort port : node.getPorts()) {
            for (LEdge edge : port.getConnectedEdges()) {
                opposite = getOpposite(port, edge).getNode();
                if (!nodeVisited[opposite.id]) {
                    connectedComponentsDFS(opposite);
                }
            }
        }
    }

    /**
     * Helper method for the network simplex layerer. It instantiates all necessary attributes for
     * the execution of the network simplex layerer and initializes them with their default values.
     * All edges in the connected component given by the input argument will be determined, as well
     * as the number of incoming and outgoing edges of each node ( {@code inDegree}, respectively
     * {@code outDegree}). All sinks and source nodes in the connected component identified in this
     * step will be added to {@code sinks}, respectively {@code sources}.
     * 
     * @param theNodes
     *            a {@code Collection} containing all nodes of the graph
     */
    private NGraph initialize(final List<LNode> theNodes) {

        final Map<LNode, NNode> nodeMap = new HashMap<>();
        
        // transform nodes
        NGraph graph = new NGraph();
        for (LNode lNode : theNodes) {
            NNode nNode = NNode.of()
                               .origin(lNode)
                               .create(graph);
            nodeMap.put(lNode, nNode);
        }
        
        // transform edges
        for (LNode lNode : theNodes) {
            for (LEdge lEdge : lNode.getOutgoingEdges()) {
                
                // ignore self-loops
                if (lEdge.isSelfLoop()) {
                    continue;
                }
                
                NEdge.of(lEdge)
                     .weight(1 * Math.max(1, lEdge.getProperty(LayeredOptions.PRIORITY_SHORTNESS)))
                     .delta(1)
                     .source(nodeMap.get(lEdge.getSource().getNode()))
                     .target(nodeMap.get(lEdge.getTarget().getNode()))
                     .create();
            }
        }
        
        return graph;
    }


    /**
     * Release all created resources so the GC can reap them.
     */
    private void dispose() {
        this.componentNodes = null;
        this.layeredGraph = null;
        this.nodeVisited = null;
    }

    // ============================== Network-Simplex Algorithm ===================================

    /** factor by which the maximal number of iterations is multiplied. */
    private static final int ITER_LIMIT_FACTOR = 4;

    /**
     * The main method of the network simplex layerer. It determines an optimal layering of all
     * nodes in the graph concerning a minimal length of all edges by using the network simplex
     * algorithm described in {@literal Emden R. Gansner, Eleftherios Koutsofios, Stephen
     * C. North, Kiem-Phong Vo: "A Technique for Drawing Directed Graphs", AT&T Bell Laboratories.
     * Note that the execution time of this implemented algorithm has not been proven quadratic yet.
     * 
     * @param theLayeredGraph
     *            a layered graph which initially only contains layerless nodes and is
     *            then filled with layers
     * @param monitor
     *            the progress monitor
     */
    public void process(final LGraph theLayeredGraph, final IElkProgressMonitor monitor) {
        monitor.begin("Network simplex layering", 1);

        layeredGraph = theLayeredGraph;
        int thoroughness = theLayeredGraph.getProperty(LayeredOptions.THOROUGHNESS) * ITER_LIMIT_FACTOR;

        List<LNode> theNodes = layeredGraph.getLayerlessNodes();
        if (theNodes.size() < 1) {
            monitor.done();
            return;
        }

        // layer graph, each connected component separately
        List<List<LNode>> connectedComponents = connectedComponents(theNodes);
        int[] previousLayeringNodeCounts = null;
        for (List<LNode> connComp : connectedComponents) {
            
            // determine a limit on the number of iterations
            int iterLimit = thoroughness * (int) Math.sqrt(connComp.size());

            NGraph graph = initialize(connComp);

            // execute the network simplex algorithm on the (sub-)graph
            NetworkSimplex.forGraph(graph).withIterationLimit(iterLimit)
                    .withPreviousLayering(previousLayeringNodeCounts)
                    .withBalancing(true)
                    .execute(monitor.subTask(1));

            // the layers are store in the NNode's layer field.
            List<Layer> layers = layeredGraph.getLayers();
            for (NNode nNode : graph.nodes) {
                // add additional layers to match required number
                while (layers.size() <= nNode.layer) {
                    layers.add(layers.size(), new Layer(layeredGraph));
                }
                LNode lNode = (LNode) nNode.origin;
                lNode.setLayer(layers.get(nNode.layer));
            }
            
            if (connectedComponents.size() > 1) {
                previousLayeringNodeCounts = new int[layeredGraph.getLayers().size()];
                int layerIdx = 0;
                for (Layer l : layeredGraph) {
                    previousLayeringNodeCounts[layerIdx++] = l.getNodes().size();
                }
            }
        }

        // empty the list of unlayered nodes
        theNodes.clear();

        // release the created resources
        dispose();
        monitor.done();
    }

    /**
     * Helper method for the network simplex layerer. It returns the port that is connected to the
     * opposite side of the specified edge from the viewpoint of the input port.
     * 
     * @param port
     *            the port to get the opposite port from
     * @param edge
     *            the edge to consider when determining the opposite port
     * @return the opposite port from the viewpoint of the given port
     * 
     * @throws IllegalArgumentException
     *             if the input edge is not connected to the input port
     */
    private LPort getOpposite(final LPort port, final LEdge edge) {
        if (edge.getSource().equals(port)) {
            return edge.getTarget();
        } else if (edge.getTarget().equals(port)) {
            return edge.getSource();
        }
        throw new IllegalArgumentException("Input edge is not connected to the input port.");
    }

}
