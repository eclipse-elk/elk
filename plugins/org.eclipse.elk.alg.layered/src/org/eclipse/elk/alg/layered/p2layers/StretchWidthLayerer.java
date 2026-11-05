/*******************************************************************************
 * Copyright (c) 2016, 2020 Kiel University and others.
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
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Set;
import java.util.ArrayList;
import java.util.HashSet;

import org.eclipse.elk.alg.layered.LayeredPhases;
import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LNode.NodeType;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.intermediate.IntermediateProcessorStrategy;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.core.alg.ILayoutPhase;
import org.eclipse.elk.core.alg.LayoutProcessorConfiguration;
import org.eclipse.elk.core.util.IElkProgressMonitor;


/**
 * This class implements the StretchWidth layering algorithm described in In Search for Efficient
 * Heuristics for Minimum-Width Graph Layering with Consideration of Dummy Nodes written by Nikolas
 * S. Nikolov, Alexandre Tarassov, and Jürgen Branke It is designed to create a layering as small as
 * possible.
 *
 * Precondition: the graph has no cycles, but may contain self-loops Postcondition: all nodes have
 * been assigned to a layer such that edges connect only nodes from one layer to another layer
 */
public class StretchWidthLayerer implements ILayoutPhase<LayeredPhases, LGraph> {
    /** Indicates the width of the currently built layer. */
    private double widthCurrent = 0;
    /** Estimated width of the next layer. */
    private double widthUp = 0;
    /** Indicates the width of the layering. */
    // Will be dynamically increased if the layering fails for lower values
    // The pseudo code says initialize with 0 but I think you can start with 1 or maybe
    // higher/average out-degree
    // I choose,later, the widest normalized node, as the algorithm fails until it can be placed
    private double maxWidth = 0;
    /** The graph the layering is done for. */
    private LGraph currentGraph;
    /** Factor which determines the bounding of the incoming edges with widthUp and maxWidth. */
    private double upperLayerInfluence;
    /** Sorted list of layerless nodes. */
    private List<LNode> sortedLayerlessNodes;
    /** Set of nodes placed in the current layer. */
    private Set<LNode> alreadyPlacedNodes = new HashSet<>();
    /** Set of nodes in all layers except the current. */
    private Set<LNode> alreadyPlacedInOtherLayers = new HashSet<>();
    /**
     * List of layerless nodes to be used in one layering approach, will be sorted after
     * initialization.
     */
    private List<LNode> tempLayerlessNodes;
    /**
     * List of sets of successors for every node, can be accessed with successors.get(node.id).
     */
    private List<Set<LNode>> successors;
    /**
     * Array of out-degrees for every node, access with outDegree[node.id].
     */
    private int[] outDegree;
    /**
     * Initialized as a copy of {@link #outDegree}, used to determine if all successors of a node
     * are already placed.
     */
    private int[] remainingOutGoing;
    /**
     * Array of in-degrees for every node, access with inDegree[node.id].
     */
    private int[] inDegree;
    /** Selected node to be placed. */
    private LNode selectedNode;
    /** Minimum node size to normalize the other sizes. */
    private double minimumNodeSize;
    /** Maximum node size to normalize the other sizes. */
    private double maximumNodeSize;
    /**
     * Normalized Size of every node, access with normSize[node.id].
     */
    private double[] normSize;
    /** the size of the dummy nodes, that should be considered. */
    private double dummySize;

    @Override
    public LayoutProcessorConfiguration<LayeredPhases, LGraph> getLayoutProcessorConfiguration(final LGraph graph) {
        return LayoutProcessorConfiguration.<LayeredPhases, LGraph>create()
                .addBefore(LayeredPhases.P1_CYCLE_BREAKING,
                        IntermediateProcessorStrategy.EDGE_AND_LAYER_CONSTRAINT_EDGE_REVERSER)
                .addBefore(LayeredPhases.P2_LAYERING,
                        IntermediateProcessorStrategy.LAYER_CONSTRAINT_PREPROCESSOR)
                .addBefore(LayeredPhases.P3_NODE_ORDERING,
                        IntermediateProcessorStrategy.LAYER_CONSTRAINT_POSTPROCESSOR);
    }

    @Override
    public void process(final LGraph layeredGraph, final IElkProgressMonitor progressMonitor) {
        progressMonitor.begin("StretchWidth layering", 1);

        // if no nodes need to be placed we can stop right here
        if (layeredGraph.getLayerlessNodes().isEmpty()) {
            progressMonitor.done();
            return;
        }

        // set graph
        currentGraph = layeredGraph;
        // reset variables
        widthCurrent = 0;
        widthUp = 0;
        // initialize min and max node sizes to 1,
        // since 0 or negative values to not make sense
        minimumNodeSize = Double.POSITIVE_INFINITY;
        maximumNodeSize = Double.NEGATIVE_INFINITY;

        // initialize the dummy size with the spacing properties
        dummySize = layeredGraph.getProperty(LayeredOptions.SPACING_EDGE_EDGE).doubleValue();
        // Sort the nodes at beginning, since the rank will not change.
        // The list is sorted in descending order by the rank.
        // It has to be computed first
        computeSortedNodes();
        // Compute a list with the successors of every node
        // Can't be swapped with the computation of the out- and in-degrees
        computeSuccessors();
        // Compute two arrays of out- and in-degree for every node
        computeDegrees();
        // Compute minimum node size
        minMaxNodeSize();
        // Compute normalized size of every node
        computeNormalizedSize();
        // make sure the values are reasonable
        minimumNodeSize = Math.max(1, minimumNodeSize);
        maximumNodeSize = Math.max(1, maximumNodeSize);

        // normalize dummy size
        dummySize = dummySize / minimumNodeSize;
        maxWidth = maximumNodeSize / minimumNodeSize;
        // To combine dummy nodes and real node size, you could also add the average normalized node
        // size
        // to this criterion and adjust the condition go up at widthUp
        upperLayerInfluence = getAverageOutDegree();
        // Layer currently worked on
        Layer currentLayer = new Layer(currentGraph);
        // add first layer to the graph
        currentGraph.getLayers().add(currentLayer);
        // Copy the sorted layerless nodes so we don't overwrite it in the reset case
        tempLayerlessNodes = new ArrayList<>(sortedLayerlessNodes);
        // Copy the outDegree Array
        remainingOutGoing = Arrays.copyOf(outDegree, outDegree.length);

        while (!tempLayerlessNodes.isEmpty()) {
            // Select a node to be placed
            selectedNode = selectNode();
            // The pseudo code computes the difference of the to sets as:
            // alreadyPlacedNodes\alreadyPlacedInOtherLayers
            // but since alreadyPlacedNodes will be cleared here, it's enough to check if
            // alreadyPlacedNodes is empty
            // if it is empty it would indicate an empty layer
            if (selectedNode == null || (conditionGoUp() && !alreadyPlacedNodes.isEmpty())) {
                // go to the next layer //
                // update the remaining successors of the node
                updateOutGoing(currentLayer);
                currentLayer = new Layer(currentGraph);
                currentGraph.getLayers().add(currentLayer);
                // union of alreadyPlacedInOtherLayers and alreadyPlacedNodes
                alreadyPlacedInOtherLayers.addAll(alreadyPlacedNodes);
                alreadyPlacedNodes.clear();
                // change width
                widthCurrent = widthUp;
                widthUp = 0;
            } else {
                if (conditionGoUp()) {
                    // reset layering //
                    // clear the placed layers
                    currentGraph.getLayers().clear();
                    // create the new first layer;
                    currentLayer = new Layer(currentGraph);
                    currentGraph.getLayers().add(currentLayer);
                    // reset variables
                    widthCurrent = 0;
                    widthUp = 0;
                    alreadyPlacedNodes.clear();
                    alreadyPlacedInOtherLayers.clear();
                    // increase maxWidth
                    maxWidth++;
                    // reset layerless nodes
                    tempLayerlessNodes = new ArrayList<>(sortedLayerlessNodes);
                    // reset successors
                    remainingOutGoing = Arrays.copyOf(outDegree, outDegree.length);

                } else {
                    // add node to current layer //
                    selectedNode.setLayer(currentLayer);
                    tempLayerlessNodes.remove(selectedNode);
                    alreadyPlacedNodes.add(selectedNode);
                    // compute new widthCurrent and widthUp
                    widthCurrent = widthCurrent - outDegree[selectedNode.id] * dummySize + normSize[selectedNode.id];
                    widthUp += inDegree[selectedNode.id] * dummySize;
                }
            }
        }

        // Layering done, delete original layerlessNodes
        layeredGraph.getLayerlessNodes().clear();
        // Algorithm is Bottom-Up -> reverse Layers
        Collections.reverse(layeredGraph.getLayers());

        progressMonitor.done();
    }

    /**
     * Checks the effects of the hypothetical placement of the selected node and whether the
     * algorithm should rather go up than placing the node.
     *
     * @return true, if the algorithm should go to the next layer, false otherwise
     */
    private boolean conditionGoUp() {
        boolean a = ((widthCurrent - (outDegree[selectedNode.id] * dummySize) + normSize[selectedNode.id]) > maxWidth);
        boolean b = ((widthUp + inDegree[selectedNode.id] * dummySize) > (maxWidth * upperLayerInfluence * dummySize));
        return a || b;
    }

    /**
     * Selects a node from the sorted list of layerless nodes. The selection is done according to
     * the rank of the node and only if all of its successors are already in the
     * {@link #alreadyPlacedInOtherLayers} set.
     *
     * @return node to be placed in the current layer, or null if no appropriate node was found
     */
    private LNode selectNode() {
        for (LNode node : tempLayerlessNodes) {
            // If all successors of node are in the alreadyPlacedInOtherLayers set, choose this
            // node, since the list is sorted by rank
            if (remainingOutGoing[node.id] <= 0) {
                return node;
            }

        }
        // If no candidate was found return null
        return null;
    }

    /**
     * Sorts a list of nodes by its rank. The rank is defined as max(d⁺(v),max(d⁺(u):(u,v)∈ E)),
     * where d⁺(v) is the number of outgoing edges of a node v. Since the id-field is reused later
     * on, this function should be executed first
     */
    private void computeSortedNodes() {
        List<LNode> unsortedNodes = currentGraph.getLayerlessNodes();
        sortedLayerlessNodes = new ArrayList<>(unsortedNodes);
        // the id-field is reused later on, be careful
        for (LNode node : unsortedNodes) {
            node.id = getRank(node);
        }
        // sorting the list by the node-id /rank
        Collections.sort(sortedLayerlessNodes, new Comparator<LNode>() {
            public int compare(final LNode o1, final LNode o2) {
                // descending sort
                if (o1.id < o2.id) {
                    return 1;
                } else if (o1.id > o2.id) {
                    return -1;
                }
                return 0;
            }
        });

    }

    /**
     * Computes the rank of a node. The rank is defined as max(d⁺(v),max(d⁺(u):(u,v)∈ E)), where
     * d⁺(v) is the number of outgoing edges of a node v.
     *
     * @param node
     *            to compute the rank for
     * @return rank of the node
     */
    private Integer getRank(final LNode node) {
        int max = ((int) java.util.stream.StreamSupport.stream(node.getOutgoingEdges().spliterator(), false).count());
        int temp;
        LNode pre;

        // compute max of predecessors out-degree and out-degree of the current node
        for (LEdge preEdge : node.getIncomingEdges()) {
            pre = preEdge.getSource().getNode();
            temp = ((int) java.util.stream.StreamSupport.stream(pre.getOutgoingEdges().spliterator(), false).count());
            max = Math.max(max, temp);
        }

        return max;
    }

    /**
     * Computes the successors of every node and saves them in a list of sets of nodes. This list is
     * used in selectNode. Needs to be executed after {@link #computeSortedNodes()}.
     */
    private void computeSuccessors() {
        int i = 0;
        successors = new ArrayList<>();
        Set<LNode> currSucc = new HashSet<>();

        for (LNode node : sortedLayerlessNodes) {
            node.id = i;
            for (LEdge edge : node.getOutgoingEdges()) {
                currSucc.add(edge.getTarget().getNode());
            }
            // remove current node from accSucc to deal with self-loops
            currSucc.remove(node);
            successors.add(new HashSet<>(currSucc));
            currSucc.clear();
            i++;
        }
    }

    /**
     * Computes the in- and out-degree of every node. the values are used to determine widthUp and
     * widthCurrent. Needs to be executed after {@link #computeSuccessors()}.
     */
    private void computeDegrees() {
        inDegree = new int[sortedLayerlessNodes.size()];
        outDegree = new int[sortedLayerlessNodes.size()];

        for (LNode node : sortedLayerlessNodes) {
            inDegree[node.id] = ((int) java.util.stream.StreamSupport.stream(node.getIncomingEdges().spliterator(), false).count());
            outDegree[node.id] = ((int) java.util.stream.StreamSupport.stream(node.getOutgoingEdges().spliterator(), false).count());
        }
    }

    /**
     * Computes the minimum and the maximum node size. Needs to be computed after
     * {@link #computeSuccessors()}.
     */
    private void minMaxNodeSize() {
        // since in ELK all things are layered, left to right
        // a preprocessor also transposes the width and height and
        // we don't need to consider the layering direction and can take the
        // the y-size
        for (LNode node : sortedLayerlessNodes) {
            if (node.getType() != NodeType.NORMAL) {
                continue;
            }
            double size = node.getSize().y;
            minimumNodeSize = Math.min(minimumNodeSize, size);
            maximumNodeSize = Math.max(maximumNodeSize, size);
        }
    }

    /**
     * Computes the average normalized node size.
     *
     * @return average normalized node size
     */
    @SuppressWarnings("unused")
    private double avereageNodeSize() {
        double sum = 0;
        for (LNode node : sortedLayerlessNodes) {
            sum = sum + normSize[node.id];
        }
        return sum / sortedLayerlessNodes.size();
    }

    /**
     * Computes the size of every node according to the layering direction and normalizes it by the
     * minimum node size. Needs to be computed after {@link #computeSuccessors()} and
     * {@link #minimumNodeSize()}.
     */
    private void computeNormalizedSize() {
        normSize = new double[sortedLayerlessNodes.size()];
        for (LNode node : sortedLayerlessNodes) {
            normSize[node.id] = node.getSize().y / minimumNodeSize;
        }
    }

    /**
     * Computes the average out-degree of the graph. Should be computed before changing the
     * layerlessNodes list.
     *
     * @return average out-degree of the Graph
     */
    private float getAverageOutDegree() {
        float allOut = 0;
        for (LNode node : currentGraph.getLayerlessNodes()) {
            allOut += ((int) java.util.stream.StreamSupport.stream(node.getOutgoingEdges().spliterator(), false).count());
        }
        return allOut / currentGraph.getLayerlessNodes().size();
    }

    /**
     * Updates the information of the nodes, telling which has successors that are not placed. Is used when one layer is
     * finished, to eliminate edges in the same layer.
     *
     * @param currentLayer
     *            which is about to be finished
     */
    private void updateOutGoing(final Layer currentLayer) {
        for (LNode node : currentLayer.getNodes()) {
            for (LEdge edge : node.getIncomingEdges()) {
                int pos = edge.getSource().getNode().id;
                remainingOutGoing[pos] = remainingOutGoing[pos] - 1;
            }
        }
    }
}
