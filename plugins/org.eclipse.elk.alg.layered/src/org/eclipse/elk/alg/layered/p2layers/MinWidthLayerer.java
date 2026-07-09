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

import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Set;

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
import org.eclipse.elk.core.util.Pair;

import com.google.common.base.Predicate;
import com.google.common.collect.Lists;
import com.google.common.collect.Range;
import com.google.common.collect.Sets;

/**
 * Implementation of the heuristic MinWidth for solving the NP-hard minimum-width layering problem
 * with consideration of dummy nodes. MinWidth is based on the longest-path algorithm, which finds
 * layerings with the minimum height, but doesn't consider the width of the graph. MinWidth also
 * considers an upper bound on the width of a given graph. The upper bound isn't a "bound" in a
 * strict sense, as some layers might exceed its limit, if certain conditions are met.
 * 
 * Details are described in
 * <ul>
 * <li>Nikola S. Nikolov, Alexandre Tarassov, and Jürgen Branke. 2005. In search for efficient
 * heuristics for minimum-width graph layering with consideration of dummy nodes. J. Exp.
 * Algorithmics 10, Article 2.7 (December 2005). DOI=10.1145/1064546.1180618
 * http://doi.acm.org/10.1145/1064546.1180618.</li>
 * </ul>
 * 
 * MinWidth takes two additional parameters, which can be configured as a property:
 * <ul>
 * <li>Upper Bound On Width {@link LayeredOptions#UPPER_BOUND_ON_WIDTH} – Defines a loose upper bound on
 * the width of the MinWidth layerer. Defaults to -1 (special value for using both 1, 2, 3 and 4 as
 * values and choosing the narrowest layering afterwards), lower bound is 1.</li>
 * <li>Upper Layer Estimation Scaling Factor
 * {@link LayeredOptions#UPPER_LAYER_ESTIMATION_SCALING_FACTOR} – Multiplied with
 * {@link LayeredOptions#UPPER_BOUND_ON_WIDTH} for defining an upper bound on the width of layers which
 * haven't been placed yet, but whose maximum width had been (roughly) estimated by the MinWidth
 * algorithm. Compensates for too high estimations. Defaults to -1 (special value for using both 1
 * and 2 as values and choosing the narrowest layering afterwards), lower bound is 1.</li>
 * </ul>
 * 
 * This version of the algorithm, however, differs from the one described in the paper as it
 * considers the actual size of the nodes of the graph in order to handle real world use cases of
 * graphs a little bit better. The approach is based on Marc Adolf's version in his implementation
 * of the heuristic {@link StretchWidthLayerer}. Some changes include:
 * <ul>
 * <li>estimating the sizes of dummy nodes by taking the edge spacing of the {@link LGraph} into
 * account,</li>
 * <li>finding the narrowest real node of the graph and normalizing all the widths of the nodes of
 * the graph (real and dummy) in relation to this node,</li>
 * <li>computing the average size of all real nodes (we don't know the number of dummy nodes in
 * advance),</li>
 * <li>using this average as a factor for the ubw-value given by the user in order to adjust the
 * boundary to our new approach (using the result of this multiplication instead of the given value
 * of ubw thus changes the condition to start a new layer from the paper slightly).</li>
 * </ul>
 * 
 * <dl>
 * <dt>Precondition:</dt>
 * <dd>the graph has no cycles, but might contain self-loops</dd>
 * <dt>Postcondition:</dt>
 * <dd>all nodes have been assigned a layer such that edges connect only nodes from layers with
 * increasing indices</dd>
 * </dl>
 */
public final class MinWidthLayerer implements ILayoutPhase<LayeredPhases, LGraph> {

    /**
     * Recommended values for the algorithm suggested bei Nikolov et al. after a parameter study,
     * see:
     * 
     * <ul>
     * <li>Alexandre Tarassov, Nikola S. Nikolov, and Jürgen Branke. 2004. A Heuristic for
     * Minimum-Width Graph Layering with Consideration of Dummy Nodes. Experimental and Efficient
     * Algorithms, Third International Workshop, WEA 2004, Lecture Notes in Computer Science 3059.
     * Springer-Verlag, New York, 570-583. DOI=10.1007/978-3-540-24838-5_42
     * http://dx.doi.org/10.1007/978-3-540-24838-5_42.</li>
     * </ul>
     */
    private static final Range<Integer> UPPERBOUND_ON_WIDTH_RANGE = Range.closed(1, 4);
    private static final Range<Integer> COMPENSATOR_RANGE = Range.closed(1, 2);

    // Some variables used for considering real node sizes, see below.
    private double dummySize;
    private double minimumNodeSize;
    private double[] normSize;
    private double avgSize;

    /**
     * Checks whether an edge is a self-loop (i.e. source node == target node). Used internally, as
     * ELK Layered doesn't remove self-loops in its cycle breaking phase, but MinWidth expects only
     * real directed acyclic graphs (DAGs).
     */
    private SelfLoopPredicate isSelfLoopTest = new SelfLoopPredicate();

    // Degrees of nodes without self loops, indexed by node.id
    private int[] inDegree;
    private int[] outDegree;

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
        progressMonitor.begin("MinWidth layering", 1);

        List<Layer> layers = layeredGraph.getLayers();
        List<LNode> notInserted = layeredGraph.getLayerlessNodes();

        // The algorithm requires DAG G = (V, E). In this version self-loops are allowed (as we're
        // going to filter them). Additional properties as described above (called UBW and c in the
        // original paper):
        final int upperBoundOnWidth = layeredGraph.getProperty(LayeredOptions.LAYERING_MIN_WIDTH_UPPER_BOUND_ON_WIDTH);
        final int compensator =
                layeredGraph.getProperty(LayeredOptions.LAYERING_MIN_WIDTH_UPPER_LAYER_ESTIMATION_SCALING_FACTOR);

        // First step to consider the real size of nodes: Initialize the dummy size with the spacing
        // properties
        dummySize = layeredGraph.getProperty(LayeredOptions.SPACING_EDGE_EDGE).doubleValue();

        // Compute the minimum nodes size (of the real nodes). We're going to use this value in the
        // next step to normalize the different node sizes.
        minimumNodeSize = Double.POSITIVE_INFINITY;
        for (LNode node : notInserted) {
            if (node.getType() != NodeType.NORMAL) {
                continue;
            }
            double size = node.getSize().y;
            minimumNodeSize = Math.min(minimumNodeSize, size);
        }
        // The minimum nodes size might be zero. If This is the case, then simply don't normalize the
        // node sizes.
        minimumNodeSize = Math.max(1, minimumNodeSize);
        
        // We initialize the nodes' id and use it to refer to its in- and out-degree stored each in
        // an array. We also compute the size of each node in relation to the smallest real node in
        // the graph (normalized size) and store it in the same way.
        int numOfNodes = notInserted.size();
        inDegree = new int[numOfNodes];
        outDegree = new int[numOfNodes];
        normSize = new double[numOfNodes];
        int i = 0;
        avgSize = 0;
        for (LNode node : notInserted) {
            // Warning: LNode.id is being redefined here!
            node.id = i++;
            inDegree[node.id] = countEdgesExceptSelfLoops(node.getIncomingEdges());
            outDegree[node.id] = countEdgesExceptSelfLoops(node.getOutgoingEdges());
            normSize[node.id] = node.getSize().y / minimumNodeSize;
            // The average size of a node will also be based on the normalized size.
            avgSize += normSize[node.id];
        }
        // normalize dummy size, too:
        dummySize /= minimumNodeSize;
        // Divide sum of normalized node sizes by the number of nodes to get an actual mean.
        avgSize /= numOfNodes;

        // Precalculate the successors of all nodes (as a Set) and put them in a list.
        List<Set<LNode>> nodeSuccessors = precalcSuccessors(notInserted);

        // Guarantee ConditionSelect from the paper, which states that nodes with maximum out-degree
        // should be preferred during layer placement, by ordering the nodes by descending maximum
        // out-degree in advance.
        notInserted.sort(Collections.reverseOrder(new MinOutgoingEdgesComparator()));

        // minimum width of a layer of maximum size in a computed layering (primary criterion used
        // for comparison, if more than one layering is computed). It's a double as it takes in
        // account the actual width based on the normalized size of the nodes.
        double minWidth = Double.POSITIVE_INFINITY;
        // minimum number of layers in a computed layering {@code minWidth} (secondary
        // criterion used for comparison, if more than one layering is computed).
        int minNumOfLayers = Integer.MAX_VALUE;
        // holding the currently chosen candidate for the final layering as a List
        List<List<LNode>> candidateLayering = null;

        // At first blindly set the parameters for the loose upper bound and the compensator to the
        // exact values, which have been configured via their respective properties, so that only
        // one layering will be computed …
        int ubwStart = upperBoundOnWidth;
        int ubwEnd = upperBoundOnWidth;
        int cStart = compensator;
        int cEnd = compensator;

        // ... then check, whether any special values (i.e. negative values, which aren't valid) have
        // been used for the properties. In that case use the recommended ranges described above …
        if (upperBoundOnWidth < 0) {
            ubwStart = UPPERBOUND_ON_WIDTH_RANGE.lowerEndpoint();
            ubwEnd = UPPERBOUND_ON_WIDTH_RANGE.upperEndpoint();
        }
        if (compensator < 0) {
            cStart = COMPENSATOR_RANGE.lowerEndpoint();
            cEnd = COMPENSATOR_RANGE.upperEndpoint();
        }

        // … Depending on the start- and end-values, this nested for-loop will last for up to 8
        // iterations resulting in one, two, four or eight different layerings.
        for (int ubw = ubwStart; ubw <= ubwEnd; ubw++) {
            for (int c = cStart; c <= cEnd; c++) {
                Pair<Double, List<List<LNode>>> result =
                        computeMinWidthLayering(ubw, c, notInserted, nodeSuccessors);
                double newWidth = result.getFirst();
                List<List<LNode>> layering = result.getSecond();

                // Important if more than one layering is computed: replace the current candidate
                // layering with a newly computed one, if it is narrower or has the same maximum
                // width but less layers.
                int newNumOfLayers = layering.size();
                if (newWidth < minWidth
                        || (newWidth == minWidth && newNumOfLayers < minNumOfLayers)) {
                    minWidth = newWidth;
                    minNumOfLayers = newNumOfLayers;
                    candidateLayering = layering;
                }
            }
        }

        // Finally, add the winning layering to the Klay layered data structures.
        for (List<LNode> layerList : candidateLayering) {
            Layer currentLayer = new Layer(layeredGraph);
            for (LNode node : layerList) {
                node.setLayer(currentLayer);
            }
            layers.add(currentLayer);
        }

        // The algorithm constructs the layering bottom up, but ElkLayered expects the list of
        // layers to be ordered top down.
        Collections.reverse(layers);
        // After the algorithm, there should be no nodes left to be put in a layer, so we're gonna
        // delete them.
        notInserted.clear();

        progressMonitor.done();
    }

    /**
     * Calculates for a given Collection of {@link LNode} all its successors (i.e. a Set of
     * {@link LNode}s) without self-loops.
     * 
     * @param nodes
     *            a Collection of {@link LNode}
     * @return List of Set of successor {@link LNode}s in order of the given nodes
     */
    private List<Set<LNode>> precalcSuccessors(final Collection<LNode> nodes) {
        List<Set<LNode>> successors = Lists.newArrayListWithCapacity(nodes.size());

        for (LNode node : nodes) {

            Set<LNode> outNodes = Sets.newHashSet();
            Iterable<LEdge> outEdges = node.getOutgoingEdges();

            for (LEdge edge : outEdges) {
                if (!isSelfLoopTest.apply(edge)) {
                    outNodes.add(edge.getTarget().getNode());
                }
            }

            successors.add(outNodes);
        }

        return successors;
    }

    /**
     * 
     * Computes a layering (as a List of Lists) for a given Iterable of {@link LNode} according to
     * the MinWidth-heuristic and considering actual node sizes.
     * 
     * @param upperBoundOnWidth
     *            Defines a loose upper bound on the width of the MinWidth layerer. Uses integer
     *            values as in the original approach described in the paper, as this bound will
     *            automatically be multiplied internally with the average normalized node size as
     *            part of the new approach considering the actual sizes of nodes.
     * @param compensator
     *            Multiplied with {@code upperBoundOnWidth} for defining an upper bound on the width
     *            of layers which haven't been determined yet, but whose maximum width had been
     *            (roughly) estimated by the MinWidth algorithm. Compensates for too high
     *            estimations.
     * @param nodes
     *            Iterable of all nodes of the Graph. The {@code id} of the nodes have to be set to
     *            the index where the respective Set of successor-nodes are stored in the List
     *            {@code nodeSuccessors}.
     * @param nodeSuccessors
     *            precomputed List of Set of successor-nodes of the elements in the Iterable
     *            {@code nodes}.
     * @return a pair of a double representing the maximum width of the resulting layering
     *         (normalized by the smallest real node) and the layering itself as a list of list of
     *         nodes
     */
    private Pair<Double, List<List<LNode>>> computeMinWidthLayering(final int upperBoundOnWidth,
            final int compensator, final Iterable<LNode> nodes,
            final List<Set<LNode>> nodeSuccessors) {

        List<List<LNode>> layers = Lists.newArrayList();
        Set<LNode> unplacedNodes = Sets.newLinkedHashSet(nodes);

        // One of the deviations from the paper is, that our upper bound is taking node sizes into
        // account:
        double ubwConsiderSize = upperBoundOnWidth * avgSize;

        // in- and out-degree of the currently considered node, see while-loop below
        int inDeg = 0;
        int outDeg = 0;

        // The actual algorithm from the paper begins here:
        // In the Paper the first Set contains all nodes, which have already been placed (in this
        // version we consider only the nodes already placed in the current layer), and the
        // second contains all nodes already placed in layers which have been determined before the
        // currentLayer.
        Set<LNode> alreadyPlacedInCurrentLayer = Sets.newHashSet();
        Set<LNode> alreadyPlacedInOtherLayers = Sets.newHashSet();

        // Set up the first layer (algorithm is bottom up, so the List layer is going to be reversed
        // at the end.
        List<LNode> currentLayer = Lists.newArrayList();

        // Initial values for the width of the current layer and the estimated width of the coming
        // layers
        double widthCurrent = 0;
        double widthUp = 0;

        // Parameters needed for computing the width of a layering including dummy nodes:
        double maxWidth = 0;
        double realWidth = 0;
        // Number of "started" edges that did not "finish" yet (now multiplied with the normalized
        // dummy size to consider actual node sizes)
        double currentSpanningEdges = 0;
        double goingOutFromThisLayer = 0;
        // No need for a variable "comingIntoThisLayer" as "widthUp" already gets the job done.

        while (!unplacedNodes.isEmpty()) {
            // Find a node, whose edges only point to nodes in the Set alreadyPlacedInOtherLayers;
            // will return {@code null} if such a node doesn't exist.
            LNode currentNode =
                    selectNode(unplacedNodes, nodeSuccessors, alreadyPlacedInOtherLayers);

            // If a node is found in the previous step:
            if (currentNode != null) {
                unplacedNodes.remove(currentNode);
                currentLayer.add(currentNode);
                alreadyPlacedInCurrentLayer.add(currentNode);

                outDeg = this.outDegree[currentNode.id];
                // Take node sizes in account: use the normalized size of current node and the
                // normalized dummy size for each edge
                widthCurrent += normSize[currentNode.id] - outDeg * dummySize;

                inDeg = this.inDegree[currentNode.id];
                // Take node sizes in account: use the normalized normalized dummy size for each
                // edge
                widthUp += inDeg * dummySize;

                goingOutFromThisLayer += outDeg * dummySize;

                realWidth += normSize[currentNode.id];
            }

            // Go to the next layer if,
            // 1) no current node has been selected,
            // 2) there are no unplaced nodes left (last iteration of the while-loop),
            // 3) The conditionGoUp from the paper (with the difference of ubw being multiplied with
            // the
            // average normalized node size) is satisfied, i.e.
            // 3.1) the width of the current layer is greater than the upper bound on the width and
            // the number of dummy nodes in the layer can't be reduced, as only nodes with no
            // outgoing edges are left for being considered for the current layer; or:
            // 3.2) The estimated width of the not yet determined layers is greater than the
            // scaling factor/compensator times the upper bound on the width.
            if (currentNode == null || unplacedNodes.isEmpty() 
                    || (widthCurrent >= ubwConsiderSize && normSize[currentNode.id] > outDeg * dummySize)
                    || widthUp >= compensator * ubwConsiderSize) {
                layers.add(currentLayer);
                currentLayer = Lists.newArrayList();
                alreadyPlacedInOtherLayers.addAll(alreadyPlacedInCurrentLayer);
                alreadyPlacedInCurrentLayer.clear();

                // Remove all edges from the dummy node count, which are starting at a node placed
                // in this layer …
                currentSpanningEdges -= goingOutFromThisLayer;
                // … Now we have the actual dummy node count (or rather the sum of their widths) for
                // this layer and can add it to the real nodes for comparing the width.
                maxWidth = Math.max(maxWidth, currentSpanningEdges * dummySize + realWidth);
                // In the next iteration we have to consider new dummy nodes from edges coming into
                // the layer we've just finished.
                currentSpanningEdges += widthUp;

                widthCurrent = widthUp;
                widthUp = 0;
                goingOutFromThisLayer = 0;
                realWidth = 0;
            }
        }

        return Pair.of(maxWidth, layers);
    }

    /**
     * Returns the first {@link LNode} in the given Set, whose outgoing edges end only in nodes of
     * the Set {@code targets}. Self-loops are ignored.
     * 
     * Warning: Returns {@code null}, if such a node doesn't exist.
     * 
     * @param nodes
     *            Set to choose {@link LNode} from
     * @param targets
     *            Set of {@link LNode}
     * @return chosen {@link LNode} from {@code nodes}, whose outgoing edges all end in a node
     *         contained in {@code targets}. Returns {@code null}, if such a node doesn't exist.
     */
    private LNode selectNode(final Set<LNode> nodes, final List<Set<LNode>> successors,
            final Set<LNode> targets) {

        for (LNode node : nodes) {
            if (targets.containsAll(successors.get(node.id))) {
                return node;
            }
        }
        return null;
    }

    /**
     * Returns the number of {@link LEdge} edges in the given Iterable, but ignores self-loops.
     * 
     * @param edges
     *            Iterable whose edges without self-loops are to be counted
     * @return number of {@link LEdge} edges without self-loops
     */
    private int countEdgesExceptSelfLoops(final Iterable<LEdge> edges) {
        int i = 0;
        for (LEdge edge : edges) {
            if (!isSelfLoopTest.apply(edge)) {
                i++;
            }
        }
        return i;
    }

    /**
     * Checks whether an edge is a self-loop (i.e. source node == target node).
     */
    private class SelfLoopPredicate implements Predicate<LEdge> {

        @Override
        public boolean apply(final LEdge input) {
            return input.getSource().getNode().equals(input.getTarget().getNode());
        }

    }

    /**
     * Comparator for determining whether a {@link LNode} has less outgoing edges than another one.
     * Requires the LNode property {@link LNode#id} to be set to the number of outgoing edges of the
     * node.
     * 
     * @author mic
     */
    private class MinOutgoingEdgesComparator implements Comparator<LNode> {
        @Override
        public int compare(final LNode o1, final LNode o2) {
            int outs1 = outDegree[o1.id];
            int outs2 = outDegree[o2.id];

            if (outs1 < outs2) {
                return -1;
            }
            if (outs1 == outs2) {
                return 0;
            }
            return 1;
        }
    }

}
