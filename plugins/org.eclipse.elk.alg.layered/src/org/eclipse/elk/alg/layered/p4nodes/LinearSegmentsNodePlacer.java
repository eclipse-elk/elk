/*******************************************************************************
 * Copyright (c) 2010, 2015 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.layered.p4nodes;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.ListIterator;

import org.eclipse.elk.alg.layered.DebugUtil;
import org.eclipse.elk.alg.layered.LayeredPhases;
import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LNode.NodeType;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.intermediate.IntermediateProcessorStrategy;
import org.eclipse.elk.alg.layered.options.GraphProperties;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.layered.options.Spacings;
import org.eclipse.elk.core.alg.ILayoutPhase;
import org.eclipse.elk.core.alg.LayoutProcessorConfiguration;
import org.eclipse.elk.core.util.IElkProgressMonitor;
import org.eclipse.elk.graph.properties.Property;

import com.google.common.collect.Lists;

/**
 * Node placement implementation that aligns long edges using linear segments. Inspired by Section 4 of
 * <ul>
 *   <li>Georg Sander, A fast heuristic for hierarchical Manhattan layout. In <i>Proceedings of the
 *     Symposium on Graph Drawing (GD'95)</i>, LNCS vol. 1027, pp. 447-458, Springer, 1996.</li>
 * </ul>
 *
 * <dl>
 * <dt>Precondition:</dt>
 * <dd>the graph has a proper layering with optimized nodes ordering; ports are properly arranged</dd>
 * <dt>Postcondition:</dt>
 * <dd>each node is assigned a vertical coordinate such that no two nodes overlap; the size of each
 * layer is set according to the area occupied by contained nodes; the height of the graph is set to
 * the maximal layer height</dd>
 * </dl>
 *
 * @author msp
 * @author grh
 * @author cds
 * @author ima
 */
public final class LinearSegmentsNodePlacer implements ILayoutPhase<LayeredPhases, LGraph> {

    /**
     * A linear segment contains a single regular node or all dummy nodes of a long edge.
     */
    public static class LinearSegment implements Comparable<LinearSegment> {
        /** Nodes of the linear segment. */
        private List<LNode> nodes = Lists.newArrayList();
        /** Identifier value, used as index in the segments array. */
        private int id;
        /** Index in the previous layer. Used for cycle avoidance. */
        private int indexInLastLayer = -1;
        /** The last layer where a node belonging to this segment was discovered. Used for cycle
         *  avoidance. */
        private int lastLayer = -1;
        /** The accumulated force of the contained nodes. */
        private double deflection;
        /** The current weight of the contained nodes. */
        private int weight;
        /** The reference segment, if this has been unified with another. */
        private LinearSegment refSegment;
        /** The nodetype contained in this linear segment. */
        private NodeType nodeType;

        /**
         * Determine the reference segment for the region to which this segment is associated.
         *
         * @return the region segment
         */
        LinearSegment region() {
            LinearSegment seg = this;
            while (seg.refSegment != null) {
                seg = seg.refSegment;
            }
            return seg;
        }

        /**
         * Splits this linear segment before the given node. The returned segment contains all nodes
         * from the given node onward, with their ID set to the new segment's ID. Those nodes are
         * removed from this segment.
         *
         * @param node
         *            the node to split the segment at.
         * @param newId
         *            the new segment's id.
         * @return new linear segment with ID {@code -1} and all nodes from {@code node} onward.
         */
        LinearSegment split(final LNode node, final int newId) {
            int nodeIndex = nodes.indexOf(node);

            // Prepare the new segment
            LinearSegment newSegment = new LinearSegment();
            newSegment.id = newId;

            // Move nodes to the new segment
            ListIterator<LNode> iterator = nodes.listIterator(nodeIndex);
            while (iterator.hasNext()) {
                LNode movedNode = iterator.next();
                movedNode.id = newId;
                newSegment.nodes.add(movedNode);
                iterator.remove();
            }

            return newSegment;
        }

        @Override
        public String toString() {
            return "ls" + nodes.toString();
        }
        
        @Override
        public int compareTo(final LinearSegment other) {
            return this.id - other.id;
        }

        @Override
        public boolean equals(final Object object) {
            if (object instanceof LinearSegment) {
                LinearSegment other = (LinearSegment) object;
                return this.id == other.id;
            }
            return false;
        }

        @Override
        public int hashCode() {
            return id;
        }

    }

    /** additional processor dependencies for graphs with hierarchical ports. */
    private static final LayoutProcessorConfiguration<LayeredPhases, LGraph> HIERARCHY_PROCESSING_ADDITIONS =
        LayoutProcessorConfiguration.<LayeredPhases, LGraph>create()
            .addBefore(LayeredPhases.P5_EDGE_ROUTING,
                    IntermediateProcessorStrategy.HIERARCHICAL_PORT_POSITION_PROCESSOR);

    @Override
    public LayoutProcessorConfiguration<LayeredPhases, LGraph> getLayoutProcessorConfiguration(final LGraph graph) {
        if (graph.getProperty(InternalProperties.GRAPH_PROPERTIES).contains(GraphProperties.EXTERNAL_PORTS)) {
            return HIERARCHY_PROCESSING_ADDITIONS;
        } else {
            return null;
        }
    }

    /** property for maximal priority of incoming edges. */
    private static final Property<Integer> INPUT_PRIO = new Property<Integer>(
            "linearSegments.inputPrio", 0);
    /** property for maximal priority of outgoing edges. */
    private static final Property<Integer> OUTPUT_PRIO = new Property<Integer>(
            "linearSegments.outputPrio", 0);

    /** array of sorted linear segments. */
    private LinearSegment[] linearSegments;

    /** Spacing values. */
    private Spacings spacings;

    @Override
    public void process(final LGraph layeredGraph, final IElkProgressMonitor monitor) {
        monitor.begin("Linear segments node placement", 1);

        spacings = layeredGraph.getProperty(InternalProperties.SPACINGS);

        // sort the linear segments of the layered graph
        sortLinearSegments(layeredGraph, monitor);

        // create an unbalanced placement from the sorted segments
        createUnbalancedPlacement(layeredGraph);

        // balance the placement
        balancePlacement(layeredGraph);

        // post-process the placement for small corrections
        postProcess(layeredGraph);

        // release the created resources
        linearSegments = null;
        spacings = null;
        monitor.done();
    }

    // /////////////////////////////////////////////////////////////////////////////
    // Linear Segments Creation

    /**
     * Sorts the linear segments of the given layered graph by finding a topological ordering in the
     * corresponding segment ordering graph.
     *
     * @param layeredGraph
     *            layered graph to process
     * @param monitor
     *            our progress monitor.
     * @return a sorted array of linear segments
     */
    private LinearSegment[] sortLinearSegments(final LGraph layeredGraph, final IElkProgressMonitor monitor) {
        // set the identifier and input / output priority for all nodes
        List<LinearSegment> segmentList = Lists.newArrayList();
        for (Layer layer : layeredGraph) {
            for (LNode node : layer) {
                node.id = -1;
                int inprio = Integer.MIN_VALUE, outprio = Integer.MIN_VALUE;
                for (LPort port : node.getPorts()) {
                    for (LEdge edge : port.getIncomingEdges()) {
                        int prio = edge.getProperty(LayeredOptions.PRIORITY_STRAIGHTNESS);
                        inprio = Math.max(inprio, prio);
                    }
                    for (LEdge edge : port.getOutgoingEdges()) {
                        int prio = edge.getProperty(LayeredOptions.PRIORITY_STRAIGHTNESS);
                        outprio = Math.max(outprio, prio);
                    }
                }
                node.setProperty(INPUT_PRIO, inprio);
                node.setProperty(OUTPUT_PRIO, outprio);
            }
        }

        // create linear segments for the layered graph, ignoring odd port side dummies
        int nextLinearSegmentID = 0;
        for (Layer layer : layeredGraph) {
            for (LNode node : layer) {
                // Test for the node ID; calls to fillSegment(...) may have caused the node ID
                // to be != -1.
                if (node.id < 0) {
                    LinearSegment segment = new LinearSegment();
                    segment.id = nextLinearSegmentID++;
                    fillSegment(node, segment);
                    segmentList.add(segment);
                }
            }
        }

        // create and initialize segment ordering graph
        List<List<LinearSegment>> outgoingList = Lists.newArrayListWithCapacity(segmentList.size());
        List<Integer> incomingCountList = Lists.newArrayListWithCapacity(segmentList.size());
        for (int i = 0; i < segmentList.size(); i++) {
            outgoingList.add(new ArrayList<LinearSegment>());
            incomingCountList.add(0);
        }

        // create edges for the segment ordering graph
        createDependencyGraphEdges(monitor, layeredGraph, segmentList, outgoingList, incomingCountList);

        // turn lists into arrays
        LinearSegment[] segments = segmentList.toArray(new LinearSegment[segmentList.size()]);

        @SuppressWarnings("unchecked")
        List<LinearSegment>[] outgoing = outgoingList.toArray(new List[outgoingList.size()]);

        int[] incomingCount = new int[incomingCountList.size()];
        for (int i = 0; i < incomingCount.length; i++) {
            incomingCount[i] = incomingCountList.get(i);
        }

        // gather the sources of the segment ordering graph
        int nextRank = 0;
        List<LinearSegment> noIncoming = Lists.newArrayList();
        for (int i = 0; i < segments.length; i++) {
            if (incomingCount[i] == 0) {
                noIncoming.add(segments[i]);
            }
        }

        // find a topological ordering of the segment ordering graph
        int[] newRanks = new int[segments.length];
        while (!noIncoming.isEmpty()) {
            LinearSegment segment = noIncoming.remove(0);
            newRanks[segment.id] = nextRank++;

            while (!outgoing[segment.id].isEmpty()) {
                LinearSegment target = outgoing[segment.id].remove(0);
                incomingCount[target.id]--;

                if (incomingCount[target.id] == 0) {
                    noIncoming.add(target);
                }
            }
        }

        // apply the new ordering to the array of linear segments
        linearSegments = new LinearSegment[segments.length];
        for (int i = 0; i < segments.length; i++) {
            assert outgoing[i].isEmpty();
            LinearSegment ls = segments[i];
            int rank = newRanks[i];
            linearSegments[rank] = ls;
            ls.id = rank;
            for (LNode node : ls.nodes) {
                node.id = rank;
            }
        }

        return linearSegments;
    }

    /**
     * Fills the dependency graph with dependencies. If a dependency would introduce a cycle, the
     * offending linear segment is split into two linear segments.
     *
     * @param monitor
     *            our progress monitor.
     * @param layeredGraph
     *            the layered graph.
     * @param segmentList
     *            the list of segments. Updated to include the newly created linear segments.
     * @param outgoingList
     *            the lists of outgoing dependencies for each segment. This essentially encodes the
     *            edges of the dependency graph.
     * @param incomingCountList
     *            the number of incoming dependencies for each segment.
     */
    private void createDependencyGraphEdges(final IElkProgressMonitor monitor, final LGraph layeredGraph,
            final List<LinearSegment> segmentList, final List<List<LinearSegment>> outgoingList,
            final List<Integer> incomingCountList) {

        /*
         * There's some <scaryVoice> faaaancy </scaryVoice> stuff going on here. Basically, we go
         * through all the layers, from left to right. In each layer, we go through all the nodes.
         * For each node, we retrieve the linear segment it's part of and add a dependency to the
         * next node's linear segment. So far so good.
         *
         * This works perfectly fine as long as we assume that the relative order of linear segments
         * doesn't change from one layer to the next. However, since the introduction of north /
         * south port dummies, it can. We now have to avoid creating cycles in the dependency graph.
         * This is done by remembering the indices of each linear segment in the previous layer.
         * When we encounter a segment x, we check if there is a segment y that came before x in the
         * previous layer. (that would introduce a cycle) If that's the case, we split x at the
         * current layer, resulting in two segments, x1 and x2, x2 starting at the current layer.
         * Now, we proceed as usual, adding a dependency from x2 to y. But we have avoided a cycle
         * because y does not depend on x2, but on x1.
         */

        int nextLinearSegmentID = segmentList.size();
        int layerIndex = 0;
        for (Layer layer : layeredGraph) {
            List<LNode> nodes = layer.getNodes();
            if (nodes.isEmpty()) {
                // Ignore empty layers
                continue;
            }

            Iterator<LNode> nodeIter = nodes.iterator();
            int indexInLayer = 0;

            // We carry the previous node with us for dependency management
            LNode previousNode = null;

            // Get the layer's first node
            LNode currentNode = nodeIter.next();
            LinearSegment currentSegment = null;

            while (currentNode != null) {
                // Get the current node's segment
                currentSegment = segmentList.get(currentNode.id);

                /*
                 * Check if we have a cycle. That's the case if the following holds: - The current
                 * segment appeared in the previous layer as well. - In the previous layer, we find
                 * a segment after the current segment that appears before the current segment in
                 * the current layer.
                 */
                if (currentSegment.indexInLastLayer >= 0) {
                    LinearSegment cycleSegment = null;
                    Iterator<LNode> cycleNodesIter = layer.getNodes()
                            .listIterator(indexInLayer + 1);
                    while (cycleNodesIter.hasNext()) {
                        LNode cycleNode = cycleNodesIter.next();
                        cycleSegment = segmentList.get(cycleNode.id);

                        if (cycleSegment.lastLayer == currentSegment.lastLayer
                                && cycleSegment.indexInLastLayer < currentSegment.indexInLastLayer) {

                            break;
                        } else {
                            cycleSegment = null;
                        }
                    }

                    // If we have found a cycle segment, we need to split the current linear segment
                    if (cycleSegment != null) {
                        // Update the current segment before it's split
                        if (previousNode != null) {
                            incomingCountList.set(currentNode.id,
                                    incomingCountList.get(currentNode.id) - 1);
                            outgoingList.get(previousNode.id).remove(currentSegment);
                        }

                        currentSegment = currentSegment.split(currentNode, nextLinearSegmentID++);
                        segmentList.add(currentSegment);
                        outgoingList.add(new ArrayList<LinearSegment>());

                        if (previousNode != null) {
                            outgoingList.get(previousNode.id).add(currentSegment);
                            incomingCountList.add(1);
                        } else {
                            incomingCountList.add(0);
                        }
                    }
                }

                // Now add a dependency to the next node, if any
                LNode nextNode = null;
                if (nodeIter.hasNext()) {
                    nextNode = nodeIter.next();
                    LinearSegment nextSegment = segmentList.get(nextNode.id);

                    outgoingList.get(currentNode.id).add(nextSegment);
                    incomingCountList.set(nextNode.id, incomingCountList.get(nextNode.id) + 1);
                }

                // Update segment's layer information
                currentSegment.lastLayer = layerIndex;
                currentSegment.indexInLastLayer = indexInLayer++;

                // Cycle nodes
                previousNode = currentNode;
                currentNode = nextNode;
            }

            layerIndex++;
        }

        // Write debug output graph
        // elkjs-exclude-start
        if (monitor.isLoggingEnabled()) {
            DebugUtil.logDebugGraph(monitor, layeredGraph, segmentList, outgoingList);
        }
        // elkjs-exclude-end
    }

    /**
     * Put a node into the given linear segment and check for following parts of a long edge.
     *
     * @param node
     *            the node to put into the linear segment
     * @param segment
     *            a linear segment
     * @return {@code true} if the given node was not already part of another segment and was thus
     *         added to the given segment.
     */
    private boolean fillSegment(final LNode node, final LinearSegment segment) {
        NodeType nodeType = node.getType();

        if (node.id >= 0) {
            // The node is already part of another linear segment
            return false;
        } else {
            // Add the node to the given linear segment
            node.id = segment.id;
            segment.nodes.add(node);
        }
        segment.nodeType = nodeType;

        if (nodeType == NodeType.LONG_EDGE
                || nodeType == NodeType.NORTH_SOUTH_PORT) {

            // This is a LONG_EDGE, NORTH_SOUTH_PORT or BIG_NODE dummy; check if any of its
            // successors are of one of these types too. If so, we can form a linear segment
            // with one of them. (not with more than one, though)
            // Note 1: LONG_EDGES and NORTH_SOUTH_PORTs can share a common linear segment
            // Note 2: we must take care not to make a segment out of nodes that are in the same layer
            // Note 3: for BIG_NODEs also the first BIG_NODE_INITIAL which is no actual dummy node has
            // to be considered here
            for (LPort sourcePort : node.getPorts()) {
                for (LPort targetPort : sourcePort.getSuccessorPorts()) {
                    LNode targetNode = targetPort.getNode();
                    NodeType targetNodeType = targetNode.getType();

                    if (node.getLayer() != targetNode.getLayer()) {
                        // current no bignode and next node is LONG_EDGE and NORTH_SOUTH_PORT
                        if (targetNodeType == NodeType.LONG_EDGE
                                || targetNodeType == NodeType.NORTH_SOUTH_PORT) {
                            if (fillSegment(targetNode, segment)) {
                                // We just added another node to this node's linear segment.
                                // That's quite enough.
                                return true;
                            }
                        }
                    }
                }
            }
        }

        return true;
    }

    // /////////////////////////////////////////////////////////////////////////////
    // Unbalanced Placement

    /**
     * Creates an unbalanced placement for the sorted linear segments.
     *
     * @param layeredGraph
     *            the layered graph to create an unbalanced placement for.
     */
    private void createUnbalancedPlacement(final LGraph layeredGraph) {
        // How many nodes are currently placed in each layer
        int[] nodeCount = new int[layeredGraph.getLayers().size()];

        // The type of the node most recently placed in a given layer
        NodeType[] recentNodeType = new NodeType[layeredGraph.getLayers().size()];
        LNode[] recentNode = new LNode[layeredGraph.getLayers().size()];

        // Iterate through the linear segments (in proper order!) and place them
        for (LinearSegment segment : linearSegments) {
            // Determine the uppermost placement for the linear segment
            double uppermostPlace = 0.0f;
            for (LNode node : segment.nodes) {
                int layerIndex = node.getLayer().getIndex();
                nodeCount[layerIndex]++;

                // Calculate how much space to leave between the linear segment and the last
                // node of the given layer
                double spacing = layeredGraph.getProperty(LayeredOptions.SPACING_EDGE_EDGE);
                if (nodeCount[layerIndex] > 0) {
                    if (recentNode[layerIndex] != null) {
                        //spacing = spacings.getVerticalSpacing(recentNodeType[layerIndex], nodeType);
                        spacing = spacings.getVerticalSpacing(recentNode[layerIndex], node);
                    }
                }

                uppermostPlace = Math.max(uppermostPlace, node.getLayer().getSize().y + spacing);
            }

            // Apply the uppermost placement to all elements
            for (LNode node : segment.nodes) {
                // Set the node position
                node.getPosition().y = uppermostPlace + node.getMargin().top;

                // Adjust layer size
                Layer layer = node.getLayer();
                layer.getSize().y = uppermostPlace + node.getMargin().top
                        + node.getSize().y + node.getMargin().bottom;

                recentNodeType[layer.getIndex()] = node.getType();
                recentNode[layer.getIndex()] = node;
            }
        }
    }


    // /////////////////////////////////////////////////////////////////////////////
    // Balanced Placement

    /** Definition of balancing modes. */
    private static enum Mode {
        FORW_PENDULUM,
        BACKW_PENDULUM,
        RUBBER;
    }

    /** factor for threshold after which balancing is aborted. */
    private static final double THRESHOLD_FACTOR = 20.0;
    /** the minimal number of iterations in pendulum mode. */
    private static final int PENDULUM_ITERS = 4;
    /** the number of additional iterations after the abort condition was met. */
    private static final int FINAL_ITERS = 3;

    /**
     * Balance the initial placement by force-based movement of regions. First perform <em>pendulum</em>
     * iterations, where only one direction of edges is considered, then <em>rubber</em> iterations,
     * where both incoming and outgoing edges are considered. In each iteration first determine the
     * <em>deflection</em> of each linear segment, i.e. the optimal position delta that leads to
     * a balanced placement with respect to its adjacent segments. Then merge regions that touch each
     * other, building mean values of the involved deflections, and finally apply the resulting
     * deflection values to all segments. The iterations stop when no further improvement is done.
     *
     * @param layeredGraph a layered graph
     */
    private void balancePlacement(final LGraph layeredGraph) {
        double deflectionDampening = layeredGraph.getProperty(
                LayeredOptions.NODE_PLACEMENT_LINEAR_SEGMENTS_DEFLECTION_DAMPENING).doubleValue();

        // Determine a suitable number of pendulum iterations
        int thoroughness = layeredGraph.getProperty(LayeredOptions.THOROUGHNESS);
        int pendulumIters = PENDULUM_ITERS;
        int finalIters = FINAL_ITERS;
        double threshold = THRESHOLD_FACTOR / thoroughness;

        // Iterate the balancing
        boolean ready = false;
        Mode mode = Mode.FORW_PENDULUM;
        double lastTotalDeflection = Integer.MAX_VALUE;
        do {

            // Calculate force for every linear segment
            boolean incoming = mode != Mode.BACKW_PENDULUM;
            boolean outgoing = mode != Mode.FORW_PENDULUM;
            double totalDeflection = 0;
            for (LinearSegment segment : linearSegments) {
                segment.refSegment = null;
                calcDeflection(segment, incoming, outgoing, deflectionDampening);
                totalDeflection += Math.abs(segment.deflection);
            }

            // Merge linear segments to form regions
            boolean merged;
            do {
                merged = mergeRegions(layeredGraph);
            } while (merged);

            // Move the nodes according to the deflection value of their region
            for (LinearSegment segment : linearSegments) {
                double deflection = segment.region().deflection;
                if (deflection != 0) {
                    for (LNode node : segment.nodes) {
                        node.getPosition().y += deflection;
                    }
                }
            }

            // Update the balancing mode
            if (mode == Mode.FORW_PENDULUM || mode == Mode.BACKW_PENDULUM) {
                pendulumIters--;
                if (pendulumIters <= 0 && (totalDeflection < lastTotalDeflection
                        || -pendulumIters > thoroughness)) {
                    mode = Mode.RUBBER;
                    lastTotalDeflection = Integer.MAX_VALUE;
                } else if (mode == Mode.FORW_PENDULUM) {
                    mode = Mode.BACKW_PENDULUM;
                    lastTotalDeflection = totalDeflection;
                } else {
                    mode = Mode.FORW_PENDULUM;
                    lastTotalDeflection = totalDeflection;
                }
            } else {
                ready = totalDeflection >= lastTotalDeflection
                        || lastTotalDeflection - totalDeflection < threshold;
                lastTotalDeflection = totalDeflection;
                if (ready) {
                    finalIters--;
                }
            }
        } while (!(ready && finalIters <= 0));
    }

    /**
     * Calculate the force acting on the given linear segment. The force is stored in the segment's
     * deflection field.
     *
     * @param segment the linear segment whose force is to be calculated
     * @param incoming whether incoming edges should be considered
     * @param outgoing whether outgoing edges should be considered
     * @param deflectionDampening factor by which deflections are dampened
     */
    private void calcDeflection(final LinearSegment segment, final boolean incoming,
            final boolean outgoing, final double deflectionDampening) {

        double segmentDeflection = 0;
        int nodeWeightSum = 0;
        for (LNode node : segment.nodes) {
            double nodeDeflection = 0;
            int edgeWeightSum = 0;
            int inputPrio = incoming ? node.getProperty(INPUT_PRIO) : Integer.MIN_VALUE;
            int outputPrio = outgoing ? node.getProperty(OUTPUT_PRIO) : Integer.MIN_VALUE;
            int minPrio = Math.max(inputPrio, outputPrio);

            // Calculate force for every port/edge
            for (LPort port : node.getPorts()) {
                double portpos = node.getPosition().y + port.getPosition().y + port.getAnchor().y;
                if (outgoing) {
                    for (LEdge edge : port.getOutgoingEdges()) {
                        LPort otherPort = edge.getTarget();
                        LNode otherNode = otherPort.getNode();
                        if (segment != linearSegments[otherNode.id]) {
                            int otherPrio = Math.max(otherNode.getProperty(INPUT_PRIO),
                                    otherNode.getProperty(OUTPUT_PRIO));
                            int prio = edge.getProperty(LayeredOptions.PRIORITY_STRAIGHTNESS);
                            if (prio >= minPrio && prio >= otherPrio) {
                                nodeDeflection += otherNode.getPosition().y
                                        + otherPort.getPosition().y + otherPort.getAnchor().y
                                        - portpos;
                                edgeWeightSum++;
                            }
                        }
                    }
                }

                if (incoming) {
                    for (LEdge edge : port.getIncomingEdges()) {
                        LPort otherPort = edge.getSource();
                        LNode otherNode = otherPort.getNode();
                        if (segment != linearSegments[otherNode.id]) {
                            int otherPrio = Math.max(otherNode.getProperty(INPUT_PRIO),
                                    otherNode.getProperty(OUTPUT_PRIO));
                            int prio = edge.getProperty(LayeredOptions.PRIORITY_STRAIGHTNESS);
                            if (prio >= minPrio && prio >= otherPrio) {
                                nodeDeflection += otherNode.getPosition().y
                                        + otherPort.getPosition().y + otherPort.getAnchor().y
                                        - portpos;
                                edgeWeightSum++;
                            }
                        }
                    }
                }
            }

            // Avoid division by zero
            if (edgeWeightSum > 0) {
                segmentDeflection += nodeDeflection / edgeWeightSum;
                nodeWeightSum++;
            }
        }
        if (nodeWeightSum > 0) {
            segment.deflection = deflectionDampening * segmentDeflection / nodeWeightSum;
            segment.weight = nodeWeightSum;
        } else {
            segment.deflection = 0;
            segment.weight = 0;
        }
    }

    /**
     * Factor for threshold within which node overlapping is detected. This factor should reflect
     * the epsilon used in the {@link NodePlacerTest}.
     */
    private static final double OVERLAP_DETECT = 0.0001;

    /**
     * Merge regions by testing whether they would overlap after applying the deflection.
     *
     * @param layeredGraph the layered graph
     * @return true if any two regions have been merged
     */
    private boolean mergeRegions(final LGraph layeredGraph) {

        boolean changed = false;
        double nodeSpacing = layeredGraph.getProperty(LayeredOptions.SPACING_NODE_NODE).doubleValue();
        double threshold = OVERLAP_DETECT * nodeSpacing;
        for (Layer layer : layeredGraph) {
            Iterator<LNode> nodeIter = layer.getNodes().iterator();

            // Get the first node
            LNode node1 = nodeIter.next();
            LinearSegment region1 = linearSegments[node1.id].region();

            // While there are still nodes following the current node
            while (nodeIter.hasNext()) {
                // Test whether nodes have different regions
                LNode node2 = nodeIter.next();
                LinearSegment region2 = linearSegments[node2.id].region();

                if (region1 != region2) {
                    // Calculate how much space is allowed between the nodes
                    double spacing = spacings.getVerticalSpacing(node1, node2);

                    double node1Extent = node1.getPosition().y + node1.getSize().y
                            + node1.getMargin().bottom + region1.deflection + spacing;
                    double node2Extent = node2.getPosition().y - node2.getMargin().top
                            + region2.deflection;

                    // Test if the nodes are overlapping
                    if (node1Extent > node2Extent + threshold) {
                        // Merge the first region under the second top level segment
                        int weightSum = region1.weight + region2.weight;
                        assert weightSum > 0;
                        region2.deflection = (region2.weight * region2.deflection
                                + region1.weight * region1.deflection) / weightSum;
                        region2.weight = weightSum;
                        region1.refSegment = region2;
                        changed = true;
                    }
                }

                node1 = node2;
                region1 = region2;
            }
        }
        return changed;
    }


    // /////////////////////////////////////////////////////////////////////////////
    // Post Processing for Correction

    /**
     * Post-process the balanced placement by moving linear segments where obvious improvements can
     * be made.
     *
     * @param layeredGraph
     *            the layered graph
     */
    private void postProcess(final LGraph layeredGraph) {

        // process each linear segment independently
        for (LinearSegment segment : linearSegments) {
            double minRoomAbove = Integer.MAX_VALUE, minRoomBelow = Integer.MAX_VALUE;

            for (LNode node : segment.nodes) {
                double roomAbove, roomBelow;
                int index = node.getIndex();

                // determine the amount by which the linear segment can be moved up without overlap
                if (index > 0) {
                    LNode neighbor = node.getLayer().getNodes().get(index - 1);
                    double spacing = spacings.getVerticalSpacing(node, neighbor);
                    roomAbove = node.getPosition().y - node.getMargin().top
                            - (neighbor.getPosition().y + neighbor.getSize().y
                                    + neighbor.getMargin().bottom + spacing);
                } else {
                    roomAbove = node.getPosition().y - node.getMargin().top;
                }
                minRoomAbove = Math.min(roomAbove, minRoomAbove);

                // determine the amount by which the linear segment can be moved down without
                // overlap
                if (index < node.getLayer().getNodes().size() - 1) {
                    LNode neighbor = node.getLayer().getNodes().get(index + 1);
                    double spacing = spacings.getVerticalSpacing(node, neighbor);
                    roomBelow = neighbor.getPosition().y - neighbor.getMargin().top
                            - (node.getPosition().y + node.getSize().y
                                    + node.getMargin().bottom + spacing);
                } else {
                    roomBelow = 2 * node.getPosition().y;
                }
                minRoomBelow = Math.min(roomBelow, minRoomBelow);
            }

            double minDisplacement = Integer.MAX_VALUE;
            boolean foundPlace = false;

            // determine the minimal displacement that would make one incoming edge straight
            LNode firstNode = segment.nodes.get(0);
            for (LPort target : firstNode.getPorts()) {
                double pos = firstNode.getPosition().y + target.getPosition().y + target.getAnchor().y;
                for (LEdge edge : target.getIncomingEdges()) {
                    LPort source = edge.getSource();
                    double d = source.getNode().getPosition().y + source.getPosition().y
                            + source.getAnchor().y - pos;
                    if (Math.abs(d) < Math.abs(minDisplacement)
                            && Math.abs(d) < (d < 0 ? minRoomAbove : minRoomBelow)) {
                        minDisplacement = d;
                        foundPlace = true;
                    }
                }
            }

            // determine the minimal displacement that would make one outgoing edge straight
            LNode lastNode = segment.nodes.get(segment.nodes.size() - 1);
            for (LPort source : lastNode.getPorts()) {
                double pos = lastNode.getPosition().y + source.getPosition().y + source.getAnchor().y;
                for (LEdge edge : source.getOutgoingEdges()) {
                    LPort target = edge.getTarget();
                    double d = target.getNode().getPosition().y + target.getPosition().y
                            + target.getAnchor().y - pos;
                    if (Math.abs(d) < Math.abs(minDisplacement)
                            && Math.abs(d) < (d < 0 ? minRoomAbove : minRoomBelow)) {
                        minDisplacement = d;
                        foundPlace = true;
                    }
                }
            }

            // if such a displacement could be found, apply it to the whole linear segment
            if (foundPlace && minDisplacement != 0) {
                for (LNode node : segment.nodes) {
                    node.getPosition().y += minDisplacement;
                }
            }
        }
    }

}
