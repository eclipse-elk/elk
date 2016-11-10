/*******************************************************************************
 * Copyright (c) 2010, 2015 Kiel University and others.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the Eclipse Public License v1.0
 * which accompanies this distribution, and is available at
 * http://www.eclipse.org/legal/epl-v10.html
 *
 * Contributors:
 *     Kiel University - initial API and implementation
 *******************************************************************************/
package org.eclipse.elk.alg.layered.p3order;

import java.util.ArrayDeque;
import java.util.Deque;
import java.util.List;
import java.util.Random;
import java.util.Set;

import org.eclipse.elk.alg.layered.ILayoutPhase;
import org.eclipse.elk.alg.layered.IntermediateProcessingConfiguration;
import org.eclipse.elk.alg.layered.IntermediateProcessingConfiguration.Slot;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LNode.NodeType;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.intermediate.IntermediateProcessorStrategy;
import org.eclipse.elk.alg.layered.p3order.counting.CrossMinUtil;
import org.eclipse.elk.alg.layered.properties.InternalProperties;
import org.eclipse.elk.alg.layered.properties.LayeredOptions;
import org.eclipse.elk.core.options.PortConstraints;
import org.eclipse.elk.core.options.PortSide;
import org.eclipse.elk.core.util.IElkProgressMonitor;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;

/**
 * This class minimizes crossings by sweeping through a graph, holding the order of nodes in one layer fixed and
 * switching the nodes in the other layer. After each re-sorting step, the ports in the two current layers are sorted.
 * 
 * By using a {@link LayerSweepTypeDecider} as called by the {@link GraphInfoHolder} initializing class, each graph can
 * either be dealt with <i>bottom-up</i> or <i>hierarchically</i>.
 * 
 * <ul>
 * <li>bottom-up: the nodes of a child graph are sorted, then the order of ports of hierarchy border traversing edges
 * are fixed. Then the parent graph is laid out, viewing the child graph as an atomic node.</li>
 * <li>hierarchical: When reaching a node with a child graph marked as hierarchical: First the ports are sorted on the
 * outside of the graph. Then the nodes on the inside are sorted by the order of the ports. Then the child graph is
 * swept through. Then the ports of the parent node are sorted by the order of the nodes on the last layer of the child
 * graph. Finally the sweep through the parent graph is continued.</li>
 * </ul>
 * 
 * Therefore this is a <i>hierarchical</i> processor which must have access to the root graph.
 * <p>
 * Reference for the original layer sweep:
 * <ul>
 * <li>Kozo Sugiyama, Shojiro Tagawa, and Mitsuhiko Toda. Methods for visual understanding of hierarchical system
 * structures. IEEE Transactions on Systems, Man and Cybernetics, 11(2):109–125, February 1981.
 * </ul>
 * </p>
 * <dl>
 * <dt>Precondition:</dt>
 * <dd>The graph has a proper layering, i.e. all long edges have been splitted; all nodes have at least fixed port
 * sides.</dd>
 * <dt>Postcondition:</dt>
 * <dd>The order of nodes in each layer and the order of ports in each node are optimized to yield as few edge crossings
 * as possible</dd>
 * </dl>
 *
 * @author alan
 *
 */
public class LayerSweepCrossingMinimizer implements ILayoutPhase {
    /** Collected information about each graph. */
    private List<GraphInfoHolder> graphInfoHolders;
    /** We only need to save the orders of graphs if their node order actually changed. */
    private Set<GraphInfoHolder> graphsWhoseNodeOrderChanged;
    private Random random;
    private long randomSeed;
    private final CrossMinType crossMinType;

    /**
     * Creates LayerSweepHierarchicalCrossingMinimizer using given minimizer type.
     *
     * @param cT
     *            the crossing minimizer 'type'
     */
    public LayerSweepCrossingMinimizer(final CrossMinType cT) {
        crossMinType = cT;
    }

    @Override
    public boolean operatesOnFullHierarchy() {
        return true;
    }

    @Override
    public void process(final LGraph layeredGraph, final IElkProgressMonitor progressMonitor) {
        progressMonitor.begin("Minimize Crossings " + crossMinType, 1);

        if (layeredGraph.getLayers().isEmpty()) {
            progressMonitor.done();
            return;
        }

        List<GraphInfoHolder> graphsToSweepOn = initialize(layeredGraph);

        Consumer<GraphInfoHolder> minimizingMethod = chooseMinimizingMethod(graphsToSweepOn);

        minimizeCrossings(graphsToSweepOn, minimizingMethod);

        transferNodeAndPortOrdersToGraph();

        progressMonitor.done();
    }

    private Consumer<GraphInfoHolder> chooseMinimizingMethod(final List<GraphInfoHolder> graphsToSweepOn) {
        GraphInfoHolder parent = graphsToSweepOn.get(0);
        if (!parent.crossMinDeterministic()) {
            return g -> compareDifferentRandomizedLayouts(g);
        } else if (parent.crossMinAlwaysImproves()) {
            return g -> minimizeCrossingsNoCounter(g);
        } else {
            return g -> minimizeCrossingsWithCounter(g);
        }
    }

    private void minimizeCrossings(final List<GraphInfoHolder> graphsToSweepOn,
            final Consumer<GraphInfoHolder> minimizingMethod) {
        for (GraphInfoHolder gData : graphsToSweepOn) {
            minimizingMethod.accept(gData);
            if (gData.hasParent()) {
                setPortOrderOnParentGraph(gData);
            }
        }
    }

    private void setPortOrderOnParentGraph(final GraphInfoHolder gData) {
        if (gData.hasExternalPorts()) {
            SweepCopy bestSweep = gData.getBestSweep();
            // Sort ports on left and right side of the parent node
            sortPortsByDummyPositionsInLastLayer(bestSweep.nodes(), gData.parent(), true);
            sortPortsByDummyPositionsInLastLayer(bestSweep.nodes(), gData.parent(), false);
            gData.parent().setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_ORDER);
        }
    }

    /** For use with any two-layer crossing minimizer which always improves crossings (e.g. two-sided greedy switch). */
    private void minimizeCrossingsNoCounter(final GraphInfoHolder gData) {
        boolean isForwardSweep = random.nextBoolean();
        boolean improved = true;
        while (improved) {
            improved = false;
            improved = gData.crossMinimizer().setFirstLayerOrder(gData.currentNodeOrder(), isForwardSweep);
            improved |= sweepReducingCrossings(gData, isForwardSweep, false);
            isForwardSweep = !isForwardSweep;
        }
        setCurrentlyBestNodeOrders();
    }

    private void compareDifferentRandomizedLayouts(final GraphInfoHolder gData) {
        // Reset the seed, otherwise copies of hierarchical graphs in different parent nodes are
        // layouted differently.
        random.setSeed(randomSeed);

        // In order to only copy graphs whose node order has changed, save them in a set.
        graphsWhoseNodeOrderChanged.clear();

        int bestCrossings = Integer.MAX_VALUE;
        int thouroughness = gData.lGraph().getProperty(LayeredOptions.THOROUGHNESS);
        for (int i = 0; i < thouroughness; i++) {
            int crossings = minimizeCrossingsWithCounter(gData);
            if (crossings < bestCrossings) {
                bestCrossings = crossings;
                saveAllNodeOrdersOfChangedGraphs();
                if (bestCrossings == 0) {
                    break;
                }
            }
        }
    }

    private int minimizeCrossingsWithCounter(final GraphInfoHolder gData) {
        boolean isForwardSweep = random.nextBoolean();

        gData.crossMinimizer().setFirstLayerOrder(gData.currentNodeOrder(), isForwardSweep);
        sweepReducingCrossings(gData, isForwardSweep, true);
        int crossingsInGraph = countCurrentNumberOfCrossings(gData);
        int oldNumberOfCrossings;
        do {
            setCurrentlyBestNodeOrders();

            if (crossingsInGraph == 0) {
                return 0;
            }

            isForwardSweep = !isForwardSweep;
            oldNumberOfCrossings = crossingsInGraph;
            sweepReducingCrossings(gData, isForwardSweep, false);
            crossingsInGraph = countCurrentNumberOfCrossings(gData);
        } while (oldNumberOfCrossings > crossingsInGraph);

        return oldNumberOfCrossings;
    }

    /*
     * We only need to count crossings below the current graph and also only if they are marked as to be processed
     * hierarchically.
     */
    private int countCurrentNumberOfCrossings(final GraphInfoHolder currentGraph) {
        int totalCrossings = 0;
        Deque<GraphInfoHolder> countCrossingsIn = new ArrayDeque<>();
        countCrossingsIn.push(currentGraph);
        while (!countCrossingsIn.isEmpty()) {
            GraphInfoHolder gD = countCrossingsIn.pop();
            totalCrossings +=
                    gD.crossCounter().countAllCrossings(gD.currentNodeOrder());
            for (LGraph childLGraph : gD.childGraphs()) {
                GraphInfoHolder child = graphInfoHolders.get(childLGraph.id);
                if (!child.dontSweepInto()) {
                    totalCrossings += countCurrentNumberOfCrossings(child);
                }
            }
        }

        return totalCrossings;
    }

    private boolean sweepReducingCrossings(final GraphInfoHolder graph, final boolean forward,
            final boolean firstSweep) {
        LNode[][] nodes = graph.currentNodeOrder();
        int length = nodes.length;

        boolean improved =
                graph.portDistributor().distributePortsWhileSweeping(nodes, firstIndex(forward, length), forward);
        LNode[] firstLayer = nodes[firstIndex(forward, length)];
        improved |= sweepInHierarchicalNodes(firstLayer, forward, firstSweep);
        for (int i = firstFree(forward, length); isNotEnd(length, i, forward); i += next(forward)) {
            improved |= graph.crossMinimizer().minimizeCrossings(nodes, i, forward, firstSweep);
            improved |= graph.portDistributor().distributePortsWhileSweeping(nodes, i, forward);
            improved |= sweepInHierarchicalNodes(nodes[i], forward, firstSweep);
        }

        graphsWhoseNodeOrderChanged.add(graph);
        return improved;
    }

    private boolean sweepInHierarchicalNodes(final LNode[] layer, final boolean isForwardSweep,
            final boolean isFirstSweep) {
        boolean improved = false;
        for (LNode node : layer) {
            if (hasNestedGraph(node) && !graphInfoHolders.get(nestedGraphOf(node).id).dontSweepInto()) {
                improved |= sweepInHierarchicalNode(isForwardSweep, node, isFirstSweep);
            }
        }
        return improved;
    }

    private boolean sweepInHierarchicalNode(final boolean isForwardSweep, final LNode node,
            final boolean isFirstSweep) {
        LGraph nestedLGraph = nestedGraphOf(node);
        GraphInfoHolder nestedGraph = graphInfoHolders.get(nestedLGraph.id);
        LNode[][] nestedGraphNodeOrder = nestedGraph.currentNodeOrder();
        int startIndex = firstIndex(isForwardSweep, nestedGraphNodeOrder.length);
        LNode firstNode = nestedGraphNodeOrder[startIndex][0];

        if (isExternalPortDummy(firstNode)) {
            nestedGraphNodeOrder[startIndex] = sortPortDummiesByPortPositions(node,
                    nestedGraphNodeOrder[startIndex], sideOpposedSweepDirection(isForwardSweep));
        } else {
            nestedGraph.crossMinimizer().setFirstLayerOrder(nestedGraphNodeOrder, isForwardSweep);
        }

        boolean improved = sweepReducingCrossings(nestedGraph, isForwardSweep, isFirstSweep);
        sortPortsByDummyPositionsInLastLayer(nestedGraph.currentNodeOrder(), nestedGraph.parent(),
                isForwardSweep);

        return improved;
    }

    private void sortPortsByDummyPositionsInLastLayer(final LNode[][] nodeOrder, final LNode parent,
            final boolean onRightMostLayer) {
        int endIndex = endIndex(onRightMostLayer, nodeOrder.length);
        LNode[] lastLayer = nodeOrder[endIndex];
        if (!isExternalPortDummy(lastLayer[0])) {
            return;
        }

        int j = firstIndex(onRightMostLayer, lastLayer.length);
        List<LPort> ports = parent.getPorts();
        for (int i = 0; i < ports.size(); i++) {
            LPort port = ports.get(i);
            if (isOnEndOfSweepSide(port, onRightMostLayer) && isHierarchical(port)) {
                ports.set(i, originPort(lastLayer[j]));
                j += next(onRightMostLayer);
            }
        }
    }

    private LNode[] sortPortDummiesByPortPositions(final LNode parentNode,
            final LNode[] layerCloseToNodeEdge, final PortSide side) {
        Iterable<LPort> ports = CrossMinUtil.inNorthSouthEastWestOrder(parentNode, side);

        LNode[] sortedDummies = new LNode[layerCloseToNodeEdge.length];

        int i = 0;
        for (LPort port : ports) {
            if (isHierarchical(port)) {
                sortedDummies[i++] = dummyNodeFor(port);
            }
        }

        return sortedDummies;
    }

    private void saveAllNodeOrdersOfChangedGraphs() {
        for (GraphInfoHolder graph : graphsWhoseNodeOrderChanged) {
            graph.setBestNodeNPortOrder(new SweepCopy(graph.currentlyBestNodeAndPortOrder()));
        }
    }

    private void setCurrentlyBestNodeOrders() {
        for (GraphInfoHolder graph : graphsWhoseNodeOrderChanged) {
            graph.setCurrentlyBestNodeAndPortOrder(new SweepCopy(graph.currentNodeOrder()));
        }
    }

    private int firstIndex(final boolean isForwardSweep, final int length) {
        return isForwardSweep ? 0 : length - 1;
    }

    private int endIndex(final boolean isForwardSweep, final int length) {
        return isForwardSweep ? length - 1 : 0;
    }

    private int firstFree(final boolean isForwardSweep, final int length) {
        return isForwardSweep ? 1 : length - 2;
    }

    private int next(final boolean isForwardSweep) {
        return isForwardSweep ? 1 : -1;
    }

    private boolean isNotEnd(final int length, final int freeLayerIndex,
            final boolean isForwardSweep) {
        return isForwardSweep ? freeLayerIndex < length : freeLayerIndex >= 0;
    }

    private LGraph nestedGraphOf(final LNode node) {
        return node.getProperty(InternalProperties.NESTED_LGRAPH);
    }

    private Boolean hasNestedGraph(final LNode node) {
        return nestedGraphOf(node) != null;
    }

    private PortSide sideOpposedSweepDirection(final boolean isForwardSweep) {
        return isForwardSweep ? PortSide.WEST : PortSide.EAST;
    }

    private boolean isExternalPortDummy(final LNode firstNode) {
        return firstNode.getType() == NodeType.EXTERNAL_PORT;
    }

    private LPort originPort(final LNode node) {
        return (LPort) node.getProperty(InternalProperties.ORIGIN);
    }

    private boolean isHierarchical(final LPort port) {
        return port.getProperty(InternalProperties.INSIDE_CONNECTIONS);
    }

    private LNode dummyNodeFor(final LPort port) {
        return port.getProperty(InternalProperties.PORT_DUMMY);
    }

    private boolean isOnEndOfSweepSide(final LPort port, final boolean isForwardSweep) {
        return isForwardSweep ? port.getSide() == PortSide.EAST : port.getSide() == PortSide.WEST;
    }

    /**
     * Traverses inclusion breadth-first and initializes each Graph.
     */
    private List<GraphInfoHolder> initialize(final LGraph rootGraph) {
        graphInfoHolders = Lists.newArrayList();
        random = rootGraph.getProperty(InternalProperties.RANDOM);
        randomSeed = random.nextLong();
        List<GraphInfoHolder> graphsToSweepOn = Lists.newLinkedList();
        List<LGraph> graphs = Lists.<LGraph>newArrayList(rootGraph);
        int i = 0;
        while (i < graphs.size()) {
            LGraph graph = graphs.get(i);
            graph.id = i++;
            GraphInfoHolder gData = new GraphInfoHolder(graph, crossMinType, graphInfoHolders);
            graphs.addAll(gData.childGraphs());
            graphInfoHolders.add(gData);
            if (gData.dontSweepInto()) {
                graphsToSweepOn.add(0, gData);
            }
        }
    
        graphsWhoseNodeOrderChanged = Sets.newHashSet();
    
        return graphsToSweepOn;
    }

    private void transferNodeAndPortOrdersToGraph() {
        for (GraphInfoHolder gD : graphInfoHolders) {
            gD.getBestSweep().transferNodeAndPortOrdersToGraph(gD.lGraph());
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public IntermediateProcessingConfiguration getIntermediateProcessingConfiguration(
            final LGraph graph) {
        IntermediateProcessingConfiguration configuration = IntermediateProcessingConfiguration
                .fromExisting(INTERMEDIATE_PROCESSING_CONFIGURATION);
        configuration.addBeforePhase3(IntermediateProcessorStrategy.PORT_LIST_SORTER);

        return configuration;
    }

    /**
     * @return the graphData
     */
    public List<GraphInfoHolder> getGraphData() {
        return graphInfoHolders;
    }

    /**
     * Type of crossing minimizer.
     * 
     */
    public enum CrossMinType {
        /** Use BarycenterHeuristic. */
        BARYCENTER,
        /** Use one-sided GreedySwitchHeuristic. */
        ONE_SIDED_GREEDY_SWITCH,
        /** Use two-sided GreedySwitchHeuristic. */
        TWO_SIDED_GREEDY_SWITCH
    }

    /** intermediate processing configuration. */
    private static final IntermediateProcessingConfiguration INTERMEDIATE_PROCESSING_CONFIGURATION =
            IntermediateProcessingConfiguration.createEmpty()
                    .addBeforePhase3(IntermediateProcessorStrategy.LONG_EDGE_SPLITTER)
                    .addBeforePhase4(IntermediateProcessorStrategy.IN_LAYER_CONSTRAINT_PROCESSOR)
                    .addAll(Slot.AFTER_PHASE_5,
                            Lists.newArrayList(IntermediateProcessorStrategy.LONG_EDGE_JOINER,
                                    IntermediateProcessorStrategy.HIERARCHICAL_NODE_RESIZER));
    
    /**
     * TODO Java 8 api not supported yet and Guava's Void is much much uglier.
     * 
     * @param <T>
     */
    private interface Consumer<T> {
        void accept(T input);
    }
}
