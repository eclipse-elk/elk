/*******************************************************************************
 * Copyright (c) 2016 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.alg.layered.p3order;

import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.Random;
import java.util.Set;

import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.intermediate.greedyswitch.GreedySwitchHeuristic;
import org.eclipse.elk.alg.layered.options.GraphProperties;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.layered.p3order.LayerSweepCrossingMinimizer.CrossMinType;
import org.eclipse.elk.alg.layered.p3order.counting.AllCrossingsCounter;
import org.eclipse.elk.alg.layered.p3order.counting.IInitializable;

import com.google.common.collect.Lists;

/**
 * Collects data needed for cross minimization and port distribution.
 * <p>
 * Must be initialized using {@link IInitializable#init(java.util.List)}!
 * </p>
 * 
 * @author alan
 */
public class GraphInfoHolder implements IInitializable {
    /** Raw graph data. */
    private LGraph lGraph;

    /** Saved node orders. */
    private LNode[][] currentNodeOrder;
    private SweepCopy currentlyBestNodeAndPortOrder;
    private SweepCopy bestNodeAndPortOrder;
    private int[] portPositions;

    /** Processing type information. */
    private boolean useBottomUp;

    /** Hierarchy access. */
    private List<LGraph> childGraphs;
    private boolean hasExternalPorts;
    private boolean hasParent;
    private GraphInfoHolder parentGraphData;
    private LNode parent;
    private LayerSweepTypeDecider layerSweepTypeDecider;

    /** Pre-initialized auxiliary objects. */
    private ICrossingMinimizationHeuristic crossMinimizer;
    private ISweepPortDistributor portDistributor;
    private AllCrossingsCounter crossingsCounter;
    private int nPorts;

    /**
     * Create object collecting information about a graph.
     * 
     * @param graph
     *            The graph
     * @param crossMinType
     *            The CrossMinimizer
     * @param graphs
     *            the complete list of all graphs in the hierarchy
     */
    public GraphInfoHolder(final LGraph graph, final CrossMinType crossMinType, final List<GraphInfoHolder> graphs) {
        lGraph = graph;
        currentNodeOrder = graph.toNodeArray();

        // Hierarchy information.
        parent = lGraph.getParentNode();
        hasParent = parent != null;
        parentGraphData = hasParent ? graphs.get(parent.getGraph().id) : null;
        Set<GraphProperties> graphProperties = graph.getProperty(InternalProperties.GRAPH_PROPERTIES);
        hasExternalPorts = graphProperties.contains(GraphProperties.EXTERNAL_PORTS);
        childGraphs = Lists.newArrayList();

        // Init all objects needing initialization by graph traversal.
        crossingsCounter = new AllCrossingsCounter(currentNodeOrder);
        Random random = lGraph.getProperty(InternalProperties.RANDOM);
        portDistributor = ISweepPortDistributor.create(crossMinType, random, currentNodeOrder);
        layerSweepTypeDecider = new LayerSweepTypeDecider(this);
        List<IInitializable> initializables =
                Lists.newArrayList(this, crossingsCounter, layerSweepTypeDecider, portDistributor);

        if (crossMinType == CrossMinType.BARYCENTER
                && !graph.getProperty(LayeredOptions.CROSSING_MINIMIZATION_FORCE_NODE_MODEL_ORDER)) {
            ForsterConstraintResolver constraintResolver = new ForsterConstraintResolver(currentNodeOrder);
            initializables.add(constraintResolver);
            crossMinimizer = new BarycenterHeuristic(constraintResolver, random,
                    (AbstractBarycenterPortDistributor) portDistributor, currentNodeOrder);
        } else if (crossMinType == CrossMinType.BARYCENTER
                && graph.getProperty(LayeredOptions.CROSSING_MINIMIZATION_FORCE_NODE_MODEL_ORDER)) {
            ForsterConstraintResolver constraintResolver = new ForsterConstraintResolver(currentNodeOrder);
            initializables.add(constraintResolver);
            crossMinimizer = new ModelOrderBarycenterHeuristic(constraintResolver, random,
                    (AbstractBarycenterPortDistributor) portDistributor, currentNodeOrder);
        } else if (crossMinType == CrossMinType.MEDIAN) {
            crossMinimizer = new MedianHeuristic(random);
        } else {
            crossMinimizer = new GreedySwitchHeuristic(crossMinType, this);
        }
        initializables.add(crossMinimizer);

        // Apply Initializer.
        IInitializable.init(initializables, currentNodeOrder);

        // calculate whether we need to use bottom up or sweep into this graph.
        useBottomUp = layerSweepTypeDecider.useBottomUp();
    }

    /**
     * @return the processRecursively
     */
    public boolean dontSweepInto() {
        return useBottomUp;
    }

    /**
     * @return the lGraph
     */
    public LGraph lGraph() {
        return lGraph;
    }

    /**
     * @return the currentNodeOrder
     */
    public LNode[][] currentNodeOrder() {
        return currentNodeOrder;
    }

    /**
     * @return the currentlyBestNodeAndPortOrder
     */
    public SweepCopy currentlyBestNodeAndPortOrder() {
        return currentlyBestNodeAndPortOrder;
    }

    /**
     * @param currentlyBestNodeAndPortOrder
     *            the currentlyBestNodeAndPortOrder to set
     */
    public void setCurrentlyBestNodeAndPortOrder(final SweepCopy currentlyBestNodeAndPortOrder) {
        this.currentlyBestNodeAndPortOrder = currentlyBestNodeAndPortOrder;
    }

    /**
     * @return the bestNodeNPortOrder
     */
    public SweepCopy bestNodeNPortOrder() {
        return bestNodeAndPortOrder;
    }

    /**
     * @param bestNodeNPortOrder
     *            the bestNodeNPortOrder to set
     */
    public void setBestNodeNPortOrder(final SweepCopy bestNodeNPortOrder) {
        this.bestNodeAndPortOrder = bestNodeNPortOrder;
    }

    /**
     * @return the crossCounter
     */
    public AllCrossingsCounter crossCounter() {
        return crossingsCounter;
    }

    /**
     * @return the crossMinimizer
     */
    public ICrossingMinimizationHeuristic crossMinimizer() {
        return crossMinimizer;
    }

    /**
     * @return the portDistributor
     */
    public ISweepPortDistributor portDistributor() {
        return portDistributor;
    }

    /**
     * @return the parent
     */
    public LNode parent() {
        return parent;
    }

    /**
     * @return the hasParent
     */
    public boolean hasParent() {
        return hasParent;
    }

    /**
     * @return child graphs
     */
    public Collection<LGraph> childGraphs() {
        return childGraphs;
    }

    /**
     * @return whether this graph's parent node has external ports.
     */
    public boolean hasExternalPorts() {
        return hasExternalPorts;
    }

    @Override
    public String toString() {
        return Arrays.deepToString(currentNodeOrder);
    }

    /**
     * @return Copy of node order for currently best sweep.
     */
    public SweepCopy getBestSweep() {
        return crossMinimizer.isDeterministic() ? currentlyBestNodeAndPortOrder() : bestNodeNPortOrder();
    }

    /**
     * @return the parentGraphData
     */
    public GraphInfoHolder parentGraphData() {
        return parentGraphData;
    }

    /**
     * @return
     */
    public boolean crossMinDeterministic() {
        return crossMinimizer.isDeterministic();
    }

    /**
     * @return whether this CrossingMinimizer always improves
     */
    public boolean crossMinAlwaysImproves() {
        return crossMinimizer.alwaysImproves();
    }

    /**
     * Return port position array of the necessary length. Each algorithm must initialize the actual positions
     * separately.
     * 
     * @return port position array
     */
    public int[] portPositions() {
        return portPositions;
    }

    @Override
    public void initAtNodeLevel(final int l, final int n, final LNode[][] nodeOrder) {
        LNode node = nodeOrder[l][n];
        LGraph nestedGraph = node.getNestedGraph();
        if (nestedGraph != null) {
            childGraphs.add(nestedGraph);
        }
    }

    @Override
    public void initAtPortLevel(final int l, final int n, final int p, final LNode[][] nodeOrder) {
        nPorts++;
    }

    @Override
    public void initAfterTraversal() {
        portPositions = new int[nPorts];
    }

}
