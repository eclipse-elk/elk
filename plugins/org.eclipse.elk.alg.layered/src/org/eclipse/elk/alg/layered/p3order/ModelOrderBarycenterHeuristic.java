/*******************************************************************************
 * Copyright (c) 2021 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 * 
 * SPDX-License-Identifier: EPL-2.0 
 *******************************************************************************/
package org.eclipse.elk.alg.layered.p3order;

import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Random;

import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.intermediate.preserveorder.CMGroupModelOrderCalculator;
import org.eclipse.elk.alg.layered.options.GroupOrderStrategy;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.LayerConstraint;
import org.eclipse.elk.alg.layered.options.LayeredOptions;

/**
 * Minimizes crossing with the barycenter method. However, the node order given by the order in the model
 * does not change. I.e. only the dummy nodes are sorted in the already sorted real nodes.
 */
public class ModelOrderBarycenterHeuristic extends BarycenterHeuristic {
    
    /**
     * Each node has an entry of nodes for which it is bigger.
     */
    private HashMap<LNode, HashSet<LNode>> biggerThan = new HashMap<>();
    /**
     * Each node has an entry of nodes for which it is smaller.
     */
    private HashMap<LNode, HashSet<LNode>> smallerThan = new HashMap<>();

    /**
     * Constructs a model order barycenter heuristic for crossing minimization.
     * After sorting {@link #clearTransitiveOrdering()} should be called. 
     * 
     * @param constraintResolver
     *            the constraint resolver
     * @param random
     *            the random number generator
     * @param portDistributor
     *            calculates the port ranks for the barycenter heuristic.
     * @param graph
     *            current node order
     */
    public ModelOrderBarycenterHeuristic(final ForsterConstraintResolver constraintResolver, final Random random,
            final AbstractBarycenterPortDistributor portDistributor, final LNode[][] graph) {
        super(constraintResolver, random, portDistributor, graph);
        // This may have the problem that the real nodes need to be compared first such that the dummy nodes do not
        // overrule model order using the barycenter method.
        barycenterStateComparator = 
                (n1, n2) -> {
                    if (n1.hasProperty(LayeredOptions.LAYERING_LAYER_CONSTRAINT) && (n1.getProperty(LayeredOptions.LAYERING_LAYER_CONSTRAINT) == LayerConstraint.FIRST_SEPARATE
                            || n1.getProperty(LayeredOptions.LAYERING_LAYER_CONSTRAINT) == LayerConstraint.LAST_SEPARATE)
                        || n2.hasProperty(LayeredOptions.LAYERING_LAYER_CONSTRAINT) && (n2.getProperty(LayeredOptions.LAYERING_LAYER_CONSTRAINT) == LayerConstraint.FIRST_SEPARATE
                                || n2.getProperty(LayeredOptions.LAYERING_LAYER_CONSTRAINT) == LayerConstraint.LAST_SEPARATE)) {
                        return 0;
                    }
                    LGraph lgraph = n1.getGraph();
                    // First check whether the transitive dependencies already determine the ordering.
                    int transitiveComparison = compareBasedOnTansitiveDependencies(n1, n2);
                    if (transitiveComparison != 0) {
                        return transitiveComparison;
                    }
                    // If this is not the case, consider the different ordering modes for group model order for the comparator.
                    if (n1.hasProperty(InternalProperties.MODEL_ORDER)
                            && n2.hasProperty(InternalProperties.MODEL_ORDER)) {
                        // The value is either a comparison both model orders or both group model orders with
                        // model order as secondary criterion.
                        int value = Integer.compare(
                                CMGroupModelOrderCalculator.calculateModelOrderOrGroupModelOrder(n1, n2, lgraph, lgraph.getProperty(InternalProperties.MAX_MODEL_ORDER_NODES)),
                                CMGroupModelOrderCalculator.calculateModelOrderOrGroupModelOrder(n2, n1, lgraph, lgraph.getProperty(InternalProperties.MAX_MODEL_ORDER_NODES)));
                        if (lgraph.getProperty(LayeredOptions.CONSIDER_MODEL_ORDER_GROUP_MODEL_ORDER_CM_GROUP_ORDER_STRATEGY) == GroupOrderStrategy.ONLY_WITHIN_GROUP) {
                            // Check list of enforced 
                            // If two nodes are in separate model order groups, do not enforce their ordering and just use barycenter.
                            if (n1.getProperty(LayeredOptions.CONSIDER_MODEL_ORDER_GROUP_MODEL_ORDER_CROSSING_MINIMIZATION_ID)
                                != n2.getProperty(LayeredOptions.CONSIDER_MODEL_ORDER_GROUP_MODEL_ORDER_CROSSING_MINIMIZATION_ID)) {
                                value = 0;
                            }
                        }
                        if (value < 0) {
                            updateBiggerAndSmallerAssociations(n1, n2);
                            return value;
                        } else if (value > 0) {
                            updateBiggerAndSmallerAssociations(n2, n1);
                            return value;
                        }
                    }
                    // If one of the nodes has no model order, fall back to the barycenter method.
                    return compareBasedOnBarycenter(n1, n2);
                };
    }

    /**
     * Don't use!
     * Only public to be accessible by a test.
     */
    @Override
    public void minimizeCrossings(final List<LNode> layer, final boolean preOrdered,
            final boolean randomize, final boolean forward) {

        if (randomize) {
            // Randomize barycenters (we don't need to update the edge count in this case;
            // there are no edges of interest since we're only concerned with one layer)
            randomizeBarycenters(layer);
        } else {
            // Calculate barycenters and assign barycenters to barycenterless node groups
            calculateBarycenters(layer, forward);
            fillInUnknownBarycenters(layer, preOrdered);
        }

        if (layer.size() > 1) {
            // Sort the vertices according to their barycenters
            if (layer.get(0).getGraph().getProperty(LayeredOptions.CROSSING_MINIMIZATION_FORCE_NODE_MODEL_ORDER)) {
                ModelOrderBarycenterHeuristic.insertionSort(layer, barycenterStateComparator,
                        (ModelOrderBarycenterHeuristic) this);
            } else {
                Collections.sort(layer, barycenterStateComparator);
            }

            // Resolve ordering constraints
            if (!layer.get(0).getGraph().getProperty(LayeredOptions.CROSSING_MINIMIZATION_FORCE_NODE_MODEL_ORDER)) {
                constraintResolver.processConstraints(layer);
            }
        }
    }
    
    /**
     * Compare two nodes based of the already collected transitive dependencies.
     * 
     * @param n1 The first node
     * @param n2 The second node
     * @return Returns -1 if n1 < n2, 0 if no ordering could be found, and 1 if n1 > n2
     */
    private int compareBasedOnTansitiveDependencies(final LNode n1, final LNode n2) {
        if (!biggerThan.containsKey(n1)) {
            biggerThan.put(n1, new HashSet<>());
        } else if (biggerThan.get(n1).contains(n2)) {
            return 1;
        }
        if (!biggerThan.containsKey(n2)) {
            biggerThan.put(n2, new HashSet<>());
        } else if (biggerThan.get(n2).contains(n1)) {
            return -1;
        }
        if (!smallerThan.containsKey(n1)) {
            smallerThan.put(n1, new HashSet<>());
        } else if (smallerThan.get(n1).contains(n2)) {
            return -1;
        }
        if (!smallerThan.containsKey(n2)) {
            smallerThan.put(n2, new HashSet<>());
        } else if (smallerThan.get(n2).contains(n1)) {
            return 1;
        }
        return 0;
    }
    
    /**
     * Order nodes by their barycenter and update their transitive dependencies.
     * @param n1 The first node
     * @param n2 The second node
     * @return -1, 0, 1, see {@link Comparator} for details.
     */
    private int compareBasedOnBarycenter(final LNode n1, final LNode n2) {
        BarycenterState s1 = stateOf(n1);
        BarycenterState s2 = stateOf(n2);
        if (s1.barycenter != null && s2.barycenter != null) {
            int value = s1.barycenter.compareTo(s2.barycenter);
            if (value < 0) {
                updateBiggerAndSmallerAssociations(n1, n2);
            } else if (value > 0) {
                updateBiggerAndSmallerAssociations(n2, n1);
            }
            return value;
        } else if (s1.barycenter != null) {
            updateBiggerAndSmallerAssociations(n1, n2);
            return -1;
        } else if (s2.barycenter != null) {
            updateBiggerAndSmallerAssociations(n2, n1);
            return 1;
        }
        return 0;
    }
    
    /**
     * Taken from {@code ModelOrderNodeComparator}. Safes transitive ordering dependencies of nodes for later use.
     * @param bigger
     * @param smaller
     */
    private void updateBiggerAndSmallerAssociations(final LNode bigger, final LNode smaller) {
        HashSet<LNode> biggerNodeBiggerThan = biggerThan.get(bigger);
        HashSet<LNode> smallerNodeBiggerThan = biggerThan.get(smaller);
        HashSet<LNode> biggerNodeSmallerThan = smallerThan.get(bigger);
        HashSet<LNode> smallerNodeSmallerThan = smallerThan.get(smaller);
        biggerNodeBiggerThan.add(smaller);
        smallerNodeSmallerThan.add(bigger);
        for (LNode verySmall : smallerNodeBiggerThan) {
            biggerNodeBiggerThan.add(verySmall);
            smallerThan.get(verySmall).add(bigger);
            smallerThan.get(verySmall).addAll(biggerNodeSmallerThan);
        }
        

        for (LNode veryBig : biggerNodeSmallerThan) {
            smallerNodeSmallerThan.add(veryBig);
            biggerThan.get(veryBig).add(smaller);
            biggerThan.get(veryBig).addAll(smallerNodeBiggerThan);
        }
    }
    
    /**
     * Sorts the layer by the given comparator with insertion sort to ensure that transitive ordering constraints are
     * respected. The usual quick sorting algorithm used by java does not do well here since the partial ordering
     * defined by the model order and the ordering of the barycenter values is prone to create a conflicting ordering.
     * Insertion sort can deal with that since no real node is moved past another real node if the model order forbids
     * it. However, the placement of dummy nodes and, therefore, long edges might suffer a little.
     * 
     * @param layer The layer to order. It is assumed that all "real" nodes which are here all nodes that have a model
     *      order are already ordered as intended.
     * @param comparator A comparator that uses a full ordering (barycenter values) and a partial ordering, which has
     *      a higher priority, (model order of real nodes) at the same time.
     * @param barycenterHeuristic The instance of this class to reset the transitive ordering associations.
     */
    public static void insertionSort(final List<LNode> layer, final Comparator<LNode> comparator,
            final ModelOrderBarycenterHeuristic barycenterHeuristic) {
        LNode temp;
        for (int i = 1; i < layer.size(); i++) {
            temp = layer.get(i);
            int j = i;
            while (j > 0 && comparator.compare(layer.get(j - 1), temp) > 0) {
                layer.set(j, layer.get(j - 1));
                j--;
            }
            layer.set(j, temp);
        }
        barycenterHeuristic.clearTransitiveOrdering();
    }
    
    /**
     * Clears the transitive ordering.
     */
    public void clearTransitiveOrdering() {
        this.biggerThan = new HashMap<>();
        this.smallerThan = new HashMap<>();
    }

}
