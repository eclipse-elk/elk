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
package org.eclipse.elk.alg.layered.p3order;

import java.util.ArrayList;
import java.util.List;
import java.util.ListIterator;

import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LNode.NodeType;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.p3order.BarycenterHeuristic.BarycenterState;
import org.eclipse.elk.alg.layered.p3order.counting.IInitializable;
import org.eclipse.elk.core.util.Pair;

import com.google.common.collect.LinkedHashMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;

/**
 * Detects and resolves violated constraints. Inspired by
 * <ul>
 * <li>Michael Forster. A fast and simple heuristic for constrained two-level crossing reduction. In <i>Graph
 * Drawing</i>, volume 3383 of LNCS, pp. 206-216. Springer, 2005.</li>
 * </ul>
 * This constraint resolver relies on the assumption that all node groups have well-defined barycenter values and are
 * already sorted by these barycenter values.
 * 
 * <p>
 * Must be initialized using {@link IInitializable#init(java.util.List)}!
 * </p>
 */
public class ForsterConstraintResolver implements IInitializable {

    /** Whether there are successor constraints between non-dummies. */
    private boolean constraintsBetweenNonDummies = false;
    /** the layout units for handling dummy nodes for north / south ports. */
    private Multimap<LNode, LNode> layoutUnits;
    /** the barycenter values of every node in the graph, indexed by layer.id and node.id. */
    private BarycenterState[][] barycenterStates;
    /** the constraint groups, indexed by layer.id and node.id. */
    private ConstraintGroup[][] constraintGroups;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Initialisation

    /**
     * Creates a new instance which must be initialized using {@link IInitializable#init(java.util.List)}.
     * 
     * @param currentNodeOrder
     *            the current node order
     */
    public ForsterConstraintResolver(final LNode[][] currentNodeOrder) {
        if (currentNodeOrder.length > 0 && currentNodeOrder[0].length > 0) {
            constraintsBetweenNonDummies = currentNodeOrder[0][0].getGraph()
                    .getProperty(InternalProperties.IN_LAYER_SUCCESSOR_CONSTRAINTS_BETWEEN_NON_DUMMIES);
        }
        barycenterStates = new BarycenterState[currentNodeOrder.length][];
        constraintGroups = new ConstraintGroup[currentNodeOrder.length][];
        layoutUnits = LinkedHashMultimap.create();
    }

    @Override
    public void initAtLayerLevel(final int l, final LNode[][] nodeOrder) {
        barycenterStates[l] = new BarycenterState[nodeOrder[l].length];
        constraintGroups[l] = new ConstraintGroup[nodeOrder[l].length];
    }
    
    @Override
    public void initAtNodeLevel(final int l, final int n, final LNode[][] nodeOrder) {
        initAtNodeLevel(nodeOrder[l][n], true);
    }
    
    /**
     * Inits data structures for the given node.
     * 
     * @param node
     *            the node for which to initialize the data structures.
     * @param fullInit
     *            if {@code true}, this will also initialize barycenter states and layout units. This only needs to be
     *            done once, when the resolver is first initialized. The layout units remain stable, and the barycenter
     *            states are updated while the resolver is running.
     */
    private void initAtNodeLevel(final LNode node, final boolean fullInit) {
        int layerIndex = node.getLayer().id;
        int nodeIndex = node.id;
        
        constraintGroups[layerIndex][nodeIndex] = new ConstraintGroup(node);
        
        if (fullInit) {
            barycenterStates[layerIndex][nodeIndex] = new BarycenterState(node);
            
            LNode layoutUnit = node.getProperty(InternalProperties.IN_LAYER_LAYOUT_UNIT);
            if (layoutUnit != null) {
                layoutUnits.put(layoutUnit, node);
            }
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Accessors

    /**
     * Returns the barycenter states.
     */
    public BarycenterState[][] getBarycenterStates() {
        return barycenterStates;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Constraint Processing

    /**
     * Finds and handles violated in-layer successor constraints.
     * 
     * @param nodes
     *            list of nodes to handle constraints for.
     */
    public void processConstraints(final List<LNode> nodes) {
        // If there are successor constraints between regular (or normal) nodes, we have to apply a two-stage process:
        // in stage 1, we resolve constraints between the normal nodes. In stage 2, we ensure that all of their north /
        // south port dummies are in order and that no layout units overlap. For the latter to work properly, we need
        // the order of normal nodes to be fixed, hence the first stage.
        if (constraintsBetweenNonDummies) {
            processConstraints(nodes, true);
            nodes.stream().forEach(node -> initAtNodeLevel(node, false));
        }
        
        processConstraints(nodes, false);
    }
    
    private void processConstraints(final List<LNode> nodes, final boolean onlyBetweenNormalNodes) {
        List<ConstraintGroup> groups = new ArrayList<>(nodes.size());
        for (LNode node : nodes) {
            groups.add(constraintGroups[node.getLayer().id][node.id]);
        }

        // Build the constraints graph
        buildConstraintsGraph(groups, onlyBetweenNormalNodes);

        // Find violated vertices
        Pair<ConstraintGroup, ConstraintGroup> violatedConstraint = null;
        while ((violatedConstraint = findViolatedConstraint(groups)) != null) {
            handleViolatedConstraint(violatedConstraint.getFirst(), violatedConstraint.getSecond(), groups);
        }

        // Apply the determined order
        nodes.clear();
        for (ConstraintGroup group : groups) {
            for (LNode node : group.getNodes()) {
                nodes.add(node);
                stateOf(node).barycenter = group.getBarycenter();
            }
        }
    }

    /**
     * Build the constraint graph for the given vertices. The constraint graph is created from the predefined
     * <em>in-layer successor constraints</em> and the <em>layout units</em>.
     * 
     * @param groups
     *            the array of single-node vertices sorted by their barycenter values.
     */
    private void buildConstraintsGraph(final List<ConstraintGroup> groups, final boolean onlyBetweenNormalNodes) {
        // Reset the constraint fields
        for (ConstraintGroup group : groups) {
            group.resetOutgoingConstraints();
            group.incomingConstraintsCount = 0;
        }

        // Iterate through the vertices, adding the necessary constraints
        LNode lastNonDummyNode = null;
        for (ConstraintGroup group : groups) {
            // at this stage all groups should consist of a single node
            LNode node = group.getNode();
            
            // We may want to skip this
            if (onlyBetweenNormalNodes && node.getType() != NodeType.NORMAL) {
                continue;
            }

            // Add the constraints given by the vertex's node
            for (LNode successor : node.getProperty(InternalProperties.IN_LAYER_SUCCESSOR_CONSTRAINTS)) {
                if (!onlyBetweenNormalNodes || successor.getType() == NodeType.NORMAL) {
                    group.getOutgoingConstraints().add(groupOf(successor));
                    groupOf(successor).incomingConstraintsCount++;
                }
            }

            // Check if we're processing a a normal, none-dummy node. If this is the case, we may want to insert
            // additional constraints between the previous normal node's layout unit and the current node's layout
            // unit to be sure they won't mix. This requires that the normal nodes won't change order, however.
            if (!onlyBetweenNormalNodes && node.getType() == NodeType.NORMAL) {
                // If we already processed another normal, non-dummy node, we need to add
                // constraints from all of that other node's layout unit's vertices to this
                // node's layout unit's vertices
                if (lastNonDummyNode != null) {
                    for (LNode lastUnitNode : layoutUnits.get(lastNonDummyNode)) {
                        for (LNode currentUnitNode : layoutUnits.get(node)) {
                            groupOf(lastUnitNode).getOutgoingConstraints().add(groupOf(currentUnitNode));
                            groupOf(currentUnitNode).incomingConstraintsCount++;
                        }
                    }
                }

                lastNonDummyNode = node;
            }
        }
    }

    /**
     * Returns a violated constraint, if any is left. Constraint violation detection is based on the barycenter values
     * of the node groups, hence it is a critical requirement that the node groups are sorted by their barycenter
     * values.
     * 
     * @param groups
     *            list of vertices.
     * @return the two vertices whose constraint is violated, or {@code null} if none could be found. The two vertices
     *         are returned in the order they should appear in, not in the order that violates their constraint.
     */
    private Pair<ConstraintGroup, ConstraintGroup> findViolatedConstraint(final List<ConstraintGroup> groups) {
        List<ConstraintGroup> activeGroups = null;

        // Iterate through the constrained vertices
        double lastValue = Short.MIN_VALUE;
        for (ConstraintGroup group : groups) {
            assert group.getBarycenter() != null && group.getBarycenter() >= lastValue;
            lastValue = group.getBarycenter();
            group.resetIncomingConstraints();

            // Find sources of the constraint graph to start the constraints check
            if (group.hasOutgoingConstraints() && group.incomingConstraintsCount == 0) {
                if (activeGroups == null) {
                    activeGroups = Lists.newArrayList();
                }
                activeGroups.add(group);
            }
        }

        // Iterate through the active node groups to find one with violated constraints
        if (activeGroups != null) {
            while (!activeGroups.isEmpty()) {
                ConstraintGroup group = activeGroups.remove(0);

                // See if we can find a violated constraint
                if (group.hasIncomingConstraints()) {
                    for (ConstraintGroup predecessor : group.getIncomingConstraints()) {
                        if (predecessor.getBarycenter().floatValue() == group.getBarycenter().floatValue()) {
                            if (groups.indexOf(predecessor) > groups.indexOf(group)) {
                                // The predecessor has equal barycenter, but higher index
                                return Pair.of(predecessor, group);
                            }
                        } else if (predecessor.getBarycenter() > group.getBarycenter()) {
                            // The predecessor has greater barycenter and thus also higher index
                            return Pair.of(predecessor, group);
                        }
                    }
                }

                // No violated constraints; add outgoing constraints to the respective incoming list
                for (ConstraintGroup successor : group.getOutgoingConstraints()) {
                    List<ConstraintGroup> successorIncomingList = successor.getIncomingConstraints();
                    successorIncomingList.add(0, group);

                    if (successor.incomingConstraintsCount == successorIncomingList.size()) {
                        activeGroups.add(successor);
                    }
                }
            }
        }

        // No violated constraints found
        return null;
    }

    /** Delta that two barycenters can differ by to still be considered equal. */
    private static final float BARYCENTER_EQUALITY_DELTA = 0.0001F;

    /**
     * Handles the case of a violated constraint. The node groups must be sorted by their barycenter values. After this
     * method has finished, the list of node groups is smaller by one element, since two node groups have been unified,
     * but the list is still correctly sorted by barycenter values.
     * 
     * @param firstNodeGroup
     *            the node group with violated outgoing constraint
     * @param secondNodeGroup
     *            the node group with violated incoming constraint
     * @param nodeGroups
     *            the list of vertices
     */
    private void handleViolatedConstraint(final ConstraintGroup firstNodeGroup, final ConstraintGroup secondNodeGroup,
            final List<ConstraintGroup> nodeGroups) {

        // Create a new vertex from the two constrain-violating vertices; this also
        // automatically calculates the new vertex's barycenter value
        ConstraintGroup newNodeGroup = new ConstraintGroup(firstNodeGroup, secondNodeGroup);
        assert newNodeGroup.getBarycenter() + BARYCENTER_EQUALITY_DELTA >= secondNodeGroup.getBarycenter();
        assert newNodeGroup.getBarycenter() - BARYCENTER_EQUALITY_DELTA <= firstNodeGroup.getBarycenter();

        // Iterate through the vertices. Remove the old vertices. Insert the new one
        // according to the barycenter value, thereby keeping the list sorted. Along
        // the way, constraint relationships will be updated
        ListIterator<ConstraintGroup> nodeGroupIterator = nodeGroups.listIterator();
        boolean alreadyInserted = false;
        while (nodeGroupIterator.hasNext()) {
            ConstraintGroup nodeGroup = nodeGroupIterator.next();

            if (nodeGroup == firstNodeGroup || nodeGroup == secondNodeGroup) {
                // Remove the two node groups with violated constraint from the list
                nodeGroupIterator.remove();
            } else if (!alreadyInserted && nodeGroup.getBarycenter() > newNodeGroup.getBarycenter()) {
                // If we haven't inserted the new node group into the list already, do that now.
                // Note: we're not calling next() again. This means that during the next iteration,
                // we will again be looking at the current node group. But then, alreadyInserted will
                // be true and we can look at that node group's outgoing constraints.
                nodeGroupIterator.previous();
                nodeGroupIterator.add(newNodeGroup);

                alreadyInserted = true;
            } else if (nodeGroup.hasOutgoingConstraints()) {
                // Check if the vertex has any constraints with the former two vertices
                boolean firstNodeGroupConstraint = nodeGroup.getOutgoingConstraints().remove(firstNodeGroup);
                boolean secondNodeGroupConstraint = nodeGroup.getOutgoingConstraints().remove(secondNodeGroup);

                if (firstNodeGroupConstraint || secondNodeGroupConstraint) {
                    nodeGroup.getOutgoingConstraints().add(newNodeGroup);
                    newNodeGroup.incomingConstraintsCount++;
                }
            }
        }

        // If we haven't inserted the new node group already, add it to the end
        if (!alreadyInserted) {
            nodeGroups.add(newNodeGroup);
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Utilities

    private ConstraintGroup groupOf(final LNode node) {
        return constraintGroups[node.getLayer().id][node.id];
    }

    private BarycenterState stateOf(final LNode node) {
        return barycenterStates[node.getLayer().id][node.id];
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Constraint Groups

    /**
     * A constraint group contains one or more nodes. Constraint groups are used to model sets of nodes that are placed
     * next to each other. A constraint group contains methods to calculate its barycenter value, to merge with another
     * group and to generally do cool stuff.
     */
    public final class ConstraintGroup {

        /** The sum of the node weights. Each node weight is the sum of the weights of adjacent ports. */
        private double summedWeight;
        /** The number of ports relevant to the barycenter calculation. */
        private int degree;
        /** List of nodes this vertex consists of. */
        private final LNode[] nodes;
        /** List of outgoing constraints. */
        private List<ConstraintGroup> outgoingConstraints;
        /** List of incoming constraints. */
        private List<ConstraintGroup> incomingConstraints;
        /** The number of incoming constraints. */
        private int incomingConstraintsCount;
        
        ///////////////////////////////////////////////////////////////////////////////////////////////////
        // Construction

        /**
         * Constructs a new instance containing the given node.
         * 
         * @param node
         *            the node the group should contain
         */
        public ConstraintGroup(final LNode node) {
            nodes = new LNode[] { node };
        }

        /**
         * Constructs a new group that is the concatenation of the given two groups. The incoming constraints count
         * is set to zero, while the lists of successors are merged, updating the successors' incoming count
         * appropriately if both groups are predecessors. The new barycenter is derived from the barycenters of the
         * given node groups.
         * 
         * @param nodeGroup1
         *            the first vertex
         * @param nodeGroup2
         *            the second vertex
         */
        public ConstraintGroup(final ConstraintGroup nodeGroup1, final ConstraintGroup nodeGroup2) {
            // create a combined nodes array
            int length1 = nodeGroup1.nodes.length;
            int length2 = nodeGroup2.nodes.length;
            nodes = new LNode[length1 + length2];
            for (int i = 0; i < length1; i++) {
                nodes[i] = nodeGroup1.nodes[i];
            }
            for (int i = 0; i < length2; i++) {
                nodes[length1 + i] = nodeGroup2.nodes[i];
            }

            // Add constraints, taking care not to add any constraints to vertex1 or vertex2
            // and to decrement the incoming constraints count of those that are successors to both
            if (nodeGroup1.outgoingConstraints != null) {
                this.outgoingConstraints = Lists.newLinkedList(nodeGroup1.outgoingConstraints);
                this.outgoingConstraints.remove(nodeGroup2);
                if (nodeGroup2.outgoingConstraints != null) {
                    for (ConstraintGroup candidate : nodeGroup2.outgoingConstraints) {
                        if (candidate == nodeGroup1) {
                            continue;
                        } else if (outgoingConstraints.contains(candidate)) {
                            // The candidate was in both vertices' successor list
                            candidate.incomingConstraintsCount--;
                        } else {
                            outgoingConstraints.add(candidate);
                        }
                    }
                }
            } else if (nodeGroup2.outgoingConstraints != null) {
                this.outgoingConstraints = Lists.newLinkedList(nodeGroup2.outgoingConstraints);
                this.outgoingConstraints.remove(nodeGroup1);
            }

            summedWeight = nodeGroup1.summedWeight + nodeGroup2.summedWeight;
            degree = nodeGroup1.degree + nodeGroup2.degree;

            if (degree > 0) {
                setBarycenter(summedWeight / degree);
            } else if (nodeGroup1.getBarycenter() != null && nodeGroup2.getBarycenter() != null) {
                setBarycenter((nodeGroup1.getBarycenter() + nodeGroup2.getBarycenter()) / 2);
            } else if (nodeGroup1.getBarycenter() != null) {
                setBarycenter(nodeGroup1.getBarycenter());
            } else if (nodeGroup2.getBarycenter() != null) {
                setBarycenter(nodeGroup2.getBarycenter());
            }
        }
        
        ///////////////////////////////////////////////////////////////////////////////////////////////////
        // Accessors

        /**
         * @param barycenter
         *            the barycenter value
         */
        public void setBarycenter(final Double barycenter) {
            for (LNode node : nodes) {
                stateOf(node).barycenter = barycenter;
            }
        }

        /**
         * @return barycenter of current constraint group.
         */
        public Double getBarycenter() {
            return stateOf(nodes[0]).barycenter;
        }

        /**
         * Returns the list of outgoing constraints, creating it if not yet done before.
         * 
         * @return the outgoing constraints list of the node group
         */
        public List<ConstraintGroup> getOutgoingConstraints() {
            if (outgoingConstraints == null) {
                outgoingConstraints = Lists.newArrayList();
            }
            return outgoingConstraints;
        }

        /**
         * Reset the list of outgoing constraints to {@code null}.
         */
        public void resetOutgoingConstraints() {
            outgoingConstraints = null;
        }

        /**
         * Determine whether there are any outgoing constraints.
         * 
         * @return true if there are outgoing constraints
         */
        public boolean hasOutgoingConstraints() {
            return outgoingConstraints != null && outgoingConstraints.size() > 0;
        }

        /**
         * Returns the list of incoming constraints, creating it if not yet done before.
         * 
         * @return the incoming constraints list of the node group
         */
        public List<ConstraintGroup> getIncomingConstraints() {
            if (incomingConstraints == null) {
                incomingConstraints = Lists.newArrayList();
            }
            return incomingConstraints;
        }

        /**
         * Reset the list of incoming constraints to {@code null}.
         */
        public void resetIncomingConstraints() {
            incomingConstraints = null;
        }

        /**
         * Determine whether there are any incoming constraints.
         * 
         * @return true if there are incoming constraints
         */
        public boolean hasIncomingConstraints() {
            return incomingConstraints != null && incomingConstraints.size() > 0;
        }

        /**
         * Returns the array of nodes.
         * 
         * @return the contained nodes of the node group
         */
        public LNode[] getNodes() {
            return nodes;
        }

        /**
         * Returns the contained node. This may only be used for node groups with exactly one node.
         * 
         * @return the contained node
         */
        public LNode getNode() {
            assert nodes.length == 1;
            return nodes[0];
        }
        
        ///////////////////////////////////////////////////////////////////////////////////////////////////
        // Utilities

        @Override
        public String toString() {
            StringBuilder sb = new StringBuilder();
            sb.append('[');
            for (int i = 0; i < nodes.length; i++) {
                sb.append(nodes[i].toString());
                if (getBarycenter() != null) {
                    sb.append("<").append(getBarycenter().toString()).append(">");
                }
                if (i < nodes.length - 1) {
                    sb.append(", ");
                }
            }
            return sb.append(']').toString();
        }
    }
}
