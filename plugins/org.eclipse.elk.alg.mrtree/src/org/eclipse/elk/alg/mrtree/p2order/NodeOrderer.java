/*******************************************************************************
 * Copyright (c) 2013 - 2022 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.alg.mrtree.p2order;

import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;
import java.util.ListIterator;
import java.util.stream.Collectors;
import org.eclipse.elk.alg.mrtree.TreeLayoutPhases;
import org.eclipse.elk.alg.mrtree.graph.TEdge;
import org.eclipse.elk.alg.mrtree.graph.TGraph;
import org.eclipse.elk.alg.mrtree.graph.TNode;
import org.eclipse.elk.alg.mrtree.intermediate.IntermediateProcessorStrategy;
import org.eclipse.elk.alg.mrtree.options.InternalProperties;
import org.eclipse.elk.alg.mrtree.options.MrTreeOptions;
import org.eclipse.elk.alg.mrtree.options.OrderWeighting;
import org.eclipse.elk.core.alg.ILayoutPhase;
import org.eclipse.elk.core.alg.LayoutProcessorConfiguration;
import org.eclipse.elk.core.util.IElkProgressMonitor;
import org.eclipse.elk.graph.properties.IProperty;
import org.eclipse.elk.graph.properties.PropertyHolderComparator;

/**
 * This phase orders the nodes of each level by separating the nodes into leaves and inner nodes.
 * It then fills whitespaces in the levels with corresponding leaves.
 * 
 * @author sor
 * @author sgu
 * @author jnc
 * @author sdo
 */
public class NodeOrderer implements ILayoutPhase<TreeLayoutPhases, TGraph> {

    /** intermediate processing configuration. */
    private static final LayoutProcessorConfiguration<TreeLayoutPhases, TGraph> INTERMEDIATE_PROCESSING_CONFIG =
            LayoutProcessorConfiguration.<TreeLayoutPhases, TGraph>create()
                    .before(TreeLayoutPhases.P2_NODE_ORDERING)
                        .add(IntermediateProcessorStrategy.ROOT_PROC)
                        .add(IntermediateProcessorStrategy.FAN_PROC)
                        .add(IntermediateProcessorStrategy.LEVEL_PROC);
    
    private OrderWeighting weighting;
    private boolean debug = false;

    /**
     * {@inheritDoc}
     */
    @Override
    public LayoutProcessorConfiguration<TreeLayoutPhases, TGraph> getLayoutProcessorConfiguration(final TGraph graph) {
        return INTERMEDIATE_PROCESSING_CONFIG;
    }

    /**
     * {@inheritDoc}
     */
    public void process(final TGraph tGraph, final IElkProgressMonitor progressMonitor) {

        progressMonitor.begin("Processor arrange node", 1);
        
        debug = tGraph.getProperty(MrTreeOptions.DEBUG_MODE);

        // elkjs-exclude-start
        if (debug) {
            progressMonitor.log("NodeOrderer!");
        }
        // elkjs-exclude-end

        // find the root node and add it to the root list, assuming that: 1. There is a root, 2. There is only one root
        TNode root = tGraph.getNodes().stream().
                filter(x -> x.getProperty(InternalProperties.ROOT)).
                findFirst().get();
        
        // call the right orderLevel on the first level, which is the root
        weighting = tGraph.getProperty(MrTreeOptions.WEIGHTING);
        if (weighting.equals(OrderWeighting.FAN) || weighting.equals(OrderWeighting.DESCENDANTS)) {
            orderLevelFanDescendants(Arrays.asList(root), progressMonitor.subTask(1.0f));
        } else if (weighting.equals(OrderWeighting.CONSTRAINT)) {
            orderLevelConstraint(Arrays.asList(root), progressMonitor.subTask(1.0f));
        } 

        progressMonitor.done();

    }

    /**
     * Order each level by separating the nodes into leaves and inner nodes. And then fill gaps with
     * corresponding leaves.
     * This ordering will be applied when using orderWighting FAN / DESCENDANTS
     * 
     * @param currentLevel
     * @param progressMonitor
     */
    private void orderLevelFanDescendants(final List<TNode> currentLevel, final IElkProgressMonitor progressMonitor) {
        progressMonitor.begin("Processor arrange level", 1);
        
        IProperty<Integer> sortProperty = InternalProperties.FAN;
        if (weighting.equals(OrderWeighting.DESCENDANTS)) {
            sortProperty = InternalProperties.DESCENDANTS;
        }
        
        int pos = 0;

        // sort all nodes in this level by their fan out
        // so the leaves are at the end of the list
        Collections.sort(currentLevel, PropertyHolderComparator.with(sortProperty));

        // find the first occurrence of a leave in the list
        int firstOcc = currentLevel.size();
        ListIterator<TNode> it = currentLevel.listIterator(currentLevel.size());
        boolean notNull = true;
        while (notNull && it.hasPrevious()) {
            TNode tNode = (TNode) it.previous();
            if ((tNode.getProperty(sortProperty) == 0)) {
                firstOcc--;
            } else {
                notNull = false;
            }
        }

        // Separate the level into leaves and inner nodes
        List<TNode> tmp = currentLevel.subList(0, firstOcc);
        LinkedList<TNode> inners = new LinkedList<TNode>(tmp);
        tmp = currentLevel.subList(firstOcc, currentLevel.size());
        LinkedList<TNode> leaves = new LinkedList<TNode>(tmp);

        // check if their are inner nodes left
        if (inners.isEmpty()) {
            // leave the leaves in their order
            for (TNode tENode : leaves) {
                tENode.setProperty(InternalProperties.POSITION, pos++);
            }
        } else {
            // order each level of descendants of the inner nodes
            int size = inners.size();
            for (TNode tPNode : inners) {
                tPNode.setProperty(InternalProperties.POSITION, pos++);

                // set the position of the children and set them in order
                LinkedList<TNode> children = tPNode.getChildrenCopy();
                orderLevelFanDescendants(children, progressMonitor.subTask(1 / size));

                // order the children by their reverse position
                Collections.sort(children,
                        Collections.reverseOrder(PropertyHolderComparator.with(InternalProperties.POSITION)));

                // reset the list of children with the new order
                List<TEdge> sortedOutEdges = new LinkedList<TEdge>();

                for (TNode tNode : children) {
                    for (TEdge tEdge : tPNode.getOutgoingEdges()) {
                        if (tEdge.getTarget() == tNode) {
                            sortedOutEdges.add(tEdge);
                        }
                    }
                }
                tPNode.getOutgoingEdges().clear();
                tPNode.getOutgoingEdges().addAll(sortedOutEdges);

                // fill gaps with leafs
                it = leaves.listIterator(leaves.size());
                int fillGap = tPNode.getOutgoingEdges().size();
                notNull = true;
                while ((0 < fillGap) && notNull && it.hasPrevious()) {
                    TNode tNode = (TNode) it.previous();
                    if ((tNode.getProperty(sortProperty) == 0)) {
                        tNode.setProperty(InternalProperties.POSITION, pos++);
                        fillGap--;
                        it.remove();
                    } else {
                        notNull = false;
                    }
                }
            }
        }
        progressMonitor.done();
    }
    
    /**
     * Order each level by moving the nodes with a position constraint option to that position.
     * 
     * @param currentLevel
     * @param progressMonitor
     */
    private void orderLevelConstraint(final List<TNode> currentLevel, final IElkProgressMonitor progressMonitor) {
        progressMonitor.begin("Processor arrange level", 1);

        // elkjs-exclude-start
        if (debug) {
            progressMonitor.log("OrderLevelConstraint!");
        }
        // elkjs-exclude-end
        
        List<TNode> undefinedNodes = currentLevel.stream().
                filter(x -> x.getProperty(MrTreeOptions.POSITION_CONSTRAINT) < 0).
                collect(Collectors.toList());
        List<TNode> inBoundNodes = currentLevel.stream().
                filter(x -> x.getProperty(MrTreeOptions.POSITION_CONSTRAINT) < currentLevel.size()
                         && x.getProperty(MrTreeOptions.POSITION_CONSTRAINT) >= 0).
                collect(Collectors.toList());
        List<TNode> outOfBoundNodes = currentLevel.stream().
                filter(x -> x.getProperty(MrTreeOptions.POSITION_CONSTRAINT) >= currentLevel.size()).
                collect(Collectors.toList());
        
        TNode[] sortedNodes = new TNode[currentLevel.size()];
        // Priority 1: Set non duplicate constraints
        for (int i = 0; i < inBoundNodes.size(); i++) {
            TNode curNode = inBoundNodes.get(i);
            int targetPos = curNode.getProperty(MrTreeOptions.POSITION_CONSTRAINT);
            if (targetPos >= 0 && targetPos < inBoundNodes.size() && sortedNodes[targetPos] == null) {
                sortedNodes[targetPos] = curNode;
                inBoundNodes.remove(i);
                i--;
            }
        }
        // Priority 2: Set duplicate constraints
        for (int i = 0; i < inBoundNodes.size(); i++) {
            TNode curNode = inBoundNodes.get(i);
            int targetPos = curNode.getProperty(MrTreeOptions.POSITION_CONSTRAINT);
            for (int j = 0;; j++) {
                int newTargetPos = targetPos + j;
                if (newTargetPos < sortedNodes.length && newTargetPos >= 0 && sortedNodes[newTargetPos] == null) {
                    sortedNodes[newTargetPos] = curNode;
                    inBoundNodes.remove(i);
                    i--;
                    break;
                }
                newTargetPos = targetPos - j;
                if (newTargetPos < sortedNodes.length && newTargetPos >= 0 && sortedNodes[newTargetPos] == null) {
                    sortedNodes[newTargetPos] = curNode;
                    inBoundNodes.remove(i);
                    i--;
                    break;
                }
            }
        }
        // Priority 3: Set out of bounds constraints
        outOfBoundNodes.sort((x, y) -> 
            -Integer.compare(
                    x.getProperty(MrTreeOptions.POSITION_CONSTRAINT), 
                    y.getProperty(MrTreeOptions.POSITION_CONSTRAINT)));
        for (int i = sortedNodes.length - 1; i >= 0; i--) {
            if (sortedNodes[i] == null && !outOfBoundNodes.isEmpty()) {
                sortedNodes[i] = outOfBoundNodes.get(0);
                outOfBoundNodes.remove(0);
            }
        }
        // Priority 4: Set no constraint nodes
        for (int i = 0; i < sortedNodes.length; i++) {
            if (sortedNodes[i] == null && !undefinedNodes.isEmpty()) {
                sortedNodes[i] = undefinedNodes.get(0);
                undefinedNodes.remove(0);
            }
        }
        
        // Set final node positions
        for (int i = 0; i < sortedNodes.length; i++) {
            sortedNodes[i].setProperty(InternalProperties.POSITION, i);
        }
        
        // Recursive calls, apply node positions
        TNode[] inners = currentLevel.stream().
                filter(x -> x.getProperty(InternalProperties.FAN) != 0).toArray(TNode[]::new);
        for (TNode tPNode : inners) {
            // Recursive call
            LinkedList<TNode> children = tPNode.getChildrenCopy();
            orderLevelConstraint(children, progressMonitor.subTask(1 / inners.length));
            
            // Order the children by their position
            Collections.sort(children, PropertyHolderComparator.with(InternalProperties.POSITION));

            // Reset the list of children with the new order
            List<TEdge> sortedOutEdges = new LinkedList<TEdge>();
            for (TNode tNode : children) {
                for (TEdge tEdge : tPNode.getOutgoingEdges()) {
                    if (tEdge.getTarget() == tNode) {
                        sortedOutEdges.add(tEdge);
                    }
                }
            }
            tPNode.getOutgoingEdges().clear();
            tPNode.getOutgoingEdges().addAll(sortedOutEdges);
        }
        
        progressMonitor.done();
    }
}
