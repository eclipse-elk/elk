/*******************************************************************************
 * Copyright (c) 2013 - 2022 Kiel University and others.
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
package org.eclipse.elk.alg.mrtree.p2order;

import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import org.eclipse.elk.alg.mrtree.TreeLayoutPhases;
import org.eclipse.elk.alg.mrtree.TreeUtil;
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

/**
 * This phase orders the nodes of each level by separating the children of nodes into leaves and
 * inner nodes. It then fill whitespaces in the levels with corresponding leaves.
 * 
 * It starts two levels above the deepest level, because the deepest level contains only nodes and
 * therefore no reordering is necessary. And the level above the deepest level contains only
 * children of the level above, which are ordered by the their parents.
 * 
 * @author sor
 * @author sgu
 * @author sdo
 */
public class OrderBalance implements ILayoutPhase<TreeLayoutPhases, TGraph> {

    /** intermediate processing configuration. */
    private static final LayoutProcessorConfiguration<TreeLayoutPhases, TGraph> INTERMEDIATE_PROCESSING_CONFIG =
            LayoutProcessorConfiguration.<TreeLayoutPhases, TGraph>create()
                    .before(TreeLayoutPhases.P2_NODE_ORDERING)
                        .add(IntermediateProcessorStrategy.ROOT_PROC)
                        .add(IntermediateProcessorStrategy.FAN_PROC)
                        .add(IntermediateProcessorStrategy.NEIGHBORS_PROC);

    /**
     * Tells the node order which weighting it should use.
     */
    private IProperty<Integer> weighting;

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

        // elkjs-exclude-start
        if (tGraph.getProperty(MrTreeOptions.DEBUG_MODE)) {
            progressMonitor.log("OrderBalance!");
        }
        // elkjs-exclude-end

        /** get the weighting from the userinterface */
        if (!tGraph.getProperty(MrTreeOptions.WEIGHTING).equals(OrderWeighting.FAN)) {
            weighting = InternalProperties.DESCENDANTS;
        } else {
            weighting = InternalProperties.FAN;
        }

        /** find the root of the component expected only one root exists */
        TNode root = null;
        Iterator<TNode> it = tGraph.getNodes().iterator();
        while (root == null && it.hasNext()) {
            TNode tNode = it.next();
            if (tNode.getProperty(InternalProperties.ROOT)) {
                root = tNode;
            }
        }

        /** check if a root is there */
        if (root != null) {

            /** start two levels above the deepest level at the leftmost node */
            TNode lM = TreeUtil.getLeftMost(root.getChildren());

            /** if there are only the root and one level or less no reordering is necessary */
            if (lM != null && lM.getParent() != root) {
                TNode parent = lM.getParent().getParent();
                TNode leftMost = parent;
                /** go to the leftmost node in this level */
                while (leftMost.getProperty(InternalProperties.LEFTNEIGHBOR) != null) {
                    leftMost = leftMost.getProperty(InternalProperties.LEFTNEIGHBOR);
                }
                /** start the order at the leftmost node */
                orderLevel(leftMost, false);

                /**
                 * reset the structure properties of each node to null, because the order of the
                 * graph has changed
                 */
                for (TNode tNode : tGraph.getNodes()) {
                    tNode.setProperty(InternalProperties.RIGHTNEIGHBOR, null);
                    tNode.setProperty(InternalProperties.LEFTNEIGHBOR, null);
                    tNode.setProperty(InternalProperties.RIGHTSIBLING, null);
                    tNode.setProperty(InternalProperties.LEFTSIBLING, null);
                }
            }
        }

    }

    /**
     * Order each level by seperating the children of the nodes into leaves and inner nodes. And
     * then fill gaps with corresponding leaves.
     * 
     * @param leftMost
     *            the leftmost node in a level
     */
    private void orderLevel(final TNode leftMost, final boolean odd) {
        if (leftMost != null) {

            /** copy current to iterate over the copy */
            TNode currentNode = leftMost;

            while (currentNode != null) {
                /**
                 * sort all children of this node by their fan out so the leaves are at the end of
                 * the list
                 */
                List<TEdge> outgoing = currentNode.getOutgoingEdges();

                Collections.sort(outgoing, new SortTEdgeTargetProperty(weighting));

                /**
                 * Add each child to a balanced list where the fat child are in the middle and the
                 * thin child are at the borders. Leaves fill the places between the inner child,
                 * also starting at fattest node in the middle.
                 * 
                 * eg. bigger number means fatter nodes, zero means leaf
                 * 
                 * unbalanced: 0 1 0 0 4 0 0 9 0 7 2 0 3 0 0 2 3 6 0 0 
                 * balanced  : 2 3 4 6 0 0 0 0 0 9 0 0 0 0 0 7 0 3 2 1
                 */
                List<TEdge> balanced = new LinkedList<TEdge>();

                boolean innerOdd = odd;
                while (!outgoing.isEmpty()) {
                    int gaps = outgoing.get(0).getTarget().getProperty(InternalProperties.FAN);
                    int index;

                    if (innerOdd) {
                        index = balanced.size();
                        balanced.add(outgoing.get(0));
                    } else {
                        index = 0;
                        balanced.add(index, outgoing.get(0));
                    }
                    outgoing.remove(0);
                    innerOdd = !innerOdd;

                    int indexEnd = outgoing.size();
                    boolean leavesOdd = odd;
                    while (0 < gaps && 0 < indexEnd) {
                        indexEnd--;
                        if (outgoing.get(indexEnd).getTarget().isLeaf()) {
                            gaps--;
                            if (leavesOdd) {
                                balanced.add(outgoing.get(indexEnd));
                            } else {
                                balanced.add(index, outgoing.get(indexEnd));
                            }
                            outgoing.remove(indexEnd);
                            leavesOdd = !leavesOdd;
                        } else {
                            gaps = 0;
                        }
                    }
                }

                /** reset the list of children with the new order */
                currentNode.getOutgoingEdges().addAll(balanced);

                /** go on with the next node to the right */
                currentNode = currentNode.getProperty(InternalProperties.RIGHTNEIGHBOR);
            }
            /** this level has been ordered, go on with the next level above */
            orderLevel(leftMost.getParent(), !odd);
        }
    }
    
    /**
     * A comparator for edge targets that uses the given property.
     */
    public static class SortTEdgeTargetProperty implements Comparator<TEdge> {
        private IProperty<Integer> property;

        SortTEdgeTargetProperty(final IProperty<Integer> property) {
            this.property = property;
        }
        
        /**
         * The compare method.
         * @param t1 The first edge
         * @param t2 The second edge
         * @return The difference for the property
         */
        public int compare(final TEdge t1, final TEdge t2) {
            return t2.getTarget().getProperty(property) - t1.getTarget().getProperty(property);
        }
    }
    
}
