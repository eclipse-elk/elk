/*******************************************************************************
 * Copyright (c) 2025 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 * 
 * SPDX-License-Identifier: EPL-2.0 
 *******************************************************************************/
package org.eclipse.elk.alg.layered.p4nodes;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import org.eclipse.elk.alg.layered.LayeredPhases;
import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.core.alg.ILayoutPhase;
import org.eclipse.elk.core.alg.LayoutProcessorConfiguration;
import org.eclipse.elk.core.util.IElkProgressMonitor;

import com.google.common.collect.Iterables;

/**
 * This algorithm was described by John Q. Walker II in John Q.Walker II, A Node-Positioning Algorithm for General
 * Trees, Software: Practice and Experience 20(7), pp. 685-705, July 1990. The implementation is largely based on the
 * NodePlacer class from the MrTree strategy in ELK.
 *
 */
public class WalkerNodePlacer implements ILayoutPhase<LayeredPhases, LGraph> {
    /*
     * In the layered strategy, the graph is always laid out horizontally (left-right). Therefore, we use the variables
     * prelim (preliminary indentation), modifier (indentation of a subtree beginning with that node), and depth
     * (basically height*spacing)
     */

    double spacingNodeNode;

    /*
     * (non-Javadoc)
     * 
     * @see org.eclipse.elk.core.alg.ILayoutProcessor#process(java.lang.Object,
     * org.eclipse.elk.core.util.IElkProgressMonitor)
     */
    @Override
    public void process(LGraph graph, IElkProgressMonitor monitor) {
        // TODO: deal with long edges: the dummy node is placed as if it were a real node, resulting in an angled edge.
        // This would be correct if it were a normal node, but looks unaesthetic for long edges. This can probably be
        // solved using post-processing, as it would be quite complicated to catch within the algorithm. Suggest looking
        // at how Brandes-Koepf and other placers deal with this problem, if they use post-processing, try to reuse that
        // post-processor.

        // remove cycles by removing certain edges (in effect constructing a spanning tree)
        List<LEdge> removedEdges = removeCycles(graph);

        monitor.begin("Walker node placement", 1);
        spacingNodeNode = graph.getProperty(LayeredOptions.SPACING_NODE_NODE);

        // TODO find a better way to get the root node of the graph
        // Maybe: create a dummy root node which is connected to all other root nodes such that there is only one root
        // node, compare RootProcessor in MrTree
        
        LNode root = graph.getLayers().get(0).getNodes().get(0);
        //root = getSuperRoot(graph);
        // here, only one root node is assumed to be the origin of the whole graph.

        // Do the preliminary positioning with a postorder walk.
        firstWalk(root, 0);
        // Apply the final positioning with a preorder walk
        secondWalk(root, 0);
        monitor.done();
        // restore removed edges
        restoreCycles(graph, removedEdges);
    }

    /**
     * @param graph
     * @return
     */
    private LNode getSuperRoot(LGraph graph) {
        List<LNode> rootCandidates = new ArrayList<>();
        for (Layer l : graph.getLayers()) {
            for (LNode n : l.getNodes()) {
                if (Iterables.size(n.getIncomingEdges()) == 0) {
                    rootCandidates.add(n);
                }
            }
        }
        if (rootCandidates.size() == 1) {
            return rootCandidates.get(0);
        }
        if (rootCandidates.size() == 0) {
            //then the graph is empty, or not a tree
            // TODO return dummy node??
        }
        // there is more than one root node
        // TODO create root node with all nodes in rootCandidates as children
        return null;
    }

    /**
     * @param graph
     * @param removedEdges
     */
    private void restoreCycles(LGraph graph, List<LEdge> removedEdges) {
        for (LEdge e : removedEdges) {
            e.getSource().getOutgoingEdges().add(e);
            e.getTarget().getIncomingEdges().add(e);
        }

    }

    /**
     * @param graph
     */
    private List<LEdge> removeCycles(LGraph graph) {
        List<LEdge> removedEdges = new ArrayList<>();
        Set<Set<LNode>> connectedNodesSet = new HashSet<>();
        for (Layer l : graph.getLayers()) {
            List<LNode> nodes = l.getNodes();
            for (LNode n : nodes) {
                LEdge[] edges = Iterables.toArray(n.getOutgoingEdges(), LEdge.class);
                LEdge current = null;
                for (int i = 0; i < edges.length; i++) {
                    current = edges[i]; // for debugging
                    LNode other = current.getOther(n);
                    // Maybe: find a method that gives me the "real" target, so no
                    // longEdge nodes, but only normal node (relevant when edges skip a layer)
                    // or deal with long edges otherwise
                    // -> Since I have to place the longEdge dummy nodes as well, leave them in

                    Set<LNode> nSet = findNodeInSetOfSets(connectedNodesSet, n);
                    Set<LNode> otherSet = findNodeInSetOfSets(connectedNodesSet, other);
                    if (nSet == null) {
                        if (otherSet == null) {
                            // neither of these nodes are connected to anything
                            Set<LNode> newSet = new HashSet<>();
                            newSet.add(n);
                            newSet.add(other);
                            connectedNodesSet.add(newSet);
                        } else {
                            // only other is connected, then add n to other's set
                            otherSet.add(n);
                        }
                    } else {
                        // n is connected
                        if (otherSet == null) {
                            // n is connected, but other is not
                            nSet.add(other);
                        } else {
                            // both n and other are in some set
                            if (nSet == otherSet) {
                                // there is a path between n and other
                                current.getSource().getOutgoingEdges().remove(current);
                                // current.setSource(null);
                                current.getTarget().getIncomingEdges().remove(current);
                                removedEdges.add(current);
                            } else {
                                // both n and other are in some set, but they are not connected
                                // keep the edge and unify the sets
                                connectedNodesSet.remove(otherSet);
                                nSet.addAll(otherSet);
                            }
                        }
                    }

                }
            }

        }
        return removedEdges;

    }

    /**
     * @param connectedNodesSet
     * @param n
     */
    private Set<LNode> findNodeInSetOfSets(Set<Set<LNode>> connectedNodesSet, LNode n) {
        Iterator<Set<LNode>> iter = connectedNodesSet.iterator();
        while (iter.hasNext()) {
            Set<LNode> next = iter.next();
            if (next.contains(n)) {
                return next;
            }
        }
        return null;
    }

    /*
     * (non-Javadoc)
     * 
     * @see org.eclipse.elk.core.alg.ILayoutPhase#getLayoutProcessorConfiguration(java.lang.Object)
     */
    @Override
    public LayoutProcessorConfiguration<LayeredPhases, LGraph> getLayoutProcessorConfiguration(LGraph graph) {
        return null;
    }

    /**
     * In this first postorder walk, every node of the tree is assigned a preliminary indentation. In addition, internal
     * nodes are given modifiers, which will be used to move their offspring to the right.
     * 
     * @param cN
     *            the root level of the tree
     * @param level
     *            the index of the passed level
     */
    private void firstWalk(final LNode cN, final int level) {
        cN.setProperty(InternalProperties.WALKER_MODIFIER, 0d);
        // get left sibling of cN
        LNode lS = leftSibling(cN);

        Iterable<LEdge> outEdges = cN.getOutgoingEdges();
        if (Iterables.size(outEdges) == 0) {
            // cN is a leaf
            if (lS == null) {
                // no left sibling -> set prelim to 0
                cN.setProperty(InternalProperties.WALKER_PRELIM, 0d);
            } else {// (lS != null)
                // cN has left sibling lS
                // calculate prelim of cN as follows:
                // cN.prelim = lS.prelim + nodeNodeSpacing + meanNodeWidth
                double prelimCN =
                        lS.getProperty(InternalProperties.WALKER_PRELIM) + spacingNodeNode + meanNodeWidth(cN, lS);
                cN.setProperty(InternalProperties.WALKER_PRELIM, prelimCN);
                // DEBUG here, only lS is respected, but we need to preserve spacingNodeNode even to non-sibling
                // neighbours

            }
        } else {
            // cN is not a leaf
            for (LEdge e : outEdges) {
                // recursion: call again on children of cN
                firstWalk(e.getTarget().getNode(), level + 1);
            }
            LNode leftMostChild = Iterables.getFirst(outEdges, null).getTarget().getNode();
            LNode rightMostChild = Iterables.getLast(outEdges, null).getTarget().getNode();
            double midPoint = (leftMostChild.getProperty(InternalProperties.WALKER_PRELIM)
                    + rightMostChild.getProperty(InternalProperties.WALKER_PRELIM)) / 2.0;
            if (lS == null) {
                // node has no left sibling
                cN.setProperty(InternalProperties.WALKER_PRELIM, midPoint);
            } else {
                // this node's children are shifted to the right as there is a left sibling
                double prelimCN =
                        lS.getProperty(InternalProperties.WALKER_PRELIM) + spacingNodeNode + meanNodeWidth(cN, lS);
                cN.setProperty(InternalProperties.WALKER_PRELIM, prelimCN);
                double modCN = prelimCN - midPoint;
                cN.setProperty(InternalProperties.WALKER_MODIFIER, modCN);
                apportion(cN, level);
            }
        }
        // TODO replace recursion with iterative construct
    }

    private void apportion(final LNode cN, final int level) {

        LNode leftMostChild = Iterables.getFirst(cN.getOutgoingEdges(), null).getTarget().getNode();
        // get the leftMostChild's left neighbour, if it exists
        LNode lMCNeighbor = null;
        if (leftMostChild != null) {
            int lMCIndex = leftMostChild.getIndex();
            if (lMCIndex > 0) {
                lMCNeighbor = leftMostChild.getLayer().getNodes().get(lMCIndex - 1);
            }
        }
        int compareDepth = 1;
        while (leftMostChild != null && lMCNeighbor != null) {
            double leftModSum = 0;
            double rightModSum = 0;
            LNode ancestorLeftmost = leftMostChild;
            LNode ancestorNeighbor = lMCNeighbor;
            // sum the modifiers of all ancestors according to the current level
            for (int i = 0; i < compareDepth; i++) {
                LEdge e = Iterables.getFirst(ancestorLeftmost.getIncomingEdges(), null);
                ancestorLeftmost = (e != null) ? e.getSource().getNode() : null;
                LEdge f = Iterables.getFirst(ancestorNeighbor.getIncomingEdges(), null);
                ancestorNeighbor = (f != null) ? f.getSource().getNode() : null;
                rightModSum +=
                        (ancestorLeftmost != null) ? ancestorLeftmost.getProperty(InternalProperties.WALKER_MODIFIER)
                                : 0;
                leftModSum +=
                        (ancestorNeighbor != null) ? ancestorNeighbor.getProperty(InternalProperties.WALKER_MODIFIER)
                                : 0;
                // added the (x != null) ? y : 0 in the two previous lines to prevent accessing property of null value

            }

            double prN = lMCNeighbor.getProperty(InternalProperties.WALKER_PRELIM);
            double prL = leftMostChild.getProperty(InternalProperties.WALKER_PRELIM);
            double mean = meanNodeWidth(leftMostChild, lMCNeighbor);
            double moveDistance = prN + leftModSum + spacingNodeNode + mean - prL - rightModSum;

            if (moveDistance > 0) {
                // Count interior sibling subtrees in LeftSiblings
                LNode leftSibling = cN;
                int numLeftSiblings = 0;
                while (leftSibling != null && leftSibling != ancestorNeighbor) {
                    numLeftSiblings++;
                    leftSibling = leftSibling(leftSibling);
                }
                // Apply portions to appropriate left sibling subtrees.
                if (leftSibling != null) {
                    double portion = moveDistance / (double) numLeftSiblings;
                    leftSibling = cN;
                    while (leftSibling != ancestorNeighbor) {
                        double newPr = leftSibling.getProperty(InternalProperties.WALKER_PRELIM) + moveDistance;
                        leftSibling.setProperty(InternalProperties.WALKER_PRELIM, newPr);
                        double newMod = leftSibling.getProperty(InternalProperties.WALKER_MODIFIER) + moveDistance;
                        leftSibling.setProperty(InternalProperties.WALKER_MODIFIER, newMod);
                        moveDistance -= portion;
                        leftSibling = leftSibling(leftSibling);
                    }
                } else {
                    // do nothing, some ancestor will do the apportioning
                    return;
                }
            }

            compareDepth++;
            if (Iterables.size(leftMostChild.getOutgoingEdges()) == 0) {
                // leftMostChild is leaf
                // set leftMostChild to the leftmost child of leftMostChild at level compareDepth
                LEdge e = Iterables.getFirst(leftMostChild.getOutgoingEdges(), null);
                if (e != null) {
                    leftMostChild = e.getTarget().getNode();
                } else {
                    leftMostChild = null;
                }

            } else {
                leftMostChild = Iterables.getFirst(leftMostChild.getOutgoingEdges(), null).getTarget().getNode();
            }
            if (leftMostChild != null) {
                int lMCIndex = leftMostChild.getIndex();
                if (lMCIndex > 0) {
                    lMCNeighbor = leftMostChild.getLayer().getNodes().get(lMCIndex - 1);
                }
            }
        }
    }

    /**
     * Calculate a node's left sibling. If the node has no left neighbour or that neighbour is not a sibling, return
     * null
     * 
     * @param cN
     *            the current Node
     * @return the left sibling or null
     */
    private LNode leftSibling(LNode cN) {
        LNode lS = null;
        if (cN.getIndex() > 0) {
            lS = cN.getLayer().getNodes().get(cN.getIndex() - 1);
            // if lS is only a neighbour, not a sibling, set it to null
            // assume that 1. siblings are ordered contiguously and
            // that 2. the graph is a tree, so no semicycles (nor cycles) exist
            if (Iterables.isEmpty(lS.getIncomingEdges())) {
                return null;
            }
            if (Iterables.isEmpty(cN.getIncomingEdges())) {// should not happen, but check for safety
                return null;
            }
            if (Iterables.getFirst(lS.getIncomingEdges(), null).getSource().getNode() != Iterables
                    .getFirst(cN.getIncomingEdges(), null).getSource().getNode()) {
                lS = null;
            }
        }
        return lS;
    }

    /**
     * This is a preorder walk. Set indentation by summing a nodes's preliminary indentation and the sum of its
     * ancestors' modifiers.
     * 
     * @param node
     *            the root of the currently examined subtree
     * @param modsum
     *            the sum of the node's ancestors' modifiers
     */
    private void secondWalk(final LNode node, final double modsum) {
        if (node != null) {
            double indent = node.getProperty(InternalProperties.WALKER_PRELIM) + modsum;
            // we do not deal with y coordinates here
            // in ELK Layered, this is done by a later stage
            // no property exists like LEVELHEIGHT from MrTree
            node.getPosition().y = indent - node.getSize().y / 2; // is this enough to respect node sizes?
            // no, if nodes are wider than the nodeNodeSpacing, they overlap
            if (!Iterables.isEmpty(node.getOutgoingEdges())) { // if node has no children
                secondWalk(Iterables.getFirst(node.getOutgoingEdges(), null).getTarget().getNode(),
                        modsum + node.getProperty(InternalProperties.WALKER_MODIFIER));
            }
            Layer layer = node.getLayer();
            if (node.getIndex() < layer.getNodes().size() - 1) {
                secondWalk(layer.getNodes().get(node.getIndex() + 1),
                        modsum + node.getProperty(InternalProperties.WALKER_MODIFIER));
            }
        }
    }

    /**
     * 
     * @param n1
     * @param n2
     * @return the mean width of nodes n1 and n2, including their relevant margins
     */
    private double meanNodeWidth(LNode n1, LNode n2) {
        double nodeWidth = 0;
        if (n1 != null) {
            nodeWidth += (n1.getSize().y + n1.getMargin().top + n1.getMargin().bottom) / 2.0;
        }
        if (n2 != null) {
            nodeWidth += (n2.getSize().y + n1.getMargin().top + n1.getMargin().bottom) / 2.0;
        }
        return nodeWidth;
    }

}
