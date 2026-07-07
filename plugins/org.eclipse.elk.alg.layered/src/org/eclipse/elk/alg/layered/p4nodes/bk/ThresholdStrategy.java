/*******************************************************************************
 * Copyright (c) 2015 Kiel University and others.
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
package org.eclipse.elk.alg.layered.p4nodes.bk;

import java.util.Queue;
import java.util.Set;
import java.util.Stack;

import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.p4nodes.bk.BKAlignedLayout.HDirection;
import org.eclipse.elk.alg.layered.p4nodes.bk.BKAlignedLayout.VDirection;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.google.common.math.DoubleMath;

/**
 * 
 */
public abstract class ThresholdStrategy {

    // TODO make this an option?!
    private static final double THRESHOLD = Double.MAX_VALUE;
    
    private static final double EPSILON = 0.0001d; 
    
    // SUPPRESS CHECKSTYLE NEXT 24 VisibilityModifier
    /**
     * The currently processed layout with its iteration directions.
     */
    protected BKAlignedLayout bal;
    
    /**
     * The precalculated neighborhood information.
     */
    protected NeighborhoodInformation ni;
    
    /**
     * We keep track of which blocks have been completely finished.
     */
    protected Set<LNode> blockFinished = Sets.newHashSet();

    /**
     * A queue with blocks that are postponed during compaction.
     */
    protected Queue<Postprocessable> postProcessablesQueue = Lists.newLinkedList();
    
    /** 
     * A stack that is used to treat postponed nodes in reversed order.
     */
    protected Stack<Postprocessable> postProcessablesStack = new Stack<>();
    
    /**
     * Resets the internal state.
     * 
     * @param theBal
     *            The currently processed layout with its iteration directions.
     * @param theNi
     *            The precalculated neighborhood information of the graph.
     */
    public void init(final BKAlignedLayout theBal, final NeighborhoodInformation theNi) {
        this.bal = theBal;
        this.ni = theNi;
        blockFinished.clear();
        postProcessablesQueue.clear();
        postProcessablesStack.clear();
    }
    
    /**
     * Marks the block of which {@code n} is the root to be completely placed.
     * 
     * @param n
     *            the root of a block.
     */
    public void finishBlock(final LNode n) {
        blockFinished.add(n);
    }
    
    
    /* ------------------------------------------------
     *  Methods to be implemented by deriving classes.
     * ------------------------------------------------
     */
    
    /**
     * @param oldThresh
     *            an old, previously calculated threshold value
     * @param blockRoot
     *            the root node of the current block being placed
     * @param currentNode
     *            the currently processed node of a block. This can be equal to
     *            {@code blockRoot}.
     * @return a threshold value representing a bound that would allow an additional edge to be
     *         drawn straight.
     */
    public abstract double calculateThreshold(double oldThresh, LNode blockRoot, LNode currentNode);
    
    
    /**
     * Handle nodes that have been marked as having potential to 
     * lead to further straight edges after all blocks were initially placed.
     */
    public abstract void postProcess();
    
    // SUPPRESS CHECKSTYLE NEXT 20 VisibilityModifier
    /**
     * @param edge
     *            the edge for which the node is requested.
     * @param n
     *            a node edge is connected to.
     * @return for an edge {@code (o,n)}, return {@code o}.
     */
    protected LNode getOther(final LEdge edge, final LNode n) {
        if (edge.getSource().getNode() == n) {
            return edge.getTarget().getNode();
        } else if (edge.getTarget().getNode() == n) {
            return edge.getSource().getNode();
        } else {
            throw new IllegalArgumentException("Node " + n
                    + " is neither source nor target of edge " + edge);
        }
    }

    
    /**
     * {@link ThresholdStrategy} for the classic compaction phase of the original bk algorithm.
     * 
     * It calculates a threshold value such that it has no effect.
     */
    public static class NullThresholdStrategy extends ThresholdStrategy {
        @Override
        public double calculateThreshold(final double oldThresh, final LNode blockRoot,
                final LNode currentNode) {
            if (bal.vdir == VDirection.UP) {
                // new value calculated using min(a,thresh) --> thresh = +infty has no effect
                return Double.POSITIVE_INFINITY;
            } else {
                return Double.NEGATIVE_INFINITY;
            }
        }

        @Override
        public void postProcess() {
        }
    }

    /**
     * <ul>
     *  <li> Only calculates threshold for the first and last node of a block.</li>
     *  <li> Picks the first edge it encounters that is valid.</li>
     * </ul>
     */
    public static class SimpleThresholdStrategy extends ThresholdStrategy {
        
        @Override
        public double calculateThreshold(final double oldThresh,
                final LNode blockRoot, final LNode currentNode) {
            
            // just the root or last node of a block
            
            // Remember that for blocks with a single node both flags can be true
            boolean isRoot = blockRoot.equals(currentNode);
            boolean isLast = bal.align[currentNode.id].equals(blockRoot);
            
            if (!(isRoot || isLast)) {
                return oldThresh;
            }
            
            // Remember two things:
            //  1) it is not guaranteed that adjacent nodes are already placed
            //  2) blocks can consist of a single node implying that the current
            //     node is both the root and the last node
    
            double t = oldThresh;
            if (bal.hdir == HDirection.RIGHT) {
                
                if (isRoot) {
                    t = getBound(blockRoot, true);
                }
                if (Double.isInfinite(t) && isLast) {
                    t = getBound(currentNode, false);
                }
                
            } else { // LEFT
                
                if (isRoot) {
                    t = getBound(blockRoot, true);
                } 
                if (Double.isInfinite(t) && isLast) {
                    t = getBound(currentNode, false);
                }
            }
            
            return t;
        }
        
        
        /**
         * Only regards for root and last nodes of a block.
         * 
         * @param bal
         * @param currentNode
         *            a node of the block
         * @param isRoot
         *            whether {@code currentNode} is considered to be the root node of the current
         *            block. For a block that consists of a single node it is important to be able
         *            to regard it as root as well as as last node of a block.
         * 
         * @return a pair with an {@link LEdge} and a {@link Boolean}. If no valid edge was picked,
         *         the pair's first element is {@code null} and the second element indicates if
         *         there are possible candidate edges that might become valid at a later stage.
         */
        private Postprocessable pickEdge(final Postprocessable pp) {
        
            Iterable<LEdge> edges;
            if (pp.isRoot) {
                edges = bal.hdir == HDirection.RIGHT
                            ? pp.free.getIncomingEdges() : pp.free.getOutgoingEdges();
            } else {
                edges = bal.hdir == HDirection.LEFT 
                            ? pp.free.getIncomingEdges() : pp.free.getOutgoingEdges();
            }
            
            boolean hasEdges = false;
            for (LEdge e : edges) {
                
                // ignore in-layer edges unless the block is solely connected by in-layer edges
                //  rationale: With self-loops and feedback edges it can happen that blocks contain only dummy nodes, 
                //  thus are not connected to other blocks by non-inlayer edges. To avoid unnecessarily long edges 
                //  in these cases, such blocks are allowed to be handled here as well (to shorten the edges)
                boolean onlyDummies = bal.od[bal.root[pp.free.id].id];
                if (!onlyDummies && e.isInLayerEdge()) {
                    continue;
                }
                
                // in order to straighten 'e' the block represented by 'pp.free'
                // would have to be moved. However, since that block is already 
                // part of a straightened edge, it cannot be moved again
                if (bal.su[bal.root[pp.free.id].id] || bal.su[bal.root[pp.free.id].id]) {
                    continue;
                }
                
                hasEdges = true;
                
                // if the other node does not have a position yet, ignore this edge
                if (blockFinished.contains(bal.root[getOther(e, pp.free).id])) {
                    pp.hasEdges = true;
                    pp.edge = e;
                    return pp;
                }
            }
            
            // no edge picked
            pp.hasEdges = hasEdges;
            pp.edge = null;
            return pp;
        }
        
    
        /**
         * 
         */
        private double getBound(final LNode blockNode,
                final boolean isRoot) {
    
            double invalid = bal.vdir == VDirection.UP 
                    ? Double.POSITIVE_INFINITY : Double.NEGATIVE_INFINITY;
    
            final Postprocessable pick = pickEdge(new Postprocessable(blockNode, isRoot));
            
            // if edges exist but we couldn't find a good one
            if (pick.edge == null && pick.hasEdges) {
                 postProcessablesQueue.add(pick);
                 return invalid;
            } else if (pick.edge != null) {
    
                double threshold;
                LPort left = pick.edge.getSource();
                LPort right = pick.edge.getTarget();
    
                if (isRoot) {
                    // We handle the root (first) node of a block here
                    LPort rootPort = bal.hdir == HDirection.RIGHT ? right : left;
                    LPort otherPort = bal.hdir == HDirection.RIGHT ? left : right;
        
                    LNode otherRoot = bal.root[otherPort.getNode().id];
                    threshold = bal.y[otherRoot.id] 
                                      + bal.innerShift[otherPort.getNode().id]
                                      + otherPort.getPosition().y 
                                      + otherPort.getAnchor().y
                                      // root node
                                      - bal.innerShift[rootPort.getNode().id] 
                                      - rootPort.getPosition().y
                                      - rootPort.getAnchor().y;
                } else {
                    
                    // ... and the last node of a block here 
                    LPort rootPort = bal.hdir == HDirection.LEFT ? right : left;
                    LPort otherPort = bal.hdir == HDirection.LEFT ? left : right;
    
                    threshold = bal.y[bal.root[otherPort.getNode().id].id]
                            + bal.innerShift[otherPort.getNode().id]
                            + otherPort.getPosition().y
                            + otherPort.getAnchor().y
                            // root node
                            - bal.innerShift[rootPort.getNode().id]
                            - rootPort.getPosition().y
                            - rootPort.getAnchor().y;
                }
                
                // we are not allowed to move this block anymore 
                // in order to straighten another edge
                bal.su[bal.root[left.getNode().id].id] = true;
                bal.su[bal.root[right.getNode().id].id] = true;
                
                return threshold;
            }
            return invalid;
        }
        
        @Override
        public void postProcess() {
            
            // try original iteration order
            while (!postProcessablesQueue.isEmpty()) {
                
                // first is the node, second whether it is regarded as root
                Postprocessable pp = postProcessablesQueue.poll();
                Postprocessable pick = pickEdge(pp);

                if (pick.edge == null) {
                    continue;
                }
                
                LEdge edge = pick.edge;
                
                // ignore in-layer edges
                boolean onlyDummies = bal.od[bal.root[pp.free.id].id];
                if (!onlyDummies && edge.isInLayerEdge()) {
                    continue;
                }

                // try to straighten the edge ...
                boolean moved = process(pp);
                // if it wasn't possible try again later in the opposite iteration direction
                if (!moved) {
                    postProcessablesStack.push(pp);
                }
                
            }

            // reversed iteration order
            while (!postProcessablesStack.isEmpty()) {
                process(postProcessablesStack.pop());
            }
        }
        
        private boolean process(final Postprocessable pp) {
            assert pp.edge != null;
            
            LEdge edge = pp.edge;
            LPort fix;
            if (edge.getSource().getNode() == pp.free) {
                fix = edge.getTarget();
            } else {
                fix = edge.getSource();
            }
            LPort block; 
            if (edge.getSource().getNode() == pp.free) {
                block = edge.getSource();
            } else {
                block = edge.getTarget();
            }

            // t has to be the root node of a different block
            double delta = bal.calculateDelta(fix, block);
            
            if (delta > 0 && delta < THRESHOLD) {
                // target y larger than source y --> shift upwards?
                double availableSpace = bal.checkSpaceAbove(block.getNode(), delta, ni);
                assert DoubleMath.fuzzyEquals(availableSpace, 0, EPSILON) || availableSpace >= 0;
                bal.shiftBlock(block.getNode(), -availableSpace);
                return availableSpace > 0;
            } else if (delta < 0 && -delta < THRESHOLD) {
                
                // direction is up, we possibly shifted some blocks too far upward 
                // for an edge to be straight, so check if we can shift down again
                double availableSpace = bal.checkSpaceBelow(block.getNode(), -delta, ni);
                assert DoubleMath.fuzzyEquals(availableSpace, 0, EPSILON) || availableSpace >= 0;
                bal.shiftBlock(block.getNode(), availableSpace);
                return availableSpace > 0;
            }
            
            return false;
        }
    }
    
    /** 
     * Represents a unit to be post-processed. 
     */
    private static class Postprocessable {
        // SUPPRESS CHECKSTYLE NEXT 8 VisibilityModifier
        /** the node whose block can potentially be moved. */
        LNode free;
        /** whether {@code free} is the root node of its block. */
        boolean isRoot;
        /** whether {@code free} has edges. */
        boolean hasEdges;
        /** the edge that was selected to be straightened. */
        LEdge edge;
        
        Postprocessable(final LNode free, final boolean isRoot) {
            this.free = free;
            this.isRoot = isRoot;
        }
    }
    
}
