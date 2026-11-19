/*******************************************************************************
 * Copyright (c) 2016 Kiel University and others.
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
package org.eclipse.elk.alg.common.networksimplex;

import java.util.ArrayDeque;
import java.util.Arrays;
import java.util.Deque;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;
import java.util.Set;

import org.eclipse.elk.core.util.BasicProgressMonitor;
import org.eclipse.elk.core.util.IElkProgressMonitor;
import org.eclipse.elk.core.util.Pair;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;

/**
 * The main class of the network simplex layerer component. It offers an algorithm to determine an
 * optimal layering of all nodes in the graph concerning a minimal length of all edges using the
 * network simplex algorithm described in
 * <ul>
 * <li>Emden R. Gansner, Eleftherios Koutsofios, Stephen C. North, Kiem-Phong Vo, A technique for
 * drawing directed graphs. <i>Software Engineering</i> 19(3), pp. 214-230, 1993.</li>
 * </ul>
 * 
 * <dl>
 * <dt>Precondition:</dt>
 * <dd>the graph has no cycles</dd>
 * <dt>Postcondition:</dt>
 * <dd>all nodes have been assigned a layer such that edges connect only nodes from layers with
 * increasing indices</dd>
 * </dl>
 */
public final class NetworkSimplex  {

    // configuration of this network simplex instance
    
    /** The number of nodes in each layer of a previous layering to be considered during {@link #normalize()} 
     * and {@link #balance(int[])}. That is, for an index {@code i}, {@code previousLayeringNodeCounts[i]} holds 
     * the number of nodes that are have previously been placed in layer {@code i}. */
    private int[] previousLayeringNodeCounts;
    /** Whether to apply {@link #balance(int[])}. */
    private boolean balance = false;
    /** A limit on the number of iterations. */
    private int iterationLimit = Integer.MAX_VALUE;
    /** Empirically determined threshold when removing subtrees pays off. */
    private static final int REMOVE_SUBTREES_THRESH = 40;
    
    /** Small value smaller than zero. Used to check whether cut values are small than zero and to deal with 
     *  imprecision of double computations. */
    private static final double FUZZY_ST_ZERO = -1e-10;
    
    /** Use {@link #forGraph(NGraph)}. */
    private NetworkSimplex() {
    }

    /**
     * @param graph
     *            the graph for which to execute the network simplex
     * @return a new instance of a {@link NetworkSimplex} algorithm.
     */
    public static NetworkSimplex forGraph(final NGraph graph) {
        NetworkSimplex ns = new NetworkSimplex();
        ns.graph = graph;
        return ns;
    }

    /**
     * It balances the layering concerning its width, i.e. the number of nodes in each layer. If the
     * graph allows multiple optimal layerings regarding a minimal edge length, this method moves
     * separate nodes to a layer with a minimal amount of currently contained nodes with respect to
     * the retention of feasibility and optimality of the given layering.
     * 
     * @param doBalance
     *            whether to apply a balancing
     * @return the {@link NetworkSimplex} instance for further configuration or execution.
     */
    public NetworkSimplex withBalancing(final boolean doBalance) {
        this.balance = doBalance;
        return this;
    }
    
    /**
     * Previously layered nodes may become relevant when moving nodes to layers with fewer nodes
     * during balancing.
     * 
     * @param considerPreviousLayering
     *            whether previously layered nodes should be considered during.
     * @return the {@link NetworkSimplex} instance for further configuration or execution.
     * @see #withBalancing(boolean)
     */
    public NetworkSimplex withPreviousLayering(final int[] considerPreviousLayering) {
        this.previousLayeringNodeCounts = considerPreviousLayering;
        return this;
    }

    /**
     * Since there is a theoretical possibility that the network simplex does not terminate Gansner
     * et al. propose to incorporate an iteration limit. However, in practice this shouldn't happen.
     * 
     * @param limit
     *            the maximum number of iterations of the network simplex algorithm.
     * @return the {@link NetworkSimplex} instance for further configuration or execution.
     */
    public NetworkSimplex withIterationLimit(final int limit) {
        this.iterationLimit = limit;
        return this;
    }
    
    // ================================== Attributes ==============================================

    /** The graph all methods in this class operate on. */
    private NGraph graph;

    /** An {@code ArrayList} containing all edges in the graph. */
    private List<NEdge> edges;
    
    /** A {@code HashSet} containing all edges that are part of the spanning tree. */
    private Set<NEdge> treeEdges;

    /**
     * A {@code LinkedList} containing all source nodes of the graph, i.e. all nodes that have no
     * incident incoming edges.
     */
    private List<NNode> sources;

    /**
     * A flag indicating whether a specified edge has been visited during DFS-traversal. This array
     * has to be filled with {@code false} each time, before a DFS-based method is invoked.
     */
    private boolean[] edgeVisited;

    /**
     * The current postorder traversal number used by {@code postorderTraversal()} to assign an
     * unique traversal ID to each node.
     * 
     * @see #postorderTraversal(NNode)
     */
    private int postOrder;

    /**
     * The postorder traversal ID of each node determined by {@code postorderTraversal()}.
     * 
     * @see #postorderTraversal(NNode)
     */
    private int[] poID;

    /**
     * The lowest postorder traversal ID of each nodes reachable through a node lower in the
     * traversal tree determined by {@code postorderTraversal}.
     * 
     * @see #postorderTraversal(NNode)
     */
    private int[] lowestPoID;

    /**
     * The cut value of every edge defined as follows: If the edge is deleted, the spanning tree
     * breaks into two connected components, the head component containing the target node of the
     * edge and the tail component containing the source node of the edge. The cut value is the sum
     * of the weight (here {@code 1}) of all edges going from the tail to the head component,
     * including the tree edge, minus the sum of the weights of all edges from the head to the tail
     * component.
     * 
     * @see #cutvalues()
     */
    private double[] cutvalue;
    
    /**
     * Nodes that are part of subtrees of the graph. They will be removed prior to the actual
     * execution of the network simplex since positioning them with minimal edge length is trivial.
     * 
     * @see #removeSubtrees()
     * @see #reattachSubtrees()
     */
    private Deque<Pair<NNode, NEdge>> subtreeNodesStack;

    // =============================== Initialization Methods =====================================

    /**
     * Helper method for the network simplex layerer. It instantiates all necessary attributes for
     * the execution of the network simplex layerer and initializes them with their default values.
     * All edges in the connected component given by the input argument will be determined, as well
     * as the number of incoming and outgoing edges of each node ( {@code inDegree}, respectively
     * {@code outDegree}). All sinks and source nodes in the connected component identified in this
     * step will be added to {@code sinks}, respectively {@code sources}.
     * 
     */
    private void initialize() {
        // initialize node attributes
        int numNodes = graph.nodes.size();
        for (NNode n : graph.nodes) {
            n.treeNode = false;
        }
        poID = new int[numNodes];
        lowestPoID = new int[numNodes];
        sources = Lists.newArrayList();

        // determine edges and re-index nodes
        int index = 0;
        List<NEdge> theEdges = Lists.newArrayList();
        for (NNode node : graph.nodes) {
            node.internalId = index++;
            // add node to sinks, resp. sources
            if (node.getIncomingEdges().size() == 0) {
                sources.add(node);
            }
            theEdges.addAll(node.getOutgoingEdges());
        }
        // re-index edges
        int counter = 0;
        for (NEdge edge : theEdges) {
            edge.internalId = counter++;
            edge.treeEdge = false;
        }
        // initialize edge attributes
        int numEdges = theEdges.size();
        if (cutvalue == null || cutvalue.length < numEdges) {
            cutvalue = new double[numEdges];
            edgeVisited = new boolean[numEdges];
        } else {
            Arrays.fill(edgeVisited, false);
        }
        edges = theEdges;
        // we iterate over this set, thus we have to use a linked hash set 
        // to get a deterministic iteration order
        treeEdges = Sets.newLinkedHashSetWithExpectedSize(edges.size());
        postOrder = 1;
    }

    /**
     * Release all created resources so the GC can reap them.
     */
    private void dispose() {
        this.cutvalue = null;
        this.edges = null;
        this.treeEdges = null;
        this.edgeVisited = null;
        this.lowestPoID = null;
        this.poID = null;
        this.sources = null;
        this.subtreeNodesStack = null;
    }

    // ============================== Network-Simplex Algorithm ===================================

    /**
     * Determine the optimal layering. 
     */
    public void execute() {
        execute(new BasicProgressMonitor());
    }
    
    /**
     * Determine the optimal layering.  
     * 
     * @param monitor
     *            a progress monitor
     */
    public void execute(final IElkProgressMonitor monitor) {
        monitor.begin("Network simplex", 1);

        if (graph.nodes.size() < 1) {
            monitor.done();
            return;
        }
        
        // reset any old layering
        for (NNode node : graph.nodes) {
            node.layer = 0;
        }
        
        // remove leafs
        boolean removeSubtrees = graph.nodes.size() >= REMOVE_SUBTREES_THRESH;
        if (removeSubtrees) {
            removeSubtrees();
        }

        // init all the data structures we use
        initialize();
        // determine an initial feasible layering
        feasibleTree();
        // improve the initial layering until it is optimal
        NEdge e = leaveEdge();
        int iter = 0;
        while (e != null && iter < iterationLimit) {
            // current layering is not optimal
            exchange(e, enterEdge(e));
            e = leaveEdge();
            iter++;
        }

        // re-attach leafs
        if (removeSubtrees) {
            reattachSubtrees();
        }
        
        // normalize and, if desired, balance
        //   both methods must work on the NNode#layer field
        if (balance) {
            balance(normalize());
        } else {
            normalize();
        }
        
        // release the created resources
        dispose();
        monitor.done();
    }
    
    
    /**
     * Recursively removes subtrees. In other words, removes leafs from the graph until no more
     * leafs are present.
     */
    private void removeSubtrees() {
        
        subtreeNodesStack = new ArrayDeque<>();
        
        // find initial leafs
        Queue<NNode> leafs = Lists.newLinkedList();
        for (NNode node : graph.nodes) {
            if (node.getConnectedEdges().size() == 1) {
                leafs.add(node);
            }
        }
        
        // remove them from the graph like there's no tomorrow
        while (!leafs.isEmpty()) {
            NNode node = leafs.poll();
            // was the edge already removed?
            if (node.getConnectedEdges().size() == 0) {
                continue;
            }
            NEdge edge = node.getConnectedEdges().get(0);
            boolean isOutEdge = node.getOutgoingEdges().size() > 0;
            
            NNode other = edge.getOther(node);
            if (isOutEdge) {
                other.getIncomingEdges().remove(edge);
            } else {
                other.getOutgoingEdges().remove(edge);
            }
            
            if (other.getConnectedEdges().size() == 1) {
                leafs.add(other);
            }
            
            Pair<NNode, NEdge> leafy = Pair.of(node, edge);
            subtreeNodesStack.push(leafy);
            // remove the node from the graph's nodes
            graph.nodes.remove(node);
        }
        
    }
    
    /**
     * Re-attaches the previously removed tree nodes. It is important that 
     * the nodes are re-attached in the opposite order than they were removed.  
     */
    private void reattachSubtrees() {
        
        while (!subtreeNodesStack.isEmpty()) {
            
            Pair<NNode, NEdge> leafy = subtreeNodesStack.pop();
            NNode node = leafy.getFirst();
            NEdge edge = leafy.getSecond();
            
            NNode placed = edge.getOther(node);
            
            if (edge.target == node) {
                placed.getOutgoingEdges().add(edge);
                node.layer = placed.layer + edge.delta;
            } else {
                placed.getIncomingEdges().add(edge);
                node.layer = placed.layer - edge.delta;
            }
            
            graph.nodes.add(node);
        }
        
    }

    /**
     * Helper method for the network simplex layerer. It determines an initial feasible spanning
     * tree of the graph. This graph will be tight by construction. For determination, an initial
     * feasible tree is being computed. If all tree edges contained are tight (i.e. their minimal
     * length corresponds with their actual length), a tight tree has already been found. If not,
     * this method iteratively determines a non-tree edge incident to the tree with a minimal amount
     * of slack (i.e. the edge with the lowest difference between its current and minimal length)
     * and shifts all tree edges accordingly to shorten the edge to its minimal size. The edge has
     * become tight and will be added to the spanning tree together with all tight edges leading to
     * non-tree nodes as well. If all nodes of the graph are contained in the spanning tree, a tight
     * tree has been found. A concluding computation of each edge's initial cut value takes place.
     * 
     * @see NetworkSimplex#tightTreeDFS(NNode) tightTreeDFS()
     */
    private void feasibleTree() {
        
        // determine initial layering
        layeringTopologicalNumbering(sources);
        
        if (edges.size() > 0) {
            Arrays.fill(edgeVisited, false);
            while (tightTreeDFS(graph.nodes.iterator().next()) < graph.nodes.size()) {
                // some nodes are still not part of the tree
                NEdge e = minimalSlack();
                int slack = e.getTarget().layer - e.getSource().layer - e.delta;
                if (e.getTarget().treeNode) {
                    slack = -slack;
                }

                // update tree
                for (NNode node : graph.nodes) {
                    if (node.treeNode) {
                        node.layer += slack;
                    }
                }
                Arrays.fill(edgeVisited, false);
            }
            // update tree-related attributes
            Arrays.fill(edgeVisited, false);
            postorderTraversal(graph.nodes.iterator().next());
            cutvalues();
        }
    }

    /**
     * Helper method for the network simplex layerer. It determines an (initial) feasible layering
     * for the graph by traversing it by a minimal topological numbering. Dependently of
     * the chosen mode indicated by {@code reverse}, this method traverses incoming edges (if
     * {@code reverse = true}), or outgoing edges, if {@code reverse = false}, only. Therefore, this
     * method should only be called with source nodes as argument in the first-mentioned case and
     * only with sink nodes in the latter case.
     * 
     * @param initialRootNodes
     *            the roots of the topological numbering (sources or sinks, depending on the direction)
     */
    private void layeringTopologicalNumbering(final List<NNode> initialRootNodes) {
        
        // initialize the number of incident edges for each node
        int[] incident = new int[graph.nodes.size()];
        for (NNode node : graph.nodes) {
            incident[node.internalId] += node.getIncomingEdges().size();
        }

        LinkedList<NNode> roots = Lists.newLinkedList(initialRootNodes);
        while (!roots.isEmpty()) {
            NNode node = roots.poll();
            
            for (NEdge edge : node.getOutgoingEdges()) {
                NNode target = edge.getTarget();
                target.layer = Math.max(target.layer, node.layer + edge.delta);
                incident[target.internalId]--;
                if (incident[target.internalId] == 0) {
                    roots.add(target);
                }
            }
        }
    }

    /**
     * Helper method for the network simplex layerer. It determines the length of the currently
     * shortest incoming or outgoing edge of the input node.
     * 
     * @param node
     *            the node to determine the length of its shortest incoming or outgoing edge
     * @return a pair containing the length of the shortest incoming (first element) and outgoing
     *         edge (second element) incident to the input node or {@code -1} as the length, if no
     *         such edge is incident
     */
    private Pair<Integer, Integer> minimalSpan(final NNode node) {
        int minSpanOut = Integer.MAX_VALUE;
        int minSpanIn = Integer.MAX_VALUE;
        int currentSpan;

        for (NEdge edge : node.getConnectedEdges()) {
            currentSpan = edge.getTarget().layer - edge.getSource().layer;
            if (edge.getTarget() == node && currentSpan < minSpanIn) {
                minSpanIn = currentSpan;
            } else if (currentSpan < minSpanOut) {
                minSpanOut = currentSpan;
            }
        }
        
        if (minSpanIn == Integer.MAX_VALUE) {
            minSpanIn = -1;
        }
        if (minSpanOut == Integer.MAX_VALUE) {
            minSpanOut = -1;
        }

        return new Pair<Integer, Integer>(minSpanIn, minSpanOut);
    }

    /**
     * Helper method for the network simplex layerer. It determines a DFS-subtree of the graph by
     * traversing tight edges only (i.e. edges whose current length matches their minimal length in
     * the layering) and returns the number of nodes in this. If this number is equal to the total
     * number of nodes in the graph, a tight spanning tree has been determined.
     * 
     * @param node
     *            the root of the DFS-subtree
     * @return the number of nodes in the determined tight DFS-tree
     */
    private int tightTreeDFS(final NNode node) {
        int nodeCount = 1;
        node.treeNode = true;
        NNode opposite = null;
        for (NEdge edge : node.getConnectedEdges()) {
            if (!edgeVisited[edge.internalId]) {
                edgeVisited[edge.internalId] = true;
                opposite = edge.getOther(node);
                if (edge.treeEdge) {
                    // edge is a tree edge already: follow this path
                    nodeCount += tightTreeDFS(opposite);
                } else if (!opposite.treeNode
                        && edge.delta == edge.getTarget().layer
                                - edge.getSource().layer) {
                    // edge is a tight non-tree edge
                    edge.treeEdge = true;
                    treeEdges.add(edge);
                    nodeCount += tightTreeDFS(opposite);
                }
            }
        }
        return nodeCount;
    }

    /**
     * Helper method for the network simplex layerer. It returns the non-tree edge incident on the
     * tree and incident to a non-tree node with a minimal amount of slack (i.e. an edge with the
     * lowest difference between its current and minimal length) or {@code null}, if no such edge
     * exists. Note, that the returned edge's slack is never {@code 0}, since otherwise, the edge
     * would be a tree-edge.
     * 
     * @return a non-tree edge incident on the tree with a minimal amount of slack or {@code null},
     *         if no such edge exists
     */
    private NEdge minimalSlack() {
        int minSlack = Integer.MAX_VALUE;
        NEdge minSlackEdge = null;
        int curSlack;
        for (NEdge edge : edges) {
            if (edge.getSource().treeNode ^ edge.getTarget().treeNode) {
                // edge is non-tree edge and incident on the tree
                curSlack = edge.getTarget().layer - edge.getSource().layer - edge.delta;
                if (curSlack < minSlack) {
                    minSlack = curSlack;
                    minSlackEdge = edge;
                }
            }
        }
        return minSlackEdge;
    }

    /**
     * Helper method for the network simplex layerer. It performs a postorder DFS-traversal of the
     * graph beginning with the input node. Each node will be assigned a unique traversal ID, which
     * will be stored in {@code poID}. Furthermore, the lowest postorder traversal ID of any node in
     * a descending path relative to the input node will be computed and stored in
     * {@code lowestPoID}, which is also the return value of this method.
     * 
     * @param node
     *            the root of the DFS-subtree
     * @return the lowest post-order ID of any descending edge in the depth-first-search
     * 
     * @see NetworkSimplex#poID poID
     * @see NetworkSimplex#lowestPoID lowestPoID
     * @see NetworkSimplex#postOrder postOrder
     */
    private int postorderTraversal(final NNode node) {
        int lowest = Integer.MAX_VALUE;
        for (NEdge edge : node.getConnectedEdges()) {
            if (edge.treeEdge && !edgeVisited[edge.internalId]) {
                edgeVisited[edge.internalId] = true;
                lowest = Math.min(lowest, postorderTraversal(edge.getOther(node)));
            }
        }
        poID[node.internalId] = postOrder;
        lowestPoID[node.internalId] = Math.min(lowest, postOrder++);
        return lowestPoID[node.internalId];
    }

    /**
     * Helper method for the the network simplex layerer. It determines, whether an node is part of
     * the head component of the given edge defined as follows: If the input edge is deleted, the
     * spanning tree breaks into to connected components. The head component is that component,
     * which contains the edge's target node, and the tail component is the component, which
     * contains the edge's source node. Note that a node either belongs to the head or tail
     * component. Therefore, if the node is not part of the head component, it must be part of the
     * tail component and vice versa.
     * 
     * @param node
     *            the node to determine, whether it belongs to the edges head (or tail) component
     * @param edge
     *            the edge to determine, whether the node is in the head (or tail) component
     * @return {@code true}, if node is in the head component or {@code false}, if the node is in
     *         the tail component of the edge
     */
    private boolean isInHead(final NNode node, final NEdge edge) {
        NNode source = edge.getSource();
        NNode target = edge.getTarget();

        if (lowestPoID[source.internalId] <= poID[node.internalId] 
                && poID[node.internalId] <= poID[source.internalId]
                && lowestPoID[target.internalId] <= poID[node.internalId] 
                && poID[node.internalId] <= poID[target.internalId]) {
            // node is in a descending path in the DFS-Tree
            if (poID[source.internalId] < poID[target.internalId]) {
                // root is in the head component
                return false;
            }
            return true;
        }
        if (poID[source.internalId] < poID[target.internalId]) {
            // root is in the head component
            return true;
        }
        return false;
    }

    /**
     * Helper method for the network simplex layerer. It determines the cut value of each tree edge,
     * which is defined as follows: If the edge is deleted, the spanning tree breaks into two
     * connected components, the head component containing the target node of the edge and the tail
     * component containing the source node of the edge. The cut value is the sum of the weights of
     * all edges going from the tail to the head component, including the tree edge itself, minus
     * the sum of the weights of all edges from the head to the tail component.
     * 
     * @see NetworkSimplex#cutvalue cutvalue
     */
    private void cutvalues() {
        // determine incident tree edges for each node
        List<NNode> leafs = Lists.newArrayList();
        int treeEdgeCount;
        for (NNode node : graph.nodes) {
            treeEdgeCount = 0;
            node.unknownCutvalues.clear();
            for (NEdge edge : node.getConnectedEdges()) {
                if (edge.treeEdge) {
                    node.unknownCutvalues.add(edge);
                    treeEdgeCount++;
                }
            }
            if (treeEdgeCount == 1) {
                leafs.add(node);
            }
        }
        
        // determine cut values
        NEdge toDetermine;
        NNode source, target;
        for (NNode node : leafs) {
            while (node.unknownCutvalues.size() == 1) {
                // one tree edge with undetermined cut value is incident
                toDetermine = node.unknownCutvalues.iterator().next();
                cutvalue[toDetermine.internalId] = toDetermine.weight;
                source = toDetermine.getSource();
                target = toDetermine.getTarget();
                for (NEdge edge : node.getConnectedEdges()) {
                    if (!edge.equals(toDetermine)) {
                        if (edge.treeEdge) {
                            // edge is tree edge
                            if (source.equals(edge.getSource())
                                    || target.equals(edge.getTarget())) {
                                // edge has not the same direction as toDetermine
                                cutvalue[toDetermine.internalId] -= cutvalue[edge.internalId] - edge.weight;
                            } else {
                                cutvalue[toDetermine.internalId] += cutvalue[edge.internalId] - edge.weight;
                            }
                        } else {
                            // edge is non-tree edge
                            if (node.equals(source)) {
                                if (edge.getSource().equals(node)) {
                                    cutvalue[toDetermine.internalId] += edge.weight;
                                } else {
                                    cutvalue[toDetermine.internalId] -= edge.weight;
                                }
                            } else {
                                if (edge.getSource().equals(node)) {
                                    cutvalue[toDetermine.internalId] -= edge.weight;
                                } else {
                                    cutvalue[toDetermine.internalId] += edge.weight;
                                }
                            }
                        }
                    }
                }
                
                // remove edge from 'unknownCutvalues'
                source.unknownCutvalues.remove(toDetermine);
                target.unknownCutvalues.remove(toDetermine);
                
                // proceed with next node
                if (source.equals(node)) {
                    node = toDetermine.getTarget();
                } else {
                    node = toDetermine.getSource();
                }
            }
        }
    }

    /**
     * Helper method for the network simplex layerer. It returns a tree edge with a negative cut
     * value or {@code null}, if no such edge exists, meaning that the current layer assignment of
     * all nodes is optimal. Note, that this method returns any edge with a negative cut value. A
     * special preference to an edge with lowest value will not be given.
     * 
     * @return a tree edge with negative cut value or {@code null}, if no such edge exists
     */
    private NEdge leaveEdge() {
        for (NEdge edge : treeEdges) {
            if (edge.treeEdge && cutvalue[edge.internalId] < FUZZY_ST_ZERO) {
                return edge;
            }
        }
        return null;
    }

    /**
     * Helper method for the network simplex layerer. It determines an non-tree edge to replace the
     * given tree edge in the spanning tree. All edges going from the head component to the tail
     * component of the edge will be considered. The edge with a minimal amount of slack (i.e. the
     * lowest difference between its current to its minimal length) will be returned.
     * 
     * @param leave
     *            the tree edge to determine a non-tree edge to be replaced with
     * @return a non-tree edge with a minimal amount of slack to replace the given edge
     * @throws IllegalArgumentException
     *             if the input edge is not a tree edge
     */
    private NEdge enterEdge(final NEdge leave) {
        if (!leave.treeEdge) {    
            throw new IllegalArgumentException("The input edge is not a tree edge.");
        }

        NEdge replace = null;
        int repSlack = Integer.MAX_VALUE;
        int slack;
        NNode source, target;
        for (NEdge edge : edges) {
            source = edge.getSource();
            target = edge.getTarget();
            if (isInHead(source, leave) && !isInHead(target, leave)) {
                // edge is to consider
                slack = target.layer - source.layer - edge.delta;
                if (slack < repSlack) {
                    repSlack = slack;
                    replace = edge;
                }
            }
        }
        return replace;
    }

    /**
     * Helper method for the network simplex layerer. It exchanges the tree-edge {@code leave} by
     * the non-tree edge {@code enter} and updates all values based on the tree (i.e. performs a new
     * postorder DFS-traversal and updates the cut values).
     * 
     * @param leave
     *            the tree-edge to be replaced
     * @param enter
     *            the non-tree edge to replace the tree edge
     * @throws IllegalArgumentException
     *             if either {@code leave} is no tree edge or {@code enter} is a tree edge already
     * 
     * @see NetworkSimplex#enterEdge(NEdge) enterEdge()
     * @see NetworkSimplex#leaveEdge() leaveEdge()
     */
    private void exchange(final NEdge leave, final NEdge enter) {
        if (!leave.treeEdge) {
            throw new IllegalArgumentException("Given leave edge is no tree edge.");
        }
        if (enter.treeEdge) {
            throw new IllegalArgumentException("Given enter edge is a tree edge already.");
        }

        // update tree
        leave.treeEdge = false;
        treeEdges.remove(leave);
        enter.treeEdge = true;
        treeEdges.add(enter);
        int delta = enter.getTarget().layer - enter.getSource().layer - enter.delta;
        if (!isInHead(enter.getTarget(), leave)) {
            delta = -delta;
        }
        for (NNode node : graph.nodes) {
            if (!isInHead(node, leave)) {
                node.layer += delta;
            }
        }
        
        // TODO it should be possible to do this incrementally right?
        // update tree-based values
        postOrder = 1;
        Arrays.fill(edgeVisited, false);
        postorderTraversal(graph.nodes.iterator().next());
        cutvalues();
    }

    /**
     * Helper method for the network simplex layerer. It normalizes the layering, i.e. determines
     * the lowest layer assigned to a node and shifts all nodes up or down in the layers
     * accordingly. After termination, the lowest layer assigned to a node will be zeroth (and
     * therefore first) layer. This method returns an integer array indicating how many nodes are
     * assigned to which layer. Note that the total number of layers necessary to layer the graph is
     * indicated thereby, which is the size if the array.
     * 
     * @return an integer array indicating how many nodes are assigned to which layer
     */
    private int[] normalize() {
        // determine lowest assigned layer and layer count
        int highest = Integer.MIN_VALUE;
        int lowest = Integer.MAX_VALUE;
        for (NNode node : graph.nodes) {
            lowest = Math.min(lowest, node.layer);
            highest = Math.max(highest, node.layer);
        }
        // normalize and determine layer filling
        int[] filling = new int[highest - lowest + 1];
        for (NNode node : graph.nodes) {
            node.layer -= lowest;
            filling[node.layer]++;
        }
        
        // also consider nodes of already layered connected components
        int layerID = 0;
        if (previousLayeringNodeCounts != null) {
            for (int nodeCntInLayer : previousLayeringNodeCounts) {
                filling[layerID++] += nodeCntInLayer;
                if (filling.length == layerID) {
                    break;
                }
            }
        }
        return filling;
    }

    /**
     * Helper method for the network simplex layerer. It balances the layering concerning its width,
     * i.e. the number of nodes in each layer. If the graph allows multiple optimal layerings
     * regarding a minimal edge length, this method moves separate nodes to a layer with a minimal
     * amount of currently contained nodes with respect to the retention of feasibility and
     * optimality of the given layering.
     * 
     * @param filling
     *            an integer array indicating how many nodes are currently assigned to each layer
     */
    private void balance(final int[] filling) {
        // determine possible layers
        int newLayer;
        Pair<Integer, Integer> range = null;
        for (NNode node : graph.nodes) {
            if (node.getIncomingEdges().size() == node.getOutgoingEdges().size()) {
                // node might get shifted
                newLayer = node.layer;
                range = minimalSpan(node);
                for (int i = node.layer - range.getFirst() + 1; i < node.layer
                        + range.getSecond(); i++) {
                    if (filling[i] < filling[newLayer]) {
                        newLayer = i;
                    }
                }
                // assign new layer
                if (filling[newLayer] < filling[node.layer]) {
                    filling[node.layer]--;
                    filling[newLayer]++;
                    node.layer = newLayer;
                }
            }
        }
    }

}
