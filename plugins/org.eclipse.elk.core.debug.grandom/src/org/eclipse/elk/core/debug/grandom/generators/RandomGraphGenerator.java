/*******************************************************************************
 * Copyright (c) 2011, 2020 Kiel University and others.
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
package org.eclipse.elk.core.debug.grandom.generators;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Deque;
import java.util.LinkedList;
import java.util.List;
import java.util.ListIterator;
import java.util.Queue;
import java.util.Random;
import java.util.Set;

import org.eclipse.elk.core.debug.grandom.gRandom.Configuration;
import org.eclipse.elk.core.debug.grandom.generators.GeneratorOptions.RandVal;
import org.eclipse.elk.core.math.ElkMath;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.options.PortConstraints;
import org.eclipse.elk.core.options.PortSide;
import org.eclipse.elk.core.util.Pair;
import org.eclipse.elk.graph.ElkConnectableShape;
import org.eclipse.elk.graph.ElkEdge;
import org.eclipse.elk.graph.ElkLabel;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.ElkPort;
import org.eclipse.elk.graph.impl.ElkNodeImpl;
import org.eclipse.elk.graph.impl.ElkPortImpl;
import org.eclipse.elk.graph.properties.IProperty;
import org.eclipse.elk.graph.util.ElkGraphUtil;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;


/**
 * The random graph generator for ELK graphs.
 */
public class RandomGraphGenerator {

    /** minimal separation of ports. */
    public static final float PORT_SEPARATION = 7.0f;

    /** the maximal number of iterations for distributing edges. */
    private static final int MAX_ITER = 12;

    /** the generator options holder. */
    private GeneratorOptions options;
    /** the random number generator used to generate the graph. */
    private Random random;
    /** the counter used to generate node labels. */
    private int nodeLabelCounter;
    /** the counter used to generate port labels. */
    private int portLabelCounter;

    private int maxHierarchyLevel;
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Constructor
    
    /**
     * Create a random graph generator with given random number generator.
     *
     * @param random
     *            a random number generator
     */
    public RandomGraphGenerator(final Random random) {
        this.random = random;
    }
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Random Graph Generation API

    /**
     * Generate a random graph based on the given generator options.
     */
    public ElkNode generate(final GeneratorOptions opts) {
        // reset the generator
        nodeLabelCounter = 0;
        portLabelCounter = 0;
        this.options = opts;
        // generate the graph
        ElkNode parent = ElkGraphUtil.createGraph();
        maxHierarchyLevel = get(GeneratorOptions.MAX_HIERARCHY_LEVEL).intVal(random);
        if (get(GeneratorOptions.SMALL_HIERARCHY)) {
            List<List<ElkNode>> atomicNodesOnLevels = addHierarchicalNodes(parent);
            connectAtomicNodesOnDifferentLevels(atomicNodesOnLevels);
            return parent;
        } else {
            return makeGraph(parent);
        }
    }
    
    /**
     * Creates a number of random graphs based on the supplied configuration.
     */
    public List<ElkNode> generate(final Configuration config) {
        List<ElkNode> generated = new ArrayList<>();
        
        // Parse configuration
        random = randomGeneratorFromConfig(config);
        GeneratorOptions opts = ConfigurationParser.parse(config, random);
        
        // Generate graphs
        int graphCount = Math.max(1, config.getSamples());
        for (int i = 0; i < graphCount; i++) {
            generated.add(generate(opts));
        }
        
        return generated;
    }
    
    /**
     * Returns a new random number generator initialized with a seed as specified by the supplied configuration.
     */
    private Random randomGeneratorFromConfig(final Configuration config) {
        Integer seed = config.getSeed();
        if (seed == null) {
            return new Random();
        } else {
            return new Random(seed);
        }
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Implementation

    private List<List<ElkNode>> addHierarchicalNodes(final ElkNode parent) {
        // List that contains all the generated children(children of the same node are in the same
        // sublist)
        List<List<ElkNode>> atomicNodes = new ArrayList<>();
        Deque<ElkNode> allGraphs = new LinkedList<ElkNode>();
        
        // the graphs in the actual level(that should get children)
        allGraphs.push(parent);
        for (int i = 0; i < maxHierarchyLevel; i++) {
            List<ElkNode> hierarchicalNodes = new ArrayList<>();
            while (!allGraphs.isEmpty()) {
                ElkNode p = allGraphs.pop();
                makeGraph(p);
                // Must copy list, because eclipse EList cannot be shuffled.
                if (!p.getChildren().isEmpty()) {
                    atomicNodes.add(new ArrayList<>(p.getChildren()));
                }
                // Do not add more hierarchical nodes on lowest level.
                if (i < maxHierarchyLevel - 1) {
                    int numHierarch = get(GeneratorOptions.NUMBER_HIERARCHICAL_NODES).intVal(random);
                    hierarchicalNodes.addAll(createIndependentSet(p, numHierarch));
                }
            }
            for (ElkNode node : hierarchicalNodes) {
                allGraphs.push(node);
            }
        }
        return atomicNodes;
    }

    private void connectAtomicNodesOnDifferentLevels(final List<List<ElkNode>> atomicNodes) {
        int numCrossHier;
        if (get(GeneratorOptions.EXACT_RELATIVE_HIER) != null) {
            Set<ElkEdge> newHashSet = Sets.newHashSet();
            for (List<ElkNode> level : atomicNodes) {
                for (ElkNode node : level) {
                    newHashSet.addAll(node.getIncomingEdges());
                    newHashSet.addAll(node.getOutgoingEdges());
                }
            }
            numCrossHier = (int) (newHashSet.size() * get(GeneratorOptions.EXACT_RELATIVE_HIER).val(random));
        } else {
            numCrossHier = get(GeneratorOptions.CROSS_HIER).intVal(random);
        }
        if (atomicNodes.size() < 2) {
            return;
        }
        for (int i = 0; i < numCrossHier; i++) {
            if (atomicNodes.size() > 1) {
                List<List<ElkNode>> levels = sample(atomicNodes, 2);
                connect(sample(levels.get(0), 1).get(0), sample(levels.get(1), 1).get(0));
            }
        }
    }

    private <T> List<T> sample(final List<T> list, final int number) {
        List<T> shuffleThis = new ArrayList<>(list);
        Collections.shuffle(shuffleThis);
        if (number > list.size()) {
            return shuffleThis;
        }
        return shuffleThis.subList(0, number);
    }

    private ElkNode makeGraph(final ElkNode graph) {
        if (!get(GeneratorOptions.ENABLE_HIERARCHY)) {
            set(GeneratorOptions.HIERARCHY_CHANCE, 0.0f);
            set(GeneratorOptions.CROSS_HIERARCHY_EDGES, false);
        }

        // initialize basic properties
        int n = get(GeneratorOptions.NUMBER_OF_NODES).intVal(random);

        int m;
        switch (get(GeneratorOptions.EDGE_DETERMINATION)) {
        case ABSOLUTE: {
            m = get(GeneratorOptions.EDGES_ABSOLUTE).intVal(random);
            break;
        }

        case RELATIVE: {
            m = (int) (n * get(GeneratorOptions.RELATIVE_EDGES).val(random));
            break;
        }

        case DENSITY: {
            double d = get(GeneratorOptions.DENSITY).val(random);
            m = (int) (Math.round(d * n * (n - 1) / 2));
            break;
        }

        case OUTGOING: {
            double edgesPerNode = get(GeneratorOptions.OUTGOING_EDGES).val(random);
            m = (int) (Math.round(n * edgesPerNode));
            break;
        }

        default:
            throw new IllegalArgumentException("Selected edge determination is not supported.");
        }

        switch (get(GeneratorOptions.GRAPH_TYPE)) {
        case CUSTOM: {
            generateCustomGraph(n, m, graph);
            break;
        }

        case BIPARTITE: {
            generateBipartite(graph, n, m, 0);
            break;
        }

        case TREE: {
            int maxDegree = get(GeneratorOptions.MAX_DEGREE);
            int maxWidth = get(GeneratorOptions.MAX_WIDTH);
            generateTree(graph, n, maxDegree, maxWidth, 0);
            break;
        }

        case BICONNECTED: {
            generateBiconnectedGraph(graph, n, m, 0);
            break;
        }

        case TRICONNECTED: {
            float p1 = random.nextFloat();
            float p2 = 1.0f - p1;
            generateTriconnectedGraph(graph, n, p1, p2, 0);
            break;
        }

        case ACYCLIC_NO_TRANSITIVE_EDGES: {
            boolean planar = get(GeneratorOptions.PLANAR);
            generateANTEGraph(graph, n, m, planar, false, 0);
            break;
        }

        default:
            throw new IllegalArgumentException("Selected graph generator is not supported.");
        }

        // remove isolated nodes if requested
        if (!get(GeneratorOptions.ISOLATED_NODES)) {
            removeIsolatedNodes(graph);
        }

        // if ports require a predefined position, assign it to all ports
        if (get(GeneratorOptions.ENABLE_PORTS) && get(GeneratorOptions.PORT_CONSTRAINTS) == PortConstraints.FIXED_POS) {
            LinkedList<ElkNode> nodeQueue = new LinkedList<ElkNode>();
            nodeQueue.add(graph);
            do {
                ElkNode node = nodeQueue.removeFirst();
                distributePorts(node);
                nodeQueue.addAll(node.getChildren());
            } while (!nodeQueue.isEmpty());
        }

        return graph;
    }

    private List<ElkNode> generateCustomGraph(final int n, final int m, final ElkNode graph) {
        List<ElkNode> generatedNodes;
        switch (get(GeneratorOptions.EDGE_DETERMINATION)) {
        case OUTGOING:
            generatedNodes = generateAnyGraph(graph, n, 0);
            break;
        default:
            generatedNodes = generateAnyGraph(graph, n, m, 0);
        }

        // TODO-alan consider removal
        if (get(GeneratorOptions.CROSS_HIERARCHY_EDGES)) {
            // create edges randomly across the whole compound graph, crossing hierarchy borders
            switch (get(GeneratorOptions.EDGE_DETERMINATION)) {
            case OUTGOING: {
                int[] outgoingEdges = determineOutgoingEdges(generatedNodes, get(GeneratorOptions.OUTGOING_EDGES));
                connectRandomlyAndConditional(generatedNodes, outgoingEdges, basicCondition);
                break;
            }
            default: {
                int createdEdges = 0;
                int iterations = 0;
                do {
                    int[] outgoingEdges = determineOutgoingEdges(generatedNodes, m - createdEdges);
                    createdEdges += connectRandomlyAndConditional(generatedNodes, outgoingEdges,
                            basicCondition);
                    iterations++;
                } while (createdEdges < m && iterations < MAX_ITER);
            }
            }
        }

        return generatedNodes;

    }

    /** the basic condition which cares for self-loops, multi-edges and cycles. */
    private final EdgeCondition basicCondition = new EdgeCondition() {
        public boolean evaluate(final ElkNode node1, final ElkNode node2) {
            if (!get(GeneratorOptions.SELF_LOOPS) && node1 == node2) {
                return false;
            }
            if (!get(GeneratorOptions.MULTI_EDGES) && connected(node1, node2)) {
                return false;
            }
            if (!get(GeneratorOptions.CYCLES) && findNodeWithDFS(node2, node1)) {
                return false;
            }
            return true;
        }
    };

    /**
     * Generates a random graph.
     *
     * @param parent
     *            the parent node
     * @param n
     *            the number of nodes
     * @param m
     *            the number of edges
     * @param hierarchyLevel
     *            the current hierarchy level
     * @return the list of created nodes
     */
    private List<ElkNode> generateAnyGraph(final ElkNode parent, final int n, final int m,
            final int hierarchyLevel) {
        // create the nodes
        List<ElkNode> nodes = createIndependentSet(parent, n);
        // connect the nodes
        if (!get(GeneratorOptions.CROSS_HIERARCHY_EDGES)) {
            // determine the number of outgoing edges for every node
            int createdEdges = 0;
            int iterations = 0;
            do {
                int[] outgoingEdges = determineOutgoingEdges(nodes, m - createdEdges);
                createdEdges += connectRandomlyAndConditional(nodes, outgoingEdges, basicCondition);
                iterations++;
            } while (createdEdges < m && iterations < MAX_ITER);
        }
        // recursively create hierarchy if applicable
        float hierarchyChance = get(GeneratorOptions.HIERARCHY_CHANCE);
        if (hierarchyChance > 0.0f && hierarchyLevel < maxHierarchyLevel) {
            for (ElkNode node : nodes.toArray(new ElkNode[nodes.size()])) {
                if (!isHypernode(node) && random.nextFloat() < hierarchyChance) {
                    // determine the number of nodes and edges in the compound node
                    float sizeFactor = random.nextFloat() * get(GeneratorOptions.HIERARCHY_NODES_FACTOR);
                    int cn = Math.round(sizeFactor * n);
                    if (cn == 0) {
                        cn = 1;
                    }
                    int cm = Math.round(sizeFactor * m);
                    List<ElkNode> childNodes = generateAnyGraph(node, cn, cm, hierarchyLevel + 1);
                    nodes.addAll(childNodes);
                }
            }
        }
        return nodes;
    }

    /**
     * Generates a random graph.
     *
     * @param parent
     *            the parent node
     * @param n
     *            the number of nodes
     * @param minOut
     *            the minimum number of outgoing edges per node
     * @param maxOut
     *            the maximum number of outgoing edges per node
     * @param hierarchyLevel
     *            the current hierarchy level
     * @return the list of created nodes
     */
    private List<ElkNode> generateAnyGraph(final ElkNode parent, final int n,
            final int hierarchyLevel) {
        // create the nodes
        List<ElkNode> nodes = createIndependentSet(parent, n);
        // connect the nodes
        if (!get(GeneratorOptions.CROSS_HIERARCHY_EDGES)) {
            // determine the number of outgoing edges for every node
            int[] outgoingEdges = determineOutgoingEdges(nodes, get(GeneratorOptions.OUTGOING_EDGES));
            connectRandomlyAndConditional(nodes, outgoingEdges, basicCondition);
        }
        // recursively create hierarchy if applicable
        float hierarchyChance = get(GeneratorOptions.HIERARCHY_CHANCE);
        if (hierarchyChance > 0.0f && hierarchyLevel < maxHierarchyLevel) {
            for (ElkNode node : nodes.toArray(new ElkNode[nodes.size()])) {
                if (!isHypernode(node) && random.nextFloat() < hierarchyChance) {
                    // determine the number of nodes in the compound node
                    float sizeFactor = random.nextFloat() * get(GeneratorOptions.HIERARCHY_NODES_FACTOR);
                    int cn = Math.round(sizeFactor * n);
                    if (cn == 0) {
                        cn = 1;
                    }
                    List<ElkNode> childNodes =
                            generateAnyGraph(node, cn, hierarchyLevel + 1, maxHierarchyLevel);
                    nodes.addAll(childNodes);
                }
            }
        }
        return nodes;
    }

    /**
     * Generates a random bipartite graph.
     *
     * @param parent
     *            the parent node
     * @param n
     *            the number of nodes
     * @param m
     *            the number of edges
     * @param minPartition
     *            the minimal fraction of nodes in the second set
     * @param maxPartition
     *            the maximal fraction of nodes in the second set
     * @param hierarchyLevel
     *            the current hierarchy level
     */
    private void generateBipartite(final ElkNode parent, final int n, final int m,
            final int hierarchyLevel) {
        int n2 = ElkMath.boundi((int) Math.round(n * get(GeneratorOptions.PARTITION_FRAC).val(random)), 1, n - 1);
        int n1 = n - n2;
        ElkNode[] nodes1 = new ElkNode[n1];
        for (int i = 0; i < n1; i++) {
            nodes1[i] = createNode(parent);
        }
        ElkNode[] nodes2 = new ElkNode[n2];
        for (int i = 0; i < n2; i++) {
            nodes2[i] = createNode(parent);
        }
        boolean allowCycles = get(GeneratorOptions.CYCLES);
        for (int j = 0; j < m; j++) {
            int source;
            if (allowCycles) {
                source = random.nextInt(n);
            } else {
                source = random.nextInt(n1);
            }
            if (source < n1) {
                int target = random.nextInt(n2);
                connectConditional(nodes1[source], nodes2[target], basicCondition);
            } else {
                int target = random.nextInt(n1);
                connectConditional(nodes2[source - n1], nodes1[target], basicCondition);
            }
        }

        // recursively create hierarchy if applicable
        float hierarchyChance = get(GeneratorOptions.HIERARCHY_CHANCE);
        if (hierarchyChance > 0.0f && hierarchyLevel < maxHierarchyLevel) {
            for (ElkNode node : parent.getChildren()) {
                if (!isHypernode(node) && random.nextFloat() < hierarchyChance) {
                    // determine the number of nodes in the compound node
                    float hierarchyNodesFactor = get(GeneratorOptions.HIERARCHY_NODES_FACTOR);
                    int cn = randomInt(1, (int) (hierarchyNodesFactor * n));
                    int cm = Math.round((float) cn / n * m);
                    generateBipartite(node, cn, cm, hierarchyLevel + 1);
                }
            }
        }
    }

    /**
     * Generates a random tree. The implementation is based upon the one used in the OGDF library.
     *
     * @param parent
     *            the parent node
     * @param n
     *            the number of nodes
     * @param maxDeg
     *            the maximum degree
     * @param maxWidth
     *            the maximum width
     * @param hierarchyLevel
     *            the current hierarchy level
     */
    private void generateTree(final ElkNode parent, final int n, final int maxDeg,
            final int maxWidth, final int hierarchyLevel) {
        int max = 0;
        int nodeIdCounter = 0;
        @SuppressWarnings("unchecked")
        Pair<ElkNode, Integer>[] possible = (Pair<ElkNode, Integer>[]) new Pair[n];
        int[] width = new int[n + 1];
        int[] level = new int[n];
        // create the root node
        ElkNode rootNode = createNode(parent);
        int rootId = nodeIdCounter++;
        possible[0] = new Pair<ElkNode, Integer>(rootNode, rootId);
        level[rootId] = 0;
        // create the tree
        for (int i = 1; i < n;) {
            // get the node to append to
            int x = randomInt(0, max);
            Pair<ElkNode, Integer> nodeInfo = possible[x];
            ElkNode node = nodeInfo.getFirst();
            int nodeId = nodeInfo.getSecond();
            // check for the width constraint
            if (maxWidth != 0 && width[level[nodeId] + 1] == maxWidth) {
                possible[x] = possible[max--];
                continue;
            }
            // check for the out-degree constraint
            if (maxDeg != 0 && node.getOutgoingEdges().size() + 1 == maxDeg) {
                possible[x] = possible[max--];
            }
            // append a new node
            ElkNode newNode = createNode(parent);
            int newNodeId = nodeIdCounter++;
            possible[++max] = new Pair<ElkNode, Integer>(newNode, newNodeId);
            connect(node, newNode);
            level[newNodeId] = level[nodeId] + 1;
            ++width[level[newNodeId]];
            ++i;
        }

        // recursively create hierarchy if applicable
        float hierarchyChance = get(GeneratorOptions.HIERARCHY_CHANCE);
        if (hierarchyChance > 0.0f && hierarchyLevel < maxHierarchyLevel) {
            for (ElkNode node : parent.getChildren()) {
                if (!isHypernode(node) && random.nextFloat() < hierarchyChance) {
                    // determine the number of nodes in the compound node
                    float hierarchyNodesFactor = get(GeneratorOptions.HIERARCHY_NODES_FACTOR);
                    int cn = randomInt(1, (int) (hierarchyNodesFactor * n));
                    generateTree(node, cn, maxDeg, maxWidth, hierarchyLevel + 1);
                }
            }
        }
    }

    /**
     * Generates a biconnected graph. The implementation is based upon the one used in the OGDF
     * library.
     *
     * @param parent
     *            the parent node
     * @param n
     *            the number of nodes
     * @param m
     *            the number of edges
     * @param hierarchyLevel
     *            the current hierarchy level
     */
    // CHECKSTYLEOFF MagicNumber
    private void generateBiconnectedGraph(final ElkNode parent, final int n, final int m,
            final int hierarchyLevel) {
        int realN = Math.max(3, n);
        int realM = Math.max(m, realN);
        // the number of split-edge operations
        int kse = realN - 3;
        // the number of add-edge operations
        int kae = realM - realN;
        ElkNode[] nodes = new ElkNode[realN];
        ElkEdge[] edges = new ElkEdge[realM];
        // start with a triangle
        nodes[0] = createNode(parent);
        nodes[1] = createNode(parent);
        nodes[2] = createNode(parent);
        edges[0] = connect(nodes[0], nodes[1]);
        edges[1] = connect(nodes[1], nodes[2]);
        edges[2] = connect(nodes[2], nodes[0]);
        int nNodes = 3;
        int nEdges = 3;
        // generate the graph
        while (kse + kae > 0) {
            int p = randomInt(1, kse + kae);
            if (p <= kse) {
                // split edge
                ElkEdge edge = edges[randomInt(0, nEdges - 1)];
                Pair<ElkNode, ElkEdge> splitInfo = split(edge);
                nodes[nNodes++] = splitInfo.getFirst();
                edges[nEdges++] = splitInfo.getSecond();
                --kse;
            } else {
                // add edge
                int i = randomInt(0, nNodes - 1);
                int j = (i + randomInt(1, nNodes - 1)) % nNodes;
                edges[nEdges++] = connect(nodes[i], nodes[j]);
                --kae;
            }
        }
        // recursively create hierarchy if applicable
        float hierarchyChance = get(GeneratorOptions.HIERARCHY_CHANCE);
        if (hierarchyChance > 0.0f && hierarchyLevel < maxHierarchyLevel) {
            for (ElkNode node : nodes) {
                if (!isHypernode(node) && random.nextFloat() < hierarchyChance) {
                    // determine the number of nodes in the compound node
                    float hierarchyNodesFactor = get(GeneratorOptions.HIERARCHY_NODES_FACTOR);
                    int cn = randomInt(1, (int) (hierarchyNodesFactor * n));
                    // preserve density for number of edges
                    float density = (float) m / (n * n);
                    int cm = (int) density * cn * cn;
                    generateBiconnectedGraph(node, cn, cm, hierarchyLevel + 1);
                }
            }
        }
    }

    // CHECKSTYLEON MagicNumber

    /**
     * Generates a triconnected graph. The implementation is based upon the one used in the OGDF
     * library.
     *
     * @param parent
     *            the parent node
     * @param n
     *            the number of nodes
     * @param p1
     *            the probability for the first additional edge to be added
     * @param p2
     *            the probability for the second additional edge to be added
     * @param hierarchyLevel
     *            the current hierarchy level
     */
    // CHECKSTYLEOFF MagicNumber
    private void generateTriconnectedGraph(final ElkNode parent, final int n, final float p1,
            final float p2, final int hierarchyLevel) {
        int realN = Math.max(n, 4);
        // start with a clique of size 4
        List<ElkNode> cliqueNodes = createClique(parent, 4);
        // array of all nodes
        ElkNode[] nodes = new ElkNode[realN];
        int i = 0;
        for (ElkNode node : cliqueNodes) {
            nodes[i++] = node;
        }
        // array of neighbors
        ElkEdge[] neighbors = new ElkEdge[realN];
        // neighbor markings
        // 0 = not marked
        // 1 = marked left
        // 2 = marked right
        // 3 = marked both
        int[] marks = new int[n];
        // generate the graph
        for (; i < n; ++i) {
            ElkNode node = nodes[randomInt(0, i - 1)];
            // create a new node to split 'node' in two
            ElkNode newNode = createNode(parent);
            nodes[i] = newNode;
            // build array of all neighbors
            int d = node.getOutgoingEdges().size() + node.getIncomingEdges().size();
            int j = 0;
            for (ElkEdge edge : node.getOutgoingEdges()) {
                neighbors[j++] = edge;
            }
            for (ElkEdge edge : node.getIncomingEdges()) {
                neighbors[j++] = edge;
            }
            // mark two distinct neighbors for left
            for (j = 2; j > 0;) {
                int r = randomInt(0, d - 1);
                if ((marks[r] & 1) == 0) {
                    marks[r] |= 1;
                    --j;
                }
            }
            // mark two distinct neighbors for right
            for (j = 2; j > 0;) {
                int r = randomInt(0, d - 1);
                if ((marks[r] & 2) == 0) {
                    marks[r] |= 2;
                    --j;
                }
            }
            // perform the node-split
            for (j = 0; j < d; ++j) {
                int mark = marks[j];
                marks[j] = 0;
                // decide to with which node each neighbor is connected
                double x = random.nextDouble();
                switch (mark) {
                case 0:
                    if (x < p1) {
                        mark = 1;
                    } else if (x < p1 + p2) {
                        mark = 2;
                    } else {
                        mark = 3;
                    }
                    break;
                case 1:
                case 2:
                    if (x >= p1 + p2) {
                        mark = 3;
                    }
                    break;
                }
                // move edge or create new one if necessary
                ElkEdge edge = neighbors[j];
                switch (mark) {
                case 2:
                    if (edge.getSources().contains(node)) {
                        moveSource(edge, node, newNode);
                    } else {
                        moveTarget(edge, node, newNode);
                    }
                    break;
                case 3:
                    if (edge.getSources().contains(node)) {
                        connect(newNode, edge.getTargets().get(0));
                    } else {
                        connect(newNode, edge.getSources().get(0));
                    }
                    break;
                }
            }
            connect(node, newNode);
        }
        // recursively create hierarchy if applicable
        float hierarchyChance = get(GeneratorOptions.HIERARCHY_CHANCE);
        if (hierarchyChance > 0.0f && hierarchyLevel < maxHierarchyLevel) {
            for (ElkNode node : parent.getChildren()) {
                if (!isHypernode(node) && random.nextFloat() < hierarchyChance) {
                    // determine the number of nodes in the compound node
                    float hierarchyNodesFactor = get(GeneratorOptions.HIERARCHY_NODES_FACTOR);
                    int cn = randomInt(1, (int) (hierarchyNodesFactor * n));
                    generateTriconnectedGraph(node, cn, p1, p2, hierarchyLevel + 1);
                }
            }
        }
    }

    /**
     * Generates an acyclic graph without transitive edges. The implementation is based upon the one
     * used in the OGDF library.
     *
     * @param parent
     *            the parent node
     * @param n
     *            the number of nodes
     * @param m
     *            the number of edges
     * @param planar
     *            whether the generated graph should be planar
     * @param singleSource
     *            whether the graph is a single source graph
     * @param hierarchyLevel
     *            the current hierarchy level
     */
    private void generateANTEGraph(final ElkNode parent, final int n, final int m,
            final boolean planar, final boolean singleSource, final int hierarchyLevel) {
        ElkNode[] nnr = new ElkNode[3 * n];
        int[] vrt = new int[3 * n];
        int[] fst = new int[n + 1];
        List<HierarchyEdge> startEdges = new LinkedList<HierarchyEdge>();
        HierarchyEdge actEdge, nextEdge;
        int act, next, n1, n2, idc = 0;
        boolean connected;
        // create the nodes
        for (int i = 0; i < n; ++i) {
            createNode(parent);
        }
        int numberOfLayers = 0, totNumber = 0, realCount = 0;
        fst[0] = 0;
        for (ElkNode node : parent.getChildren()) {
            nnr[totNumber] = node;
            vrt[totNumber++] = 0;
            realCount++;
            float r = random.nextFloat();
            if (totNumber == 1 && singleSource || realCount == n || r * r * n < 1) {
                fst[++numberOfLayers] = totNumber;
            }
        }
        // determine allowed neighbors
        int[] leftN = new int[totNumber];
        int[] rightN = new int[totNumber];
        for (int l = 1; l < numberOfLayers; l++) {
            if (planar) {
                n1 = fst[l - 1];
                n2 = fst[l];
                leftN[n2] = n1;
                while (n1 < fst[l] && n2 < fst[l + 1]) {
                    float r = random.nextFloat();
                    if (n1 != fst[l] - 1
                            && (n2 == fst[l + 1] - 1 || r < (float) (fst[l] - fst[l - 1])
                                    / (float) (fst[l + 1] - fst[l - 1]))) {
                        n1++;
                    } else {
                        rightN[n2] = n1;
                        if (++n2 < fst[l + 1]) {
                            leftN[n2] = n1;
                        }
                    }
                }
            } else {
                for (n2 = fst[l]; n2 < fst[l + 1]; n2++) {
                    leftN[n2] = fst[l - 1];
                    rightN[n2] = fst[l] - 1;
                }
            }
        }
        // insert edges
        @SuppressWarnings("unchecked")
        List<HierarchyEdge>[] edgeIn = new LinkedList[totNumber];
        @SuppressWarnings("unchecked")
        List<HierarchyEdge>[] edgeOut = new LinkedList[totNumber];
        for (int i = 0; i < totNumber; ++i) {
            edgeIn[i] = new LinkedList<HierarchyEdge>();
            edgeOut[i] = new LinkedList<HierarchyEdge>();
        }
        if (numberOfLayers != 0) {
            float x1 = m;
            float x2 = 0;
            for (n2 = fst[1]; n2 < totNumber; n2++) {
                if (vrt[n2] == 0) {
                    x2 += rightN[n2] - leftN[n2] + 1;
                }
            }
            for (n2 = fst[1]; n2 < totNumber; n2++) {
                if (vrt[n2] == 0) {
                    connected = !singleSource;
                    for (n1 = leftN[n2]; n1 <= rightN[n2] || !connected; n1++) {
                        float r = random.nextFloat();
                        if (r < x1 / x2 || n1 > rightN[n2]) {
                            next = (n1 <= rightN[n2] ? n1 : randomInt(leftN[n2], rightN[n2]));
                            act = n2;
                            nextEdge = new HierarchyEdge(next, act, idc++);
                            while (vrt[next] != 0) {
                                act = next;
                                next = randomInt(leftN[act], rightN[act]);
                                edgeOut[act].add(nextEdge);
                                nextEdge = new HierarchyEdge(next, act, idc++);
                                edgeIn[act].add(nextEdge);
                            }
                            startEdges.add(nextEdge);
                            connected = true;
                            x1 -= 1;
                        }
                        if (n1 <= rightN[n2]) {
                            x2 -= 1;
                        }
                    }
                }
            }
        }
        if (planar) {
            for (act = 0; act < totNumber; act++) {
                Collections.sort(edgeIn[act], new TailComparator());
                Collections.sort(edgeOut[act], new HeadComparator());
            }
        }
        for (act = 0; act < totNumber; act++) {
            List<HierarchyEdge> hedges = edgeIn[act];
            for (HierarchyEdge hedge : hedges) {
                nextEdge = hedge;
                nextEdge.setNext(edgeOut[act].remove(0));
            }
        }
        for (HierarchyEdge hedge : startEdges) {
            actEdge = hedge;
            nextEdge = actEdge;
            while (vrt[nextEdge.getHead()] != 0) {
                nextEdge = nextEdge.getNext();
            }
            connect(nnr[actEdge.getTail()], nnr[nextEdge.getHead()]);
        }
        // recursively create hierarchy if applicable
        float hierarchyChance = get(GeneratorOptions.HIERARCHY_CHANCE);
        if (hierarchyChance > 0.0f && hierarchyLevel < maxHierarchyLevel) {
            for (ElkNode node : parent.getChildren()) {
                if (!isHypernode(node) && random.nextFloat() < hierarchyChance) {
                    // determine the number of nodes in the compound node
                    float hierarchyNodesFactor = get(GeneratorOptions.HIERARCHY_NODES_FACTOR);
                    int cn = randomInt(1, (int) (hierarchyNodesFactor * n));
                    // preserve density for number of edges
                    float density = (float) m / (n * n);
                    int cm = (int) density * cn * cn;
                    generateANTEGraph(node, cn, cm, planar, singleSource, hierarchyLevel + 1);
                }
            }
        }
    }

    // CHECKSTYLEON MagicNumber

    /**
     * A helper class for creating hierarchical graphs.
     */
    private static class HierarchyEdge {

        /** the head, tail and id. */
        private int head, tail, id;

        /** the next edge. */
        private HierarchyEdge next;

        /**
         * Constructs a HierarchyEdge.
         *
         * @param head
         *            the head
         * @param tail
         *            the tail
         * @param id
         *            the id
         */
        HierarchyEdge(final int head, final int tail, final int id) {
            this.head = head;
            this.tail = tail;
            this.id = id;
        }

        /**
         * Returns the head.
         *
         * @return the head
         */
        public int getHead() {
            return head;
        }

        /**
         * Returns the tail.
         *
         * @return the tail
         */
        public int getTail() {
            return tail;
        }

        /**
         * Returns the id.
         *
         * @return the id
         */
        public int getId() {
            return id;
        }

        /**
         * Returns the next edge.
         *
         * @return the next edge
         */
        public HierarchyEdge getNext() {
            return next;
        }

        /**
         * Sets the next edge.
         *
         * @param next
         *            the next edge
         */
        public void setNext(final HierarchyEdge next) {
            this.next = next;
        }
    }

    /**
     * Compares hierarchy edges by the id attribute.
     *
     * @param edge1
     *            the first edge
     * @param edge2
     *            the second edge
     * @return the result of the comparison
     */
    private static int compareId(final HierarchyEdge edge1, final HierarchyEdge edge2) {
        return edge1.getId() < edge2.getId() ? -1 : (edge1.getId() > edge2.getId() ? 1 : 0);
    }

    /**
     * A helper class for comparing hierarchy edges by the head attribute.
     */
    private static class HeadComparator implements Comparator<HierarchyEdge> {

        @Override
        public int compare(final HierarchyEdge edge1, final HierarchyEdge edge2) {
            return edge1.getHead() < edge2.getHead() ? -1
                    : (edge1.getHead() > edge2.getHead() ? 1 : compareId(edge1, edge2));
        }
    }

    /**
     * A helper class for comparing hierarchy edges by the tail attribute.
     */
    private static class TailComparator implements Comparator<HierarchyEdge> {

        @Override
        public int compare(final HierarchyEdge edge1, final HierarchyEdge edge2) {
            return edge1.getTail() < edge2.getTail() ? -1
                    : (edge1.getTail() > edge2.getTail() ? 1 : compareId(edge1, edge2));
        }
    }

    /**
     * Creates an independent set of nodes.
     *
     * @param parent
     *            the parent node
     * @param n
     *            the number of nodes
     * @return the list of created nodes
     */
    private List<ElkNode> createIndependentSet(final ElkNode parent, final int n) {
        List<ElkNode> nodes = new ArrayList<ElkNode>(n);
        for (int i = 0; i < n; ++i) {
            ElkNode node = createNode(parent);
            nodes.add(node);
        }
        return nodes;
    }

    /**
     * Creates a clique.
     *
     * @param parent
     *            the parent node
     * @param n
     *            the number of nodes
     * @return the list of created nodes
     */
    private List<ElkNode> createClique(final ElkNode parent, final int n) {
        List<ElkNode> nodes = createIndependentSet(parent, n);
        for (int i = 0; i < n - 1; ++i) {
            for (int j = i + 1; j < n; ++j) {
                connect(nodes.get(i), nodes.get(j));
            }
        }
        return nodes;
    }

    /**
     * Creates a node.
     *
     * @param parent
     *            the parent node
     * @return the node
     */
    private ElkNode createNode(final ElkNode parent) {
        ElkNode node = ElkGraphUtil.createNode(null);
        float hypernodeChance = get(GeneratorOptions.HYPERNODE_CHANCE);
        if (hypernodeChance > 0.0f && random.nextFloat() < hypernodeChance) {
            node.setProperty(CoreOptions.HYPERNODE, true);
        }

        // create label and identifier
        String nodeid = String.valueOf(nodeLabelCounter++);
        if (get(GeneratorOptions.CREATE_NODE_LABELS)) {
            ElkGraphUtil.createLabel("N" + nodeid, node);
        }
        node.setIdentifier("n" + nodeid);

        // set size of the node
        if (get(GeneratorOptions.SET_NODE_SIZE)) {
            node.setWidth(get(GeneratorOptions.NODE_WIDTH).floatVal(random));
            node.setHeight(get(GeneratorOptions.NODE_HEIGHT).floatVal(random));
        }

        // set port constraints
        PortConstraints portConstraints = get(GeneratorOptions.PORT_CONSTRAINTS);
        if (portConstraints != PortConstraints.UNDEFINED) {
            node.setProperty(CoreOptions.PORT_CONSTRAINTS, portConstraints);
        }

        parent.getChildren().add(node);
        return node;
    }

    /**
     * Creates a port.
     *
     * @param node
     *            the containing node
     * @param source
     *            {@code true} if the port will be used as a source port, {@code false} if it will
     *            be used as a target port.
     * @return the port
     */
    private ElkPort createPort(final ElkNode node, final boolean source) {
        ElkPort port = ElkGraphUtil.createPort(null);
        node.getPorts().add(port);

        // create label and identifier
        String portId = String.valueOf(portLabelCounter++);
        if (get(GeneratorOptions.CREATE_PORT_LABELS)) {
            ElkGraphUtil.createLabel("P" + portId, port);
        }
        port.setIdentifier("p" + portId);

        // set size of the port
        if (get(GeneratorOptions.SET_PORT_SIZE)) {
            port.setWidth(get(GeneratorOptions.PORT_WIDTH).floatVal(random));
            port.setHeight(get(GeneratorOptions.PORT_HEIGHT).floatVal(random));
        }

        if (get(GeneratorOptions.PORT_CONSTRAINTS).isSideFixed()) {
            // determine a random node side
            int northProb, eastProb, southProb, westProb;
            if (source) {
                northProb = get(GeneratorOptions.OUTGOING_NORTH_SIDE);
                eastProb = get(GeneratorOptions.OUTGOING_EAST_SIDE);
                southProb = get(GeneratorOptions.OUTGOING_SOUTH_SIDE);
                westProb = get(GeneratorOptions.OUTGOING_WEST_SIDE);
            } else {
                northProb = get(GeneratorOptions.INCOMING_NORTH_SIDE);
                eastProb = get(GeneratorOptions.INCOMING_EAST_SIDE);
                southProb = get(GeneratorOptions.INCOMING_SOUTH_SIDE);
                westProb = get(GeneratorOptions.INCOMING_WEST_SIDE);
            }
            int p = random.nextInt(northProb + eastProb + southProb + westProb);
            PortSide portSide;
            if (p < northProb) {
                portSide = PortSide.NORTH;
            } else if (p < northProb + eastProb) {
                portSide = PortSide.EAST;
            } else if (p < northProb + eastProb + southProb) {
                portSide = PortSide.SOUTH;
            } else {
                portSide = PortSide.WEST;
            }
            port.setProperty(CoreOptions.PORT_SIDE, portSide);
        }

        return port;
    }

    /**
     * Connects two nodes with an edge.
     *
     * @param source
     *            the source node
     * @param target
     *            the target node
     * @param directed
     *            whether the edge should be directed or undirected
     * @return the edge
     */
    private ElkEdge connect(final ElkConnectableShape source, final ElkConnectableShape target) {

        ElkNode sourceNode = null;
        ElkNode targetNode = null;
        ElkPort sourcePort = null;
        ElkPort targetPort = null;

        if (source instanceof ElkNode) {
            sourceNode = (ElkNode) source;
        }
        if (source instanceof ElkPort) {
            sourcePort = (ElkPort) source;
            sourceNode = sourcePort.getParent();
        }
        if (target instanceof ElkPort) {
            targetPort = (ElkPort) target;
            targetNode = targetPort.getParent();
        }
        if (target instanceof ElkNode) {
            targetNode = (ElkNode) target;
        }

        ElkEdge edge = ElkGraphUtil.createEdge(null);

        if (get(GeneratorOptions.ENABLE_PORTS)) {
            if (source instanceof ElkNode) {
                sourcePort = retrievePort(sourceNode, true);
            }

            if (!isHypernode(sourceNode)) {
                // remove node
                if (edge.getSources().contains(sourceNode)) {
                    edge.getSources().remove(sourceNode);
                }
                edge.getSources().add(sourcePort);
                // edge.setSourcePort(sourcePort);
            }

            if (!isHypernode(targetNode)) {
                if (targetPort == null) {
                    targetPort = retrievePort(targetNode, false);
                }
                if (edge.getTargets().contains(targetNode)) {
                    edge.getTargets().remove(targetNode);
                }
                edge.getTargets().remove(target);
                edge.getTargets().add(targetPort);
                // edge.setTargetPort(targetPort);
            }

        } else

        {
            edge.getSources().add(sourceNode);
            edge.getTargets().add(targetNode);
        }

        if (

        get(GeneratorOptions.EDGE_LABELS)) {
            addEdgeLabel(edge);
        }

        ElkGraphUtil.updateContainment(edge);
        return edge;
    }

    /**
     * Creates a {@link KLabel} with random text and adds it to the edge.
     *
     * @param edge
     *            the edge to add the label to.
     */
    private void addEdgeLabel(final ElkEdge edge) {
        ElkLabel label = ElkGraphUtil.createLabel(BlindTextGenerator.generate(), edge);
        edge.getLabels().add(label);
    }

    /**
     * Retrieves a port for a new edge to connect to the given node through. This can either be a
     * newly created port, or an existing one. Which one it is depends on the chance of ports to be
     * reused.
     *
     * <p>
     * An outgoing edge will only ever try to reuse ports that only have outgoing edges connected to
     * them. The same is true for incoming edges and ports with only incoming edges.
     * </p>
     *
     * @param node
     *            the node to add the port to.
     * @param source
     *            {@code true} if the port will be used as a source port, {@code false} if it will
     *            be used as a target port.
     * @return the new or existing port.
     */
    private ElkPort retrievePort(final ElkNode node, final boolean source) {
        // We might want to reuse an existing port
        float reusePortsChance = get(GeneratorOptions.USE_EXISTING_PORTS_CHANCE).floatVal(random);
        if (reusePortsChance > 0.0f && random.nextFloat() < reusePortsChance) {
            // Collect candidate ports for reuse
            List<ElkPort> reuseCandidates = Lists.newLinkedList();

            for (ElkPort port : node.getPorts()) {
                // Two flags indicating whether the port already has edges pointing in the right
                // or wrong direction connected to it
                boolean connectedToDesiredEdges = false;
                boolean connectedToBadEdges = false;

                // for (ElkEdge edge : port.getEdges()) {
                // connectedToDesiredEdges = (source && edge.getSourcePort() == port)
                // || (!source && edge.getTargetPort() == port);
                // connectedToBadEdges = (source && edge.getTargetPort() == port)
                // || (!source && edge.getSourcePort() == port);
                // }
                // TODO -apo- not exactly the same, but makes more sense...
                connectedToDesiredEdges = (!port.getIncomingEdges().isEmpty() && !source)
                        || (!port.getOutgoingEdges().isEmpty() && source);

                connectedToBadEdges = (!port.getIncomingEdges().isEmpty() && source)
                        || (!port.getOutgoingEdges().isEmpty() && !source);
                // If there are only edges pointing in the same direction as the new edge connected
                // to
                // the port, it qualifies as a candidate for reuse
                if (connectedToDesiredEdges && !connectedToBadEdges) {
                    reuseCandidates.add(port);
                }
            }

            // If there are candidates for reuse, choose one at random
            if (!reuseCandidates.isEmpty()) {
                return reuseCandidates.get(randomInt(0, reuseCandidates.size() - 1));
            }
        }

        // We were unable to reuse an existing port, so create a new one and return that
        return createPort(node, source);
    }

    /**
     * Connects two nodes with an edge if the given condition is evaluated to true.
     *
     * @param source
     *            the source node
     * @param target
     *            the target node
     * @param condition
     *            the condition
     * @return whether the nodes have been connected
     */
    private boolean connectConditional(final ElkNode source, final ElkNode target,
            final EdgeCondition condition) {
        if (condition.evaluate(source, target)) {
            connect(source, target);
            return true;
        }
        return false;
    }

    /**
     * Changes the source of a given edge to a given node.
     *
     * @param edge
     *            the edge
     * @param node
     *            the new source node
     */
    private void moveSource(final ElkEdge edge, final ElkConnectableShape oldConnectShape,
            final ElkConnectableShape newConnectShape) {

        ElkPort port = null;
        ElkNode oldNode = null;
        ElkNode newNode = null;
        ElkPort oldPort = null;
        ElkPort newPort = null;
        if (get(GeneratorOptions.ENABLE_PORTS)) {

            if (!edge.getSources().isEmpty()) {
                if (oldConnectShape instanceof ElkNodeImpl) {
                    oldNode = (ElkNode) oldConnectShape;
                    for (ElkConnectableShape source : edge.getSources()) {
                        if (source instanceof ElkPortImpl) {
                            port = (ElkPortImpl) source;
                            if (port.getParent().equals(oldNode)) {
                                edge.getSources().remove(port);
                                if (port.getOutgoingEdges().size() == 0
                                        && port.getIncomingEdges().size() == 1) {
                                    port.getParent().getPorts().remove(port);
                                }
                            }
                        }
                    }
                } else if (oldConnectShape instanceof ElkPortImpl) {
                    oldPort = (ElkPort) oldConnectShape;
                    if (edge.getSources().contains(oldPort)) {
                        edge.getSources().remove(oldPort);
                    }
                    if (oldPort.getOutgoingEdges().size() == 0
                            && oldPort.getIncomingEdges().size() == 1) {
                        oldPort.getParent().getPorts().remove(oldPort);
                    }
                }
            }

            if (newConnectShape instanceof ElkNodeImpl) {
                newNode = (ElkNode) oldConnectShape;
                newPort = retrievePort(newNode, true);
            } else if (newConnectShape instanceof ElkPortImpl) {
                newPort = (ElkPort) oldConnectShape;
            }
            if (!isHypernode(newNode)) {
                // edge.setSourcePort(newPort);
                edge.getSources().add(newPort);
            }
        } else {

            // edge.setSource(node);
            if (oldConnectShape instanceof ElkNodeImpl) {
                oldNode = (ElkNode) oldConnectShape;
            } else if (oldConnectShape instanceof ElkPortImpl) {
                oldPort = (ElkPort) oldConnectShape;
                oldNode = oldPort.getParent();
            }
            if (edge.getSources().contains(oldNode)) {
                edge.getSources().remove(oldNode);
            }

            if (newConnectShape instanceof ElkNodeImpl) {
                newNode = (ElkNode) newConnectShape;
            } else if (newConnectShape instanceof ElkPortImpl) {
                newPort = (ElkPort) newConnectShape;
                newNode = newPort.getParent();
            }
            edge.getSources().add(newNode);
        }
        ElkGraphUtil.updateContainment(edge);
    }

    /**
     * Changes the target of a given edge to a given node.
     *
     * @param edge
     *            the edge
     * @param node
     *            the new target node
     */
    private void moveTarget(final ElkEdge edge, final ElkConnectableShape oldConnectShape,
            final ElkConnectableShape newConnectShape) {
        ElkPort port = null;
        ElkNode oldNode = null;
        ElkNode newNode = null;
        ElkPort oldPort = null;
        ElkPort newPort = null;
        if (get(GeneratorOptions.ENABLE_PORTS)) {
            // Check if we need to remove the old target port
            // if (edge.getTargetPort() != null && edge.getTargetPort().getEdges().size() == 1) {
            // edge.getTarget().getPorts().remove(edge.getTargetPort());
            // }
            if (!edge.getTargets().isEmpty()) {
                if (oldConnectShape instanceof ElkNodeImpl) {
                    oldNode = (ElkNode) oldConnectShape;
                    for (ElkConnectableShape target : edge.getTargets()) {
                        if (target instanceof ElkPortImpl) {
                            port = (ElkPortImpl) target;
                            if (port.getParent().equals(oldNode)) {
                                edge.getTargets().remove(port);
                                if (port.getOutgoingEdges().size() == 0
                                        && port.getIncomingEdges().size() == 1) {
                                    port.getParent().getPorts().remove(port);
                                }
                            }
                        }
                    }
                } else if (oldConnectShape instanceof ElkPortImpl) {
                    oldPort = (ElkPort) oldConnectShape;
                    if (edge.getTargets().contains(oldPort)) {
                        edge.getTargets().remove(oldPort);
                    }
                    if (oldPort.getOutgoingEdges().size() == 0
                            && oldPort.getIncomingEdges().size() == 1) {
                        oldPort.getParent().getPorts().remove(oldPort);
                    }
                }
            }

            if (newConnectShape instanceof ElkNodeImpl) {
                newNode = (ElkNode) oldConnectShape;
                newPort = retrievePort(newNode, true);
            } else if (newConnectShape instanceof ElkPortImpl) {
                newPort = (ElkPort) oldConnectShape;
            }
            if (!isHypernode(newNode)) {
                // edge.setSourcePort(newPort);
                edge.getTargets().add(newPort);
            }
        } else {

            // edge.setSource(node);
            if (oldConnectShape instanceof ElkNodeImpl) {
                oldNode = (ElkNode) oldConnectShape;
            } else if (oldConnectShape instanceof ElkPortImpl) {
                oldPort = (ElkPort) oldConnectShape;
                oldNode = oldPort.getParent();
            }
            if (edge.getTargets().contains(oldNode)) {
                edge.getTargets().remove(oldNode);
            }

            if (newConnectShape instanceof ElkNodeImpl) {
                newNode = (ElkNode) newConnectShape;
            } else if (newConnectShape instanceof ElkPortImpl) {
                newPort = (ElkPort) newConnectShape;
                newNode = newPort.getParent();
            }
            edge.getTargets().add(newNode);
        }
        ElkGraphUtil.updateContainment(edge);
    }

    /**
     * Splits an edge by inserting a new node and a new edge. The old edge will have all the old
     * sources and just the new node as a target and the new edge will have all the old targets and
     * just the new node as a source. The parent of the new node will be the node containing a whole
     * edge.
     *
     * @param edge
     *            the edge
     * @return a pair containing the new node and the new edge
     */
    private Pair<ElkNode, ElkEdge> split(final ElkEdge edge) {
        ElkNode newNode = createNode(edge.getContainingNode());
        ElkEdge newEdge = connect(newNode, edge.getTargets().get(0));
        ElkConnectableShape oldTarget = edge.getTargets().get(0);
        edge.getTargets().remove(0);
        newEdge.getTargets().addAll(edge.getTargets());
        edge.getTargets().clear();
        edge.getTargets().add(oldTarget);
        moveTarget(edge, edge.getTargets().get(0), newNode);
        ElkGraphUtil.updateContainment(newEdge);
        return new Pair<ElkNode, ElkEdge>(newNode, newEdge);
    }

    /**
     * Connects a source node a number of times to randomly selected nodes of a given list if the
     * condition evaluates to true for the selected node.
     *
     * @param source
     *            the source node
     * @param targets
     *            the target nodes
     * @param number
     *            the number of times the node has to be connected to random nodes
     * @param condition
     *            the condition
     * @return the number of edges which could be inserted
     */
    private int connectRandomlyAndConditional(final ElkNode source, final List<ElkNode> targets,
            final int number, final EdgeCondition condition) {
        ElkNode[] targetBuffer = targets.toArray(new ElkNode[0]);
        int edges = 0;
        int bufferEnd = targetBuffer.length - 1;
        // try connecting the source to up to 'number' nodes randomly
        while (edges < number && bufferEnd >= 0) {
            int i = randomInt(0, bufferEnd);
            ElkNode target = targetBuffer[i];
            if (connectConditional(source, target, condition)) {
                edges++;
            } else {
                // the current node does not fulfill the condition so replace it with an element
                // from the end of the buffer
                targetBuffer[i] = targetBuffer[bufferEnd];
                bufferEnd--;
            }

        }
        return edges;
    }

    /**
     * Connects every node in a list of nodes with random nodes of the same list for the specified
     * number of times.
     *
     * @param nodes
     *            the list of nodes
     * @param outgoingEdges
     *            the number of outgoing edges for every node
     * @param condition
     *            the condition
     * @return the number of edges which could be inserted
     */
    private int connectRandomlyAndConditional(final List<ElkNode> nodes, final int[] outgoingEdges,
            final EdgeCondition condition) {
        // connect every node to the specified number of other nodes
        int edges = 0;
        for (int i = 0; i < nodes.size(); i++) {
            ElkNode source = nodes.get(i);
            edges += connectRandomlyAndConditional(source, nodes, outgoingEdges[i], condition);
        }
        return edges;
    }

    private int[] determineOutgoingEdges(final List<ElkNode> nodes, final RandVal val) {
        // determine the number of outgoing edges for every node
        int n = nodes.size();
        int[] numberOfEdges = new int[n];
        for (int i = 0; i < n; ++i) {
            int c = val.intVal(random);
            numberOfEdges[i] = c;
        }
        return numberOfEdges;
    }

    /**
     * Randomly calculates a number of outgoing edges for every node in a list until a given number
     * of edges have been inserted.
     *
     * @param nodes
     *            the list of nodes
     * @param m
     *            the number of edges
     * @return the number of outgoing edges for every node
     */
    private int[] determineOutgoingEdges(final List<ElkNode> nodes, final int m) {
        // determine the number of outgoing edges for every node
        int[] outgoingEdges = new int[nodes.size()];
        for (int c = 0; c < m; ++c) {
            int i = randomInt(0, nodes.size() - 1);
            ++outgoingEdges[i];
        }
        return outgoingEdges;
    }

    /**
     * Remove all nodes from the graph that have no incoming or outgoing connections and have no
     * child nodes.
     *
     * @param parent
     *            the parent node of the graph
     */
    private static void removeIsolatedNodes(final ElkNode parent) {
        ListIterator<ElkNode> nodeIter = parent.getChildren().listIterator();
        while (nodeIter.hasNext()) {
            ElkNode node = nodeIter.next();
            removeIsolatedNodes(node);
            if (node.getIncomingEdges().isEmpty() && node.getOutgoingEdges().isEmpty()
                    && node.getChildren().isEmpty()) {
                nodeIter.remove();
            }
        }
    }

    /**
     * Returns whether a node can be reached from another node. PRECONDITION: the graph containing
     * the nodes has to be acyclic! If that condition is not met the method is likely to loop
     * infinitely!
     *
     * @param start
     *            the start node
     * @param end
     *            the end node
     * @return whether the end node can be reached from the start node
     */
    private static boolean findNodeWithDFS(final ElkNode start, final ElkNode end) {
        Queue<ElkNode> nodes = new LinkedList<ElkNode>();
        nodes.add(start);
        do {
            ElkNode node = nodes.poll();
            if (node == end) {
                return true;
            }
            for (ElkEdge edge : node.getOutgoingEdges()) {
                for (ElkConnectableShape target : edge.getTargets()) {
                    if (target instanceof ElkNode) {
                        nodes.add((ElkNode) target);
                    } else if (target instanceof ElkPort) {
                        ElkPort port = (ElkPort) target;
                        nodes.add(port.getParent());
                    }
                }
            }
        } while (!nodes.isEmpty());
        return false;
    }

    /**
     * Returns whether two nodes are connected.
     *
     * @param node1
     *            the first node
     * @param node2
     *            the second node
     * @return whether the two nodes are connected
     */
    private static boolean connected(final ElkNode node1, final ElkNode node2) {
        for (ElkEdge edge : node1.getOutgoingEdges()) {
            if (edge.getTargets().contains(node2)) {
                return true;
            }
        }
        for (ElkEdge edge : node2.getOutgoingEdges()) {
            if (edge.getTargets().contains(node1)) {
                return true;
            }
        }
        return false;
    }

    /**
     * Distribute all ports on the border of the given node.
     *
     * @param node
     *            a node
     */
    private void distributePorts(final ElkNode node) {
        // count the ports on each side of the node
        int northCount = 0, eastCount = 0, southCount = 0, westCount = 0;
        for (ElkPort port : node.getPorts()) {
            switch (port.getProperty(CoreOptions.PORT_SIDE)) {
            case NORTH:
                northCount++;
                break;
            case EAST:
                eastCount++;
                break;
            case SOUTH:
                southCount++;
                break;
            case WEST:
                westCount++;
                break;
            }
        }

        // make sure the node is big enough to contain all ports
        float portWidth = get(GeneratorOptions.PORT_WIDTH).floatVal(random);
        float portHeight = get(GeneratorOptions.PORT_HEIGHT).floatVal(random);
        if (node.getWidth() < (northCount + 1) * (portWidth + PORT_SEPARATION)) {
            node.setWidth((northCount + 1) * (portWidth + PORT_SEPARATION));
        }
        if (node.getWidth() < (southCount + 1) * (portWidth + PORT_SEPARATION)) {
            node.setWidth((southCount + 1) * (portWidth + PORT_SEPARATION));
        }
        if (node.getHeight() < (eastCount + 1) * (portHeight + PORT_SEPARATION)) {
            node.setHeight((eastCount + 1) * (portHeight + PORT_SEPARATION));
        }
        if (node.getHeight() < (westCount + 1) * (portHeight + PORT_SEPARATION)) {
            node.setHeight((westCount + 1) * (portHeight + PORT_SEPARATION));
        }

        // distribute the ports on each node side
        double northDelta = node.getWidth() / (northCount + 1);
        double eastDelta = node.getHeight() / (eastCount + 1);
        double southDelta = node.getWidth() / (southCount + 1);
        double westDelta = node.getHeight() / (westCount + 1);
        double northPos = 0, eastPos = 0, southPos = 0, westPos = 0;
        for (ElkPort port : node.getPorts()) {
            switch (port.getProperty(CoreOptions.PORT_SIDE)) {
            case NORTH:
                northPos += northDelta;
                port.setLocation(northPos - port.getWidth() / 2, -port.getHeight());
                break;
            case EAST:
                eastPos += eastDelta;
                port.setLocation(node.getWidth(), eastPos - port.getHeight() / 2);
                break;
            case SOUTH:
                southPos += southDelta;
                port.setLocation(southPos - port.getWidth() / 2, node.getHeight());
                break;
            case WEST:
                westPos += westDelta;
                port.setLocation(-port.getWidth(), westPos - port.getHeight() / 2);
                break;
            }
        }
    }

    /**
     * Returns a random integer number in the given range (including the boundaries).
     *
     * @param from
     *            the minimal number
     * @param to
     *            the maximal number
     * @return a random integer number
     */
    private int randomInt(final int from, final int to) {
        assert from <= to;
        return from + random.nextInt(to - from + 1);
    }

    /**
     * Determine whether the given node is a hypernode.
     *
     * @param node
     *            a node
     * @return true if the node is a hypernode
     */
    private static boolean isHypernode(final ElkNode node) {
        return node.getProperty(CoreOptions.HYPERNODE);
    }

    private <T> T get(final IProperty<T> property) {
        return options.getProperty(property);
    }

    private <T> void set(final IProperty<T> property, final T val) {
        options.setProperty(property, val);
    }

    /**
     * An interface for expressing conditions for creating an edge between two nodes.
     */
    private interface EdgeCondition {
        /**
         * Returns whether the condition is met for an edge from the first to the second node.
         *
         * @param node1
         *            the first node
         * @param node2
         *            the second node
         * @return true if the condition for the edge is met; false else
         */
        boolean evaluate(final ElkNode node1, final ElkNode node2);
    }

}
