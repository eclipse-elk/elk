/*******************************************************************************
 * Copyright (c) 2010, 2020 Kiel University and others.
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
package org.eclipse.elk.alg.layered.p5edges.orthogonal;

import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.eclipse.elk.alg.layered.DebugUtil;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.PortType;
import org.eclipse.elk.alg.layered.p5edges.orthogonal.direction.BaseRoutingDirectionStrategy;
import org.eclipse.elk.alg.layered.p5edges.orthogonal.direction.RoutingDirection;
import org.eclipse.elk.core.options.PortSide;
import org.eclipse.elk.core.util.IElkProgressMonitor;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

/**
 * Edge routing implementation that creates orthogonal bend points. Inspired by:
 * <ul>
 *   <li>Georg Sander. Layout of directed hypergraphs with orthogonal hyperedges. In
 *     <i>Proceedings of the 11th International Symposium on Graph Drawing (GD '03)</i>,
 *     volume 2912 of LNCS, pp. 381-386. Springer, 2004.</li>
 *   <li>Giuseppe di Battista, Peter Eades, Roberto Tamassia, Ioannis G. Tollis,
 *     <i>Graph Drawing: Algorithms for the Visualization of Graphs</i>,
 *     Prentice Hall, New Jersey, 1999 (Section 9.4, for cycle breaking in the
 *     hyperedge segment graph)
 * </ul>
 *
 * <p>This is a generic implementation that can be applied to all four routing directions.
 * Usually, edges will be routed from west to east. However, with northern and southern
 * external ports, this changes: edges are routed from south to north and north to south,
 * respectively. To support these different requirements, the routing direction-related
 * code is factored out into {@link IRoutingDirectionStrategy routing strategies}.</p>
 *
 * <p>When instantiating a new routing generator, the concrete directional strategy must be
 * specified. Once that is done, {@link #routeEdges(LGraph, List, int, List, double)}
 * is called repeatedly to route edges between given lists of nodes.</p>
 */
public final class OrthogonalRoutingGenerator {

    ///////////////////////////////////////////////////////////////////////////////
    // Constants and Variables

    /** differences below this tolerance value are treated as zero. */
    public static final double TOLERANCE = 1e-3;
    
    /** a special return value used by the conflict counting method. */
    private static final int CRITICAL_CONFLICTS_DETECTED = -1;

    /** factor for edge spacing used to determine the {@link #conflictThreshold} (determined experimentally). */
    private static final double CONFLICT_THRESHOLD_FACTOR = 0.5;
    /** factor to compute {@link #criticalConflictThreshold} (determined experimentally). */
    private static final double CRITICAL_CONFLICT_THRESHOLD_FACTOR = 0.2;
    
    /** weight penalty for (non-critical) conflicts. */
    private static final int CONFLICT_PENALTY = 1;
    /** weight penalty for crossings. */
    private static final int CROSSING_PENALTY = 16;
    
    /** we'll be using this thing to split hyper edge segments, if necessary. */
    private HyperEdgeSegmentSplitter segmentSplitter;

    /** routing direction strategy. */
    private final BaseRoutingDirectionStrategy routingStrategy;
    /** spacing between edges. */
    private final double edgeSpacing;
    
    /**
     * Threshold at which horizontal line segments are considered to be too close to one another. A conflict will cause
     * the involved hyperedges to be assigned to different routing slots in order to resolve the conflict.
     * 
     * <p>
     * This value is a constant fraction of the edge spacing.
     * </p>
     */
    private final double conflictThreshold;
    /**
     * Threshold at which horizontal line segments are considered to overlap. A critical conflict will cause the
     * involved hyperedges to be assigned to different routing slots in order to resolve the critical conflict. If there
     * is a cycle of critical conflicts, there is no assignment that will resolve all overlaps. This will result in one
     * of the involved hyperedges to be split into two vertical segments.
     * 
     * <p>
     * This value depends on how close edges from the source layer and edges from the target layer come to each other,
     * respectively. We cannot change that. The critical conflict threshold is thus a fraction of that minimum
     * proximity and is different for each pair of layers. Hence, this is not a constant.
     * </p>
     */
    private double criticalConflictThreshold;
    
    /** prefix of debug output files. */
    private final String debugPrefix;


    ///////////////////////////////////////////////////////////////////////////////
    // Constructor

    /**
     * Constructs a new instance.
     *
     * @param direction
     *            the direction edges should point at.
     * @param edgeSpacing
     *            the space between edges.
     * @param debugPrefix
     *            prefix of debug output files, or {@code null} if no debug output should be generated.
     */
    public OrthogonalRoutingGenerator(final RoutingDirection direction, final double edgeSpacing,
            final String debugPrefix) {

        this.routingStrategy = BaseRoutingDirectionStrategy.forRoutingDirection(direction);
        this.edgeSpacing = edgeSpacing;
        this.conflictThreshold = CONFLICT_THRESHOLD_FACTOR * edgeSpacing;
        this.debugPrefix = debugPrefix;
    }


    ///////////////////////////////////////////////////////////////////////////////
    // Edge Routing

    /**
     * Route edges between the given layers.
     *
     * @param monitor
     *            the progress monitor we're using.
     * @param layeredGraph
     *            the layered graph.
     * @param sourceLayerNodes
     *            the left layer. May be {@code null}.
     * @param sourceLayerIndex
     *            the source layer's index. Ignored if there is no source layer.
     * @param targetLayerNodes
     *            the right layer. May be {@code null}.
     * @param startPos
     *            horizontal position of the first routing slot
     * @return the number of routing slots for this layer
     */
    public int routeEdges(final IElkProgressMonitor monitor, final LGraph layeredGraph,
            final Iterable<LNode> sourceLayerNodes, final int sourceLayerIndex, final Iterable<LNode> targetLayerNodes,
            final double startPos) {

        // Keep track of our hyperedge segements, and which ports they were created for
        Map<LPort, HyperEdgeSegment> portToEdgeSegmentMap = Maps.newHashMap();
        List<HyperEdgeSegment> edgeSegments = Lists.newArrayList();

        // create hyperedge segments for eastern output ports of the left layer and for western output ports of the
        // right layer
        createHyperEdgeSegments(
                sourceLayerNodes, routingStrategy.getSourcePortSide(), edgeSegments, portToEdgeSegmentMap);
        createHyperEdgeSegments(
                targetLayerNodes, routingStrategy.getTargetPortSide(), edgeSegments, portToEdgeSegmentMap);

        // Our critical conflict threshold is a fraction of the minimum distance between two horizontal hyperedge
        // segments
        criticalConflictThreshold = CRITICAL_CONFLICT_THRESHOLD_FACTOR * minimumHorizontalSegmentDistance(edgeSegments);
        
        // create dependencies for the hyperedge segment ordering graph and note how many critical dependencies have
        // been created
        int criticalDependencyCount = 0;
        for (int firstIdx = 0; firstIdx < edgeSegments.size() - 1; firstIdx++) {
            HyperEdgeSegment firstSegment = edgeSegments.get(firstIdx);
            for (int secondIdx = firstIdx + 1; secondIdx < edgeSegments.size(); secondIdx++) {
                criticalDependencyCount += createDependencyIfNecessary(firstSegment, edgeSegments.get(secondIdx));
            }
        }

        // write the full dependency graph to an output file
        // elkjs-exclude-start
        if (debugPrefix != null && monitor.isLoggingEnabled()) {
            DebugUtil.logDebugGraph(
                    monitor,
                    layeredGraph,
                    sourceLayerNodes == null ? 0 : sourceLayerIndex + 1,
                    edgeSegments,
                    debugPrefix,
                    sourceLayerIndex + "-full");
        }
        // elkjs-exclude-end
        
        // if there are at least two critical dependencies, there may be critical cycles that need to be broken
        Random random = layeredGraph.getProperty(InternalProperties.RANDOM);
        if (criticalDependencyCount >= 2) {
            breakCriticalCycles(edgeSegments, random);
        }

        // break non-critical cycles
        breakNonCriticalCycles(edgeSegments, random);

        // write the acyclic dependency graph to an output file
        // elkjs-exclude-start
        if (debugPrefix != null && monitor.isLoggingEnabled()) {
            DebugUtil.logDebugGraph(
                    monitor,
                    layeredGraph,
                    sourceLayerNodes == null ? 0 : sourceLayerIndex + 1,
                    edgeSegments,
                    debugPrefix,
                    sourceLayerIndex + "-acyclic");
        }
        // elkjs-exclude-end
        
        // assign ranks to the edge segments
        topologicalNumbering(edgeSegments);

        // set bend points with appropriate coordinates
        int rankCount = -1;
        for (HyperEdgeSegment node : edgeSegments) {
            // edges that are just straight lines don't take up a slot and don't need bend points
            if (Math.abs(node.getStartCoordinate() - node.getEndCoordinate()) < TOLERANCE) {
                continue;
            }

            rankCount = Math.max(rankCount, node.getRoutingSlot());

            routingStrategy.calculateBendPoints(node, startPos, edgeSpacing);
        }

        // release the created resources
        routingStrategy.clearCreatedJunctionPoints();
        return rankCount + 1;
    }

    /**
     * Computes and returns the minimum distance between any two adjacent source connections and any two adjacent target
     * connections.
     */
    private double minimumHorizontalSegmentDistance(final List<HyperEdgeSegment> edgeSegments) {
        double minIncomingDistance = minimumDifference(
                edgeSegments.stream().flatMap(segment -> segment.getIncomingConnectionCoordinates().stream()));
        double minOutgoingDistance = minimumDifference(
                edgeSegments.stream().flatMap(segment -> segment.getOutgoingConnectionCoordinates().stream()));
        
        return Math.min(minIncomingDistance, minOutgoingDistance);
    }
    
    /**
     * Returns the smallest difference between any two numbers in the given stream of numbers. If there are less than
     * two numbers, we return {@link Double#MAX_VALUE}.
     */
    private double minimumDifference(final Stream<Double> numberStream) {
        List<Double> numbers = numberStream.sorted().distinct().collect(Collectors.toList());
        
        double minDifference = Double.MAX_VALUE;
        
        if (numbers.size() >= 2) {
            Iterator<Double> iter = numbers.iterator();
            Double currentNumber = iter.next();
            
            while (iter.hasNext()) {
                Double previousNumber = currentNumber;
                currentNumber = iter.next();
                
                // This relies on the fact that the numbers are distinct and sorted ascendingly
                minDifference = Math.min(minDifference, currentNumber.doubleValue() - previousNumber.doubleValue());
            }
        }
        
        return minDifference;
    }


    ///////////////////////////////////////////////////////////////////////////////
    // Hyper Edge Graph Creation

    /**
     * Creates hyperedge segments for the given layer.
     *
     * @param nodes
     *            the layer. May be {@code null}, in which case nothing happens.
     * @param portSide
     *            side of the output ports for whose outgoing edges hyperedge segments should be created.
     * @param hyperEdges
     *            list the created hyperedge segments should be added to.
     * @param portToHyperEdgeSegmentMap
     *            map from ports to hyperedge segments that should be filled.
     */
    private void createHyperEdgeSegments(final Iterable<LNode> nodes, final PortSide portSide,
            final List<HyperEdgeSegment> hyperEdges, final Map<LPort, HyperEdgeSegment> portToHyperEdgeSegmentMap) {

        if (nodes != null) {
            for (LNode node : nodes) {
                for (LPort port : node.getPorts(PortType.OUTPUT, portSide)) {
                    HyperEdgeSegment hyperEdge = portToHyperEdgeSegmentMap.get(port);
                    if (hyperEdge == null) {
                        hyperEdge = new HyperEdgeSegment(routingStrategy);
                        hyperEdges.add(hyperEdge);
                        hyperEdge.addPortPositions(port, portToHyperEdgeSegmentMap);
                    }
                }
            }
        }
    }

    /**
     * Create dependencies between the two given hyperedge segments, if one is needed. This method is used not just
     * here, but also by {@link HyperEdgeSegmentSplitter}.
     *
     * @param he1
     *            first hyperedge segments
     * @param he2
     *            second hyperedge segments
     * @return the number of critical dependencies that were added
     */
    int createDependencyIfNecessary(final HyperEdgeSegment he1, final HyperEdgeSegment he2) {
        // check if at least one of the two nodes is just a straight line; those don't
        // create dependencies since they don't take up a slot
        if (Math.abs(he1.getStartCoordinate() - he1.getEndCoordinate()) < TOLERANCE
                || Math.abs(he2.getStartCoordinate() - he2.getEndCoordinate()) < TOLERANCE) {
            
            return 0;
        }

        // compare number of conflicts for both variants
        int conflicts1 = countConflicts(he1.getOutgoingConnectionCoordinates(), he2.getIncomingConnectionCoordinates());
        int conflicts2 = countConflicts(he2.getOutgoingConnectionCoordinates(), he1.getIncomingConnectionCoordinates());
        
        boolean criticalConflictsDetected =
                conflicts1 == CRITICAL_CONFLICTS_DETECTED || conflicts2 == CRITICAL_CONFLICTS_DETECTED;
        int criticalDependencyCount = 0;
        
        if (criticalConflictsDetected) {
            // Check which critical dependencies have to be added
            if (conflicts1 == CRITICAL_CONFLICTS_DETECTED) {
                // hyperedge 1 MUST NOT be left of hyperedge 2, since that would cause critical conflicts
                HyperEdgeSegmentDependency.createAndAddCritical(he2, he1);
                criticalDependencyCount++;
            }
            
            if (conflicts2 == CRITICAL_CONFLICTS_DETECTED) {
                // hyperedge 2 MUST NOT be left of hyperedge 1, since that would cause critical conflicts
                HyperEdgeSegmentDependency.createAndAddCritical(he1, he2);
                criticalDependencyCount++;
            }
            
        } else {
            // we did not detect critical conflicts, so count the number of crossings for both variants
            int crossings1 = countCrossings(
                    he1.getOutgoingConnectionCoordinates(), he2.getStartCoordinate(), he2.getEndCoordinate());
            crossings1 += countCrossings(
                    he2.getIncomingConnectionCoordinates(), he1.getStartCoordinate(), he1.getEndCoordinate());
            int crossings2 = countCrossings(
                    he2.getOutgoingConnectionCoordinates(), he1.getStartCoordinate(), he1.getEndCoordinate());
            crossings2 += countCrossings(
                    he1.getIncomingConnectionCoordinates(), he2.getStartCoordinate(), he2.getEndCoordinate());
            
            // compute the penalty; crossings are deemed worse than (non-critical) conflicts
            int depValue1 = CONFLICT_PENALTY * conflicts1 + CROSSING_PENALTY * crossings1;
            int depValue2 = CONFLICT_PENALTY * conflicts2 + CROSSING_PENALTY * crossings2;
            
            if (depValue1 < depValue2) {
                // hyperedge 1 wants to be left of hyperedge 2
                HyperEdgeSegmentDependency.createAndAddRegular(he1, he2, depValue2 - depValue1);
            } else if (depValue1 > depValue2) {
                // hyperedge 2 wants to be left of hyperedge 1
                HyperEdgeSegmentDependency.createAndAddRegular(he2, he1, depValue1 - depValue2);
            } else if (depValue1 > 0 && depValue2 > 0) {
                // create two dependencies with zero weight
                HyperEdgeSegmentDependency.createAndAddRegular(he1, he2, 0);
                HyperEdgeSegmentDependency.createAndAddRegular(he2, he1, 0);
            }
        }
        
        return criticalDependencyCount;
    }

    /**
     * Counts the number of conflicts for the given lists of positions.
     *
     * @param posis1
     *            sorted list of positions
     * @param posis2
     *            sorted list of positions
     * @return number of positions that overlap, or {@link #CRITICAL_CONFLICTS_DETECTED} if a critical conflict was
     *         detected.
     */
    private int countConflicts(final List<Double> posis1, final List<Double> posis2) {
        int conflicts = 0;

        if (!posis1.isEmpty() && !posis2.isEmpty()) {
            Iterator<Double> iter1 = posis1.iterator();
            Iterator<Double> iter2 = posis2.iterator();
            double pos1 = iter1.next();
            double pos2 = iter2.next();
            boolean hasMore = true;

            do {
                if (pos1 > pos2 - criticalConflictThreshold && pos1 < pos2 + criticalConflictThreshold) {
                    // We're done as soon as we find a single critical conflict
                    return -1;
                } else if (pos1 > pos2 - conflictThreshold && pos1 < pos2 + conflictThreshold) {
                    conflicts++;
                }

                if (pos1 <= pos2 && iter1.hasNext()) {
                    pos1 = iter1.next();
                } else if (pos2 <= pos1 && iter2.hasNext()) {
                    pos2 = iter2.next();
                } else {
                    hasMore = false;
                }
            } while (hasMore);
        }

        return conflicts;
    }

    /**
     * Counts the number of crossings for a given list of positions.
     *
     * @param posis sorted list of positions
     * @param start start of the critical area
     * @param end end of the critical area
     * @return number of positions in the critical area
     */
    static int countCrossings(final List<Double> posis, final double start, final double end) {
        int crossings = 0;
        for (double pos : posis) {
            if (pos > end) {
                break;
            } else if (pos >= start) {
                crossings++;
            }
        }
        return crossings;
    }

    ///////////////////////////////////////////////////////////////////////////////
    // Cycle Breaking

    /**
     * Finds and breaks critical cycles to avoid edge overlaps. The critical cycles we find will be handled by splitting
     * an edge segment into two.
     */
    private void breakCriticalCycles(final List<HyperEdgeSegment> edgeSegments, final Random random) {
        List<HyperEdgeSegmentDependency> cycleDependencies =
                HyperEdgeCycleDetector.detectCycles(edgeSegments, true, random);
        
        // Lazy initialisation
        if (segmentSplitter == null) {
            segmentSplitter = new HyperEdgeSegmentSplitter(this);
        }
        
        segmentSplitter.splitSegments(cycleDependencies, edgeSegments, criticalConflictThreshold);
    }

    /**
     * Finds and breaks non-critical cycles by removing and reversing non-critical dependencies. This method is used by
     * the self loop routing code as well.
     */
    public static void breakNonCriticalCycles(final List<HyperEdgeSegment> edgeSegments, final Random random) {
        List<HyperEdgeSegmentDependency> cycleDependencies =
                HyperEdgeCycleDetector.detectCycles(edgeSegments, false, random);
        
        for (HyperEdgeSegmentDependency cycleDependency : cycleDependencies) {
            if (cycleDependency.getWeight() == 0) {
                // Simply remove this dependency. This assumes that only two-cycles will have dependency weight 0
                cycleDependency.remove();
            } else {
                cycleDependency.reverse();
            }
        }
    }

    ///////////////////////////////////////////////////////////////////////////////
    // Topological Ordering

    /**
     * Perform a topological numbering of the given hyperedge segments.
     *
     * @param segments list of hyperedge segments
     */
    private static void topologicalNumbering(final List<HyperEdgeSegment> segments) {
        // determine sources, targets, incoming count and outgoing count; targets are only
        // added to the list if they only connect westward ports (that is, if all their
        // horizontal segments point to the right)
        List<HyperEdgeSegment> sources = Lists.newArrayList();
        List<HyperEdgeSegment> rightwardTargets = Lists.newArrayList();
        for (HyperEdgeSegment node : segments) {
            node.setInWeight(node.getIncomingSegmentDependencies().size());
            node.setOutWeight(node.getOutgoingSegmentDependencies().size());

            if (node.getInWeight() == 0) {
                sources.add(node);
            }

            if (node.getOutWeight() == 0 && node.getIncomingConnectionCoordinates().size() == 0) {
                rightwardTargets.add(node);
            }
        }

        int maxRank = -1;

        // assign ranks using topological numbering
        while (!sources.isEmpty()) {
            HyperEdgeSegment node = sources.remove(0);
            for (HyperEdgeSegmentDependency dep : node.getOutgoingSegmentDependencies()) {
                HyperEdgeSegment target = dep.getTarget();
                target.setRoutingSlot(Math.max(target.getRoutingSlot(), node.getRoutingSlot() + 1));
                maxRank = Math.max(maxRank, target.getRoutingSlot());

                target.setInWeight(target.getInWeight() - 1);
                if (target.getInWeight() == 0) {
                    sources.add(target);
                }
            }
        }

        /* If we stopped here, hyperedge segments that don't have any horizontal segments pointing
         * leftward would be ranked just like every other hyperedge segments. This would move back
         * edges too far away from their target node. To remedy that, we move all hyperedge segments
         * with horizontal segments only pointing rightwards as far right as possible.
         */
        if (maxRank > -1) {
            // assign all target nodes with horizontal segments pointing to the right the
            // rightmost rank
            for (HyperEdgeSegment node : rightwardTargets) {
                node.setRoutingSlot(maxRank);
            }

            // let all other segments with horizontal segments pointing rightwards move as
            // far right as possible
            while (!rightwardTargets.isEmpty()) {
                HyperEdgeSegment node = rightwardTargets.remove(0);

                // The node only has connections to western ports
                for (HyperEdgeSegmentDependency dep : node.getIncomingSegmentDependencies()) {
                    HyperEdgeSegment source = dep.getSource();
                    if (source.getIncomingConnectionCoordinates().size() > 0) {
                        continue;
                    }

                    source.setRoutingSlot(Math.min(source.getRoutingSlot(), node.getRoutingSlot() - 1));

                    source.setOutWeight(source.getOutWeight() - 1);
                    if (source.getOutWeight() == 0) {
                        rightwardTargets.add(source);
                    }
                }
            }
        }
    }

}
