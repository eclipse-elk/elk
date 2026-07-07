/*******************************************************************************
 * Copyright (c) 2010, 2020 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.layered.p5edges.orthogonal;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.Random;
import java.util.Set;

import org.eclipse.elk.alg.layered.p5edges.orthogonal.HyperEdgeSegmentDependency.DependencyType;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;

/**
 * Finds a set of dependencies to remove or reverse to break cycles in the conflict graph of {@link HyperEdgeSegment}s
 * and {@link HyperEdgeSegmentDependency}s. What clients will do with these, exactly, will be up to them. The hope is to
 * return the minimum-weight set of dependencies such that their reversal or removal causes the graph to become acyclic.
 * Inspired by:
 * <ul>
 *   <li>Eades, Lin, Smyth. A fast and effective heuristic for the feedback arc set problem. In <i>Information
 *     Processing Letters</i>, 1993.</li>
 * </ul>
 */
public final class HyperEdgeCycleDetector {
    
    /**
     * No instance required.
     */
    private HyperEdgeCycleDetector() {
        assert false;
    }

    /**
     * Finds a set of dependencies whose reversal or removal will make the graph acyclic. The method can concentrate on
     * critical dependencies only, or detect cycles among all dependencies. The latter assumes that there are no cycles
     * among critical dependencies left and will never return a critical dependency in the result. (doing so would cause
     * clients to reverse or remove critical dependencies, which is a big no-no, to put it eloquently!) This code also
     * assumes that only dependencies involved in non-critical two-cycles will have a weight of zero.
     *
     * @param segments
     *            list of hyper edge segments.
     * param criticalOnly
     *            {@code true} if we should only detect cycles among critical dependencies.
     * @param random
     *            random number generator.
     * @return list of dependencies whose removal or reversal will make the graph acyclic.
     */
    public static List<HyperEdgeSegmentDependency> detectCycles(final List<HyperEdgeSegment> segments,
            final boolean criticalOnly, final Random random) {
        
        List<HyperEdgeSegmentDependency> result = new ArrayList<>();
        
        LinkedList<HyperEdgeSegment> sources = Lists.newLinkedList();
        LinkedList<HyperEdgeSegment> sinks = Lists.newLinkedList();

        // initialize values for the algorithm
        initialize(segments, sources, sinks, criticalOnly);

        // assign marks to all nodes
        computeLinearOrderingMarks(segments, sources, sinks, criticalOnly, random);

        // process edges that point left: remove those of zero weight, reverse the others
        for (HyperEdgeSegment source : segments) {
            for (HyperEdgeSegmentDependency outDependency : source.getOutgoingSegmentDependencies()) {
                // Only consider critical dependencies, if required
                if (!criticalOnly || outDependency.getType() == DependencyType.CRITICAL) {
                    if (source.mark > outDependency.getTarget().mark) {
                        result.add(outDependency);
                    }
                }
            }
        }
        
        return result;
    }

    /**
     * Initializes the mark, in weight and out weight of each hyper edge segment. Also adds all sources (segments
     * without incoming weight 0) and sinks (segments with outgoing weight 0) to their respective lists. Once this is
     * complete, all segments are marked from {@code -1} to {@code -segments.size()}.
     */
    private static void initialize(final List<HyperEdgeSegment> segments, final List<HyperEdgeSegment> sources,
            final List<HyperEdgeSegment> sinks, final boolean criticalOnly) {
        
        int nextMark = -1;
        for (HyperEdgeSegment segment : segments) {
            segment.mark = nextMark--;

            // Sum up the weights of our critical dependencies
            int criticalInWeight = segment.getIncomingSegmentDependencies().stream()
                    .filter(dep -> dep.getType() == DependencyType.CRITICAL)
                    .mapToInt(dep -> dep.getWeight())
                    .sum();

            int criticalOutWeight = segment.getOutgoingSegmentDependencies().stream()
                    .filter(dep -> dep.getType() == DependencyType.CRITICAL)
                    .mapToInt(dep -> dep.getWeight())
                    .sum();

            // If we're only considering critical dependencies, we'll ignore the others
            int inWeight = criticalInWeight;
            int outWeight = criticalOutWeight;

            if (!criticalOnly) {
                // Just sum up everything
                inWeight = segment.getIncomingSegmentDependencies().stream()
                        .mapToInt(dep -> dep.getWeight())
                        .sum();

                outWeight = segment.getOutgoingSegmentDependencies().stream()
                        .mapToInt(dep -> dep.getWeight())
                        .sum();
            }

            // Apply the weight
            segment.setInWeight(inWeight);
            segment.setCriticalInWeight(criticalInWeight);
            segment.setOutWeight(outWeight);
            segment.setCriticalOutWeight(criticalOutWeight);

            // Add the segment to either sources or sinks if the corresponding weight is zero
            if (outWeight == 0) {
                sinks.add(segment);
            } else if (inWeight == 0) {
                sources.add(segment);
            }
        }
    }

    /**
     * Computes marks for all segments based on a linear ordering or the segments. Marks will be mutually different, and
     * sinks will have higher marks than sources.
     */
    private static void computeLinearOrderingMarks(final List<HyperEdgeSegment> segments,
            final LinkedList<HyperEdgeSegment> sources, final LinkedList<HyperEdgeSegment> sinks,
            final boolean criticalOnly, final Random random) {
        
        Set<HyperEdgeSegment> unprocessed = Sets.newTreeSet(segments);
        List<HyperEdgeSegment> maxSegments = new ArrayList<HyperEdgeSegment>();
        
        // We'll mark sinks with marks < markBase and sources with marks > markBase. Sink marks will later be offset to
        // be higher than the source marks, but this way, the sink marks will reflect the order in which the sinks were
        // discovered and added.
        int markBase = segments.size();
        int nextSinkMark = markBase - 1;
        int nextSourceMark = markBase + 1;

        while (!unprocessed.isEmpty()) {
            while (!sinks.isEmpty()) {
                HyperEdgeSegment sink = sinks.removeFirst();
                unprocessed.remove(sink);
                sink.mark = nextSinkMark--;
                updateNeighbors(sink, sources, sinks, criticalOnly);
            }

            while (!sources.isEmpty()) {
                HyperEdgeSegment source = sources.removeFirst();
                unprocessed.remove(source);
                source.mark = nextSourceMark++;
                updateNeighbors(source, sources, sinks, criticalOnly);
            }

            // If any segments are still unprocessed, they are neither source nor sink. Assemble the list of segments
            // with the highest out flow (out weight - in weight), to be placed among the sources. If we're looking at
            // both, critical and non-critical dependencies, we must be sure that critical dependencies will always
            // point rightwards; we thus stop immediately once we find one
            int maxOutflow = Integer.MIN_VALUE;
            for (HyperEdgeSegment segment : unprocessed) {
                // If we're not only regarding critical dependencies, we must ensure that critical dependencies will
                // always point to the right to prevent them from being reversed later on. Thus, once we find a single
                // segment that still has an outgoing critical dependency and no incoming ones, we'll take that and
                // leave!
                if (!criticalOnly && segment.getCriticalOutWeight() > 0 && segment.getCriticalInWeight() <= 0) {
                    maxSegments.clear();
                    maxSegments.add(segment);
                    break;
                }
                
                int outflow = segment.getOutWeight() - segment.getInWeight();
                if (outflow >= maxOutflow) {
                    if (outflow > maxOutflow) {
                        maxSegments.clear();
                        maxOutflow = outflow;
                    }
                    maxSegments.add(segment);
                }
            }

            // If there are segments with maximal out flow, select one randomly; this might yield new sources and sinks
            if (!maxSegments.isEmpty()) {
                HyperEdgeSegment maxNode = maxSegments.get(random.nextInt(maxSegments.size()));
                unprocessed.remove(maxNode);
                maxNode.mark = nextSourceMark++;
                updateNeighbors(maxNode, sources, sinks, criticalOnly);
                maxSegments.clear();
            }
        }

        // shift ranks that are left of the mark base so that sinks now have higher marks than sources
        int shiftBase = segments.size() + 1;
        for (HyperEdgeSegment node : segments) {
            if (node.mark < markBase) {
                node.mark = node.mark + shiftBase;
            }
        }
    }

    /**
     * Updates in-weight and out-weight values of the neighbors of the given node, simulating its removal from the
     * graph. The sources and sinks lists are also updated.
     *
     * @param node
     *            node for which neighbors are updated
     * @param sources
     *            list of sources
     * @param sinks
     *            list of sinks
     * @param criticalOnly
     *            {@code true} if only critical dependencies should be taken into account.
     */
    private static void updateNeighbors(final HyperEdgeSegment node, final List<HyperEdgeSegment> sources,
            final List<HyperEdgeSegment> sinks, final boolean criticalOnly) {

        // process following nodes
        for (HyperEdgeSegmentDependency dep : node.getOutgoingSegmentDependencies()) {
            if (!criticalOnly || dep.getType() == DependencyType.CRITICAL) {
                HyperEdgeSegment target = dep.getTarget();
                if (target.mark < 0 && dep.getWeight() > 0) {
                    // Remove weight (and possibly critical weight) from the target
                    target.setInWeight(target.getInWeight() - dep.getWeight());
                    if (dep.getType() == DependencyType.CRITICAL) {
                        target.setCriticalInWeight(target.getCriticalInWeight() - dep.getWeight());
                    }
                    
                    if (target.getInWeight() <= 0 && target.getOutWeight() > 0) {
                        sources.add(target);
                    }
                }
            }
        }

        // process preceding nodes
        for (HyperEdgeSegmentDependency dep : node.getIncomingSegmentDependencies()) {
            if (!criticalOnly || dep.getType() == DependencyType.CRITICAL) {
                HyperEdgeSegment source = dep.getSource();
                if (source.mark < 0 && dep.getWeight() > 0) {
                    // Remove weight (and possibly critical weight) from the source
                    source.setOutWeight(source.getOutWeight() - dep.getWeight());
                    if (dep.getType() == DependencyType.CRITICAL) {
                        source.setCriticalOutWeight(source.getCriticalOutWeight() - dep.getWeight());
                    }
                    
                    if (source.getOutWeight() <= 0 && source.getInWeight() > 0) {
                        sinks.add(source);
                    }
                }
            }
        }
    }
    
}
