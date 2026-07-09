/*******************************************************************************
 * Copyright (c) 2020 Kiel University and others.
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

import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Stream;

import org.eclipse.elk.core.util.Pair;

import com.google.common.collect.Streams;

/**
 * Responsible for splitting {@link HyperEdgeSegment}s in order to avoid overlaps. This class assumes that critical
 * dependency cycles have been found and that a set of critical dependencies whose removal will break those cycles has
 * been computed. These dependencies will determine the {@link HyperEdgeSegment}s that will be split.
 * 
 * <p>
 * All throughout this code, the terminology will refer to the horizontal layout case.
 * </p>
 */
public final class HyperEdgeSegmentSplitter {
    
    /** The routing generator that created us. */
    private OrthogonalRoutingGenerator routingGenerator;
    
    /**
     * Creates a new instance to be used by the given routing generator.
     */
    public HyperEdgeSegmentSplitter(final OrthogonalRoutingGenerator routingGenerator) {
        this.routingGenerator = routingGenerator;
    }
    
    /**
     * Breaks critical dependency cycles by resolving the given dependencies. For each dependency, we split one of the
     * involved {@link HyperEdgeSegment}s, which will result in modified dependencies and additional segments, which are
     * added to the given list of segments.
     */
    public void splitSegments(final List<HyperEdgeSegmentDependency> dependenciesToResolve,
            final List<HyperEdgeSegment> segments, final double criticalConflictThreshold) {
        
        // Only start preparations if there's going to be things to do
        if (dependenciesToResolve.isEmpty()) {
            return;
        }
        
        // Collect all relevant spaces between horizontal segments that are large enough to house another horizontal
        // segment without causing additional conflicts
        List<FreeArea> freeAreas = findFreeAreas(segments, criticalConflictThreshold);
        
        // For each dependency, choose which segment to split
        Set<HyperEdgeSegment> segmentsToSplit = decideWhichSegmentsToSplit(dependenciesToResolve);
        
        // Split the segments in order from smalles to largest. The smallest ones need to be split first since they
        // have fewer options for where to put their horizontal connecting segments.
        segmentsToSplit.stream()
            .sorted((hes1, hes2) -> Double.compare(hes1.getLength(), hes2.getLength()))
            .forEach(segment -> split(segment, segments, freeAreas, criticalConflictThreshold));
    }
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Finding Space
    
    /**
     * Assembles the areas between horizontal segments that are large enough to allow another horizontal segment to slip
     * in.
     */
    private List<FreeArea> findFreeAreas(final List<HyperEdgeSegment> segments,
            final double criticalConflictThreshold) {
        
        List<FreeArea> freeAreas = new ArrayList<>();
        
        // Retrieve all positions where hyperedge segments connect to ports, and sort them
        Stream<Double> inCoordinates = segments.stream().flatMap(s -> s.getIncomingConnectionCoordinates().stream());
        Stream<Double> outCoordinates = segments.stream().flatMap(s -> s.getOutgoingConnectionCoordinates().stream());
        
        double[] sortedCoordinates = Streams.concat(inCoordinates, outCoordinates)
            .mapToDouble(coordinate -> coordinate.doubleValue())
            .sorted()
            .toArray();
        
        // Go through each pair of coordinates and create free areas for those that are at least twice the critical
        // threshold
        for (int i = 1; i < sortedCoordinates.length; i++) {
            if (sortedCoordinates[i] - sortedCoordinates[i - 1] >= 2 * criticalConflictThreshold) {
                freeAreas.add(new FreeArea(
                        sortedCoordinates[i - 1] + criticalConflictThreshold,
                        sortedCoordinates[i] - criticalConflictThreshold));
            }
        }
        
        return freeAreas;
    }
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Split Segment Decisions

    /**
     * Given a list of hyperedge dependencies, this method assembles a list of {@link HyperEdgeSegment}s which, if
     * split, will cause the critical cycles caused by the dependencies to be broken.
     */
    private LinkedHashSet<HyperEdgeSegment> decideWhichSegmentsToSplit(
            final List<HyperEdgeSegmentDependency> dependencies) {
        
        // We want to be sure to have a stable iteration order, but we also want fast containment checks
        LinkedHashSet<HyperEdgeSegment> segmentsToSplit = new LinkedHashSet<>();
        
        for (HyperEdgeSegmentDependency dependency : dependencies) {
            HyperEdgeSegment sourceSegment = dependency.getSource();
            HyperEdgeSegment targetSegment = dependency.getTarget();
            
            // If either of the involved segments were already selected for splitting because of another dependency,
            // that's sufficient
            if (segmentsToSplit.contains(sourceSegment) || segmentsToSplit.contains(targetSegment)) {
                continue;
            }
            
            // One segment will be split, and the other one will be remembered to have caused the split. The latter
            // will have to stay between the two segments the former will be split into, or else we'll have overlaps
            // again
            HyperEdgeSegment segmentToSplit = sourceSegment;
            HyperEdgeSegment segmentCausingSplit = targetSegment;
            
            // This is perfectly fine unless the segment to split represents a hyperedge and the other does not. We
            // prefer splitting regular edges since hyperedges have a higher chance of causing additional crossings.
            // Note: the original code took the layout direction into account here. We may want to investigate at some
            // point whether that's necessary.
            if (sourceSegment.representsHyperedge() && !targetSegment.representsHyperedge()) {
                segmentToSplit = targetSegment;
                segmentCausingSplit = sourceSegment;
            }
            
            segmentsToSplit.add(segmentToSplit);
            segmentToSplit.setSplitBy(segmentCausingSplit);
        }
        
        return segmentsToSplit;
    }
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Actual Splitting
    
    /**
     * Splits the given segment at the optimal position.
     */
    private void split(final HyperEdgeSegment segment, final List<HyperEdgeSegment> segments,
            final List<FreeArea> freeAreas, final double criticalConflictThreshold) {
        
        // Split the segment at the best position and add the new segment to our list
        double splitPosition = computePositionToSplitAndUpdateFreeAreas(segment, freeAreas, criticalConflictThreshold);
        segments.add(segment.splitAt(splitPosition));
        
        // Update the dependencies to reflect the new situation
        updateDependencies(segment, segments);
    }

    /**
     * Adds all necessary dependencies for the given segment and its split partner. The segment was recently split, and
     * neither it nor its split partner have any dependencies beyond the critical ones.
     */
    private void updateDependencies(final HyperEdgeSegment segment, final List<HyperEdgeSegment> segments) {
        HyperEdgeSegment splitCausingSegment = segment.getSplitBy();
        HyperEdgeSegment splitPartner = segment.getSplitPartner();
        
        // The segment currently has no dependencies at all. We first need for the segments to be ordered like this:
        //    segment ---> split-causing segment ---> split partner
        HyperEdgeSegmentDependency.createAndAddCritical(segment, splitCausingSegment);
        HyperEdgeSegmentDependency.createAndAddCritical(splitCausingSegment, splitPartner);

        // Now we just need to re-introduce dependencies to other segments
        for (HyperEdgeSegment otherSegment : segments) {
            // We already have our dependencies between our three segments involved in the conflict settled
            if (otherSegment != splitCausingSegment && otherSegment != segment && otherSegment != splitPartner) {
                routingGenerator.createDependencyIfNecessary(otherSegment, segment);
                routingGenerator.createDependencyIfNecessary(otherSegment, splitPartner);
            }
        }
    }
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Split Position Computation

    /**
     * Returns a hopefully good position to split the given segment at, and updates the list of free areas since an area
     * now either disappears or is split into two smaller areas.
     */
    private double computePositionToSplitAndUpdateFreeAreas(final HyperEdgeSegment segment,
            final List<FreeArea> freeAreas, final double criticalConflictThreshold) {
        
        // We shall iterate over the available areas (which are sorted by ascending positions) and find the index of the
        // first and the last area we could use since they are in our segment's reach
        int firstPossibleAreaIndex = -1;
        int lastPossibleAreaIndex = -1;
        
        for (int i = 0; i < freeAreas.size(); i++) {
            FreeArea currArea = freeAreas.get(i);
            
            if (currArea.startPosition > segment.getEndCoordinate()) {
                // We're past the possible areas, so stop
                break;
            } else if (currArea.endPosition >= segment.getStartCoordinate()) {
                // We've found a possible area; it might be the first
                if (firstPossibleAreaIndex < 0) {
                    firstPossibleAreaIndex = i;
                }
                
                lastPossibleAreaIndex = i;
            }
        }
        
        // Determine the position to split the segment at
        double splitPosition = center(segment);

        if (firstPossibleAreaIndex >= 0) {
            // There are areas we can use
            int bestAreaIndex = chooseBestAreaIndex(segment, freeAreas, firstPossibleAreaIndex, lastPossibleAreaIndex);
            
            // We'll use the best area's centre and update the area list
            splitPosition = center(freeAreas.get(bestAreaIndex));
            useArea(freeAreas, bestAreaIndex, criticalConflictThreshold);
        }
        
        return splitPosition;
    }
    
    /**
     * Finds the best area to use for splitting the given segment between the given two area indices (both inclusive).
     */
    private int chooseBestAreaIndex(final HyperEdgeSegment segment, final List<FreeArea> freeAreas,
            final int fromIndex, final int toIndex) {
        
        int bestAreaIndex = fromIndex;
        
        if (fromIndex < toIndex) {
            // We have more areas to choose from, so rate them and find the best one. We need to simulate splitting the
            // segment so that we can count potential crossings
            Pair<HyperEdgeSegment, HyperEdgeSegment> splitSegments = segment.simulateSplit();
            HyperEdgeSegment splitSegment = splitSegments.getFirst();
            HyperEdgeSegment splitPartner = splitSegments.getSecond();
            
            FreeArea bestArea = freeAreas.get(bestAreaIndex);
            AreaRating bestRating = rateArea(segment, splitSegment, splitPartner, bestArea);
            
            for (int i = fromIndex + 1; i <= toIndex; i++) {
                // Determine how good our current area is
                FreeArea currArea = freeAreas.get(i);
                AreaRating currRating = rateArea(segment, splitSegment, splitPartner, currArea);
                
                if (isBetter(currArea, currRating, bestArea, bestRating)) {
                    bestArea = currArea;
                    bestRating = currRating;
                    bestAreaIndex = i;
                }
            }
        }
        
        return bestAreaIndex;
    }

    /**
     * Rates what would happen if the given split segments were connected through the given area.
     */
    private AreaRating rateArea(final HyperEdgeSegment segment, final HyperEdgeSegment splitSegment,
            final HyperEdgeSegment splitPartner, final FreeArea area) {
        
        // The area's centre would be used to link the two split segments, so we need to add that to their incident
        // connections
        double areaCentre = center(area);

        splitSegment.getOutgoingConnectionCoordinates().clear();
        splitSegment.getOutgoingConnectionCoordinates().add(areaCentre);
        
        splitPartner.getIncomingConnectionCoordinates().clear();
        splitPartner.getIncomingConnectionCoordinates().add(areaCentre);
        
        // We need to count the dependencies and crossings that the split partners would cause with the original
        // segments incident dependencies
        AreaRating rating = new AreaRating(0, 0);
        
        for (HyperEdgeSegmentDependency dependency : segment.getIncomingSegmentDependencies()) {
            HyperEdgeSegment otherSegment = dependency.getSource();
            
            updateConsideringBothOrderings(rating, splitSegment, otherSegment);
            updateConsideringBothOrderings(rating, splitPartner, otherSegment);
        }
        
        for (HyperEdgeSegmentDependency dependency : segment.getOutgoingSegmentDependencies()) {
            HyperEdgeSegment otherSegment = dependency.getTarget();
            
            updateConsideringBothOrderings(rating, splitSegment, otherSegment);
            updateConsideringBothOrderings(rating, splitPartner, otherSegment);
        }
        
        // There will be two additional dependencies: splitSegment --> splitBySegment --> splitPartner. The order
        // between the three will not change, so we only have to count crossings for this concrete order
        rating.dependencies += 2;
        
        rating.crossings += countCrossingsForSingleOrdering(splitSegment, segment.getSplitBy());
        rating.crossings += countCrossingsForSingleOrdering(segment.getSplitBy(), splitPartner);
        
        return rating;
    }
    
    /**
     * Considers both orderings (s1 left of s2, and vice versa) and considers the number of crossings that would ensue.
     * Based on the crossings, determines the minimum number of dependencies and crossings we'd have to expect.
     */
    private void updateConsideringBothOrderings(final AreaRating rating, final HyperEdgeSegment s1,
            final HyperEdgeSegment s2) {
        
        int crossingsS1LeftOfS2 = countCrossingsForSingleOrdering(s1, s2);
        int crossingsS2LeftOfS1 = countCrossingsForSingleOrdering(s2, s1);
        
        if (crossingsS1LeftOfS2 == crossingsS2LeftOfS1) {
            // If the crossings are the same, we're only interested if there are more than 0
            if (crossingsS1LeftOfS2 > 0) {
                // Both orderings generate the same number of crossings (more than 0), so we have a two-cycle
                rating.dependencies += 2;
                rating.crossings += crossingsS1LeftOfS2;
            }
        } else {
            // One order is better than the other, so there will be a single dependency
            rating.dependencies += 1;
            rating.crossings += Math.min(crossingsS1LeftOfS2, crossingsS2LeftOfS1);
        }
    }
    
    /**
     * Counts the number of crossings that would ensue between a left and a right segment.
     */
    private int countCrossingsForSingleOrdering(final HyperEdgeSegment left, final HyperEdgeSegment right) {
        return OrthogonalRoutingGenerator.countCrossings(
                left.getOutgoingConnectionCoordinates(), right.getStartCoordinate(), right.getEndCoordinate())
            + OrthogonalRoutingGenerator.countCrossings(
                right.getIncomingConnectionCoordinates(), left.getStartCoordinate(), left.getEndCoordinate());
    }
    
    /**
     * Determines whether the given area with the given rating is better than another.
     */
    private boolean isBetter(final FreeArea currArea, final AreaRating currRating, final FreeArea bestArea,
            final AreaRating bestRating) {
        
        if (currRating.crossings < bestRating.crossings) {
            // First criterion: number of crossings
            return true;
            
        } else if (currRating.crossings == bestRating.crossings) {
            if (currRating.dependencies < bestRating.dependencies) {
                // Second criterion: number of dependencies
                return true;
            
            } else if (currRating.dependencies == bestRating.dependencies) {
                if (currArea.size > bestArea.size) {
                    // Third criterion: size
                    return true;
                }
            }
        }
        
        return false;
    }

    /**
     * When the given area is being used, it falls into two parts which may be usable themselves. This method does the
     * necessary bookkeeping.
     */
    private void useArea(final List<FreeArea> freeAreas, final int usedAreaIndex,
            final double criticalConflictThreshold) {
        
        FreeArea oldArea = freeAreas.get(usedAreaIndex);
        freeAreas.remove(usedAreaIndex);
        
        if (oldArea.size / 2 >= criticalConflictThreshold) {
            // We will probably insert new areas. Keep track of where to insert them
            int insertIndex = usedAreaIndex;
            
            // This area is large enough to split and still have enough space to the position we're now using for our
            // latest segment
            double oldAreaCentre = center(oldArea);

            // Create the two new areas (and be doubly sure that double precision does not bite us)
            double newEnd1 = oldAreaCentre - criticalConflictThreshold;
            if (oldArea.startPosition <= oldAreaCentre - criticalConflictThreshold) {
                FreeArea newArea1 = new FreeArea(oldArea.startPosition, newEnd1);
                freeAreas.add(insertIndex++, newArea1);
            }
            
            double newStart2 = oldAreaCentre + criticalConflictThreshold;
            if (newStart2 <= oldArea.endPosition) {
                FreeArea newArea2 = new FreeArea(newStart2, oldArea.endPosition);
                freeAreas.add(insertIndex, newArea2);
            }
        }
    }

    private static double center(final HyperEdgeSegment s) {
        return center(s.getStartCoordinate(), s.getEndCoordinate());
    }

    private static double center(final FreeArea a) {
        return center(a.startPosition, a.endPosition);
    }

    private static double center(final double p1, final double p2) {
        return (p1 + p2) / 2;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Data Holders
    
    /**
     * Represents a free area between two horizontal edge segments, which can be used to route new horizontal segments
     * that connect the two parts of a split {@link HyperEdgeSegment}.
     */
    private static final class FreeArea {
        
        /** Start of the free area. */
        private final double startPosition;
        /** End of the free area. */
        private final double endPosition;
        /** Convenience variable that stores the size of the area. */
        private final double size;
        
        private FreeArea(final double startPosition, final double endPosition) {
            assert endPosition >= startPosition;
            
            this.startPosition = startPosition;
            this.endPosition = endPosition;
            this.size = endPosition - startPosition;
        }
        
    }
    
    /**
     * An indication of what would happen if a given segment was connected to its split partner through this area.
     */
    private static final class AreaRating {
        
        /** The number of dependencies caused by using this area. */
        private int dependencies;
        /** The number of crossings using this area would cause. */
        private int crossings;
        
        private AreaRating(final int dependencies, final int crossings) {
            this.dependencies = dependencies;
            this.crossings = crossings;
        }
        
    }
    
}
