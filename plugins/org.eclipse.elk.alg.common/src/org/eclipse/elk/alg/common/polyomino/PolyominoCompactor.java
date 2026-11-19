/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
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
package org.eclipse.elk.alg.common.polyomino;

import java.util.Comparator;
import java.util.List;
import java.util.Set;
import java.util.function.BiFunction;
import java.util.function.Function;
import java.util.stream.Collectors;

import org.eclipse.elk.alg.common.compaction.options.PolyominoOptions;
import org.eclipse.elk.alg.common.polyomino.structures.Direction;
import org.eclipse.elk.alg.common.polyomino.structures.PlanarGrid;
import org.eclipse.elk.alg.common.polyomino.structures.Polyomino;
import org.eclipse.elk.alg.common.polyomino.structures.Polyominoes;
import org.eclipse.elk.alg.common.utils.UniqueTriple;
import org.eclipse.elk.core.util.Pair;

/**
 * <p>
 * Implementation of a packing algorithm based on using polyominoes (i. e. low resolution representations of shapes like
 * nodes, components, graphs, etc.).
 * </p>
 * 
 * <p>
 * The approach is based on:
 * </p>
 * 
 * <p>
 * Freivalds, K., Dogrusoz, U. and Kikusts, P. Disconnected Graph Layout and the Polyomino Packing Approach. in Mutzel,
 * P., Jünger, M. and Leipert, S. (Eds.). Graph Drawing: 9th International Symposium, GD 2001 Vienna, Austria, September
 * 23-26, 2001 Revised Papers. Springer Berlin Heidelberg, Berlin, Heidelberg, 2002, 378-39. DOI: <a
 * href=http://dx.doi.org/10.1007/3-540-45848-4_30>10.1007/3-540-45848-4_30</a>
 * </p>
 * 
 * <p>
 * A particular subset of polyominoes consisting of four cells might remind you of your sheltered childhood, when you
 * still had time to play video games on your popular, barely pocketable, japanese mobile video gaming device (see <a
 * href=https://en.wikipedia.org/wiki/Polyomino>Polyomino (Wikipedia)</a> and <a
 * href=https://en.wikipedia.org/wiki/Tetromino>Tetromino (Wikipedia)</a> for visual reminders, retrieved on
 * 2016-09-20).
 * </p>
 */
public class PolyominoCompactor {

    /**
     * Places {@link Polyomino polyominoes} close together, in order to achieve a minimum area, based on the heuristic
     * algorithm PackPolyominoes from the paper.
     */
    public <P extends Polyomino> void packPolyominoes(final Polyominoes<P> polyHolder) {
        List<P> polys = polyHolder.getPolyominoes();
        PlanarGrid grid = polyHolder.getGrid();

        // 1. Sort polyominoes by some strategy
        switch (polyHolder.getProperty(PolyominoOptions.POLYOMINO_LOW_LEVEL_SORT)) {
        case BY_SIZE:
            polys.sort(new MinPerimeterComparator().reversed());
            break;

        case BY_SIZE_AND_SHAPE:
        default:
            polys.sort(new MinPerimeterComparatorWithShape().reversed());
            break;
        }

        switch (polyHolder.getProperty(PolyominoOptions.POLYOMINO_HIGH_LEVEL_SORT)) {
        case CORNER_CASES_THAN_SINGLE_SIDE_LAST:
            polys.sort(new MinNumOfExtensionsComparator());
            polys.sort(new SingleExtensionSideGreaterThanRestComparator());
            polys.sort(new CornerCasesGreaterThanRestComparator());
            break;
        case NUM_OF_EXTERNAL_SIDES_THAN_NUM_OF_EXTENSIONS_LAST:
        default:
            polys.sort(new MinNumOfExtensionsComparator());
            polys.sort(new MinNumOfExtensionDirectionsComparator());
        }

        // 2. Initialize cost function. Implemented as a BiFunction for future interchangeability with different
        // metrics.
        BiFunction<Pair<Integer, Integer>, Polyomino, Pair<Integer, Integer>> successorBasedOnCost;
        switch (polyHolder.getProperty(PolyominoOptions.POLYOMINO_TRAVERSAL_STRATEGY)) {
        case SPIRAL:
            successorBasedOnCost = new SuccessorMaxNormWindingInMathPosSense();
            break;
        case LINE_BY_LINE:
            successorBasedOnCost = new SuccessorLineByLine();
            break;
        case MANHATTAN:
            successorBasedOnCost = new SuccessorManhattan();
            break;
        case JITTER:
            successorBasedOnCost = new SuccessorJitter();
            break;
        case QUADRANTS_MANHATTAN:
            successorBasedOnCost = new SuccessorQuadrantsGeneric(new SuccessorManhattan());
            break;
        case QUADRANTS_LINE_BY_LINE:
            successorBasedOnCost = new SuccessorQuadrantsGeneric(new SuccessorLineByLine());
            break;
        case COMBINE_LINE_BY_LINE_MANHATTAN:
            successorBasedOnCost = new SuccessorCombination(new SuccessorQuadrantsGeneric(new SuccessorLineByLine()),
                    new SuccessorQuadrantsGeneric(new SuccessorManhattan()));
            break;
        case COMBINE_JITTER_MANHATTAN:
            successorBasedOnCost = new SuccessorCombination(new SuccessorQuadrantsGeneric(new SuccessorJitter()),
                    new SuccessorQuadrantsGeneric(new SuccessorManhattan()));
            break;
        case QUADRANTS_JITTER:
        default:
            successorBasedOnCost = new SuccessorQuadrantsGeneric(new SuccessorJitter());
        }

        for (Polyomino poly : polys) {
            // 3. Start placement of each polyomino at the center of the grid
            int offX = 0;
            int offY = 0;
            Pair<Integer, Integer> next = new Pair<Integer, Integer>(offX, offY);

            // 4. Until no polyomino cell intersects with another polyomino already placed ...
            while (grid.intersectsWithCenterBased(poly, offX, offY)) {
                // ... get next trial position based on cost function
                next = successorBasedOnCost.apply(next, poly);
                offX = next.getFirst();
                offY = next.getSecond();
            }
            // 5. When a valid position is found, mark all new cells in the underlying grid and save the position of the
            // polyomino's center relative to the grid's center
            grid.addFilledCellsFrom(poly, offX, offY);
        }
    }

    ///////////////////////////////////////////////////////////////////////////////
    // Inner Classes

    /**
     * Comparator for sorting {@link Polyomino polyominoes} based on the length of the perimeter of their bounding
     * boxes.
     */
    private class MinPerimeterComparator implements Comparator<PlanarGrid> {
        @Override
        public int compare(final PlanarGrid o1, final PlanarGrid o2) {
            int halfPeri1 = o1.getWidth() + o1.getHeight();
            int halfPeri2 = o2.getWidth() + o2.getHeight();

            if (halfPeri1 < halfPeri2) {
                return -1;
            }
            if (halfPeri1 == halfPeri2) {
                return 0;
            }
            return 1;
        }
    }

    /**
     * Comparator for sorting {@link Polyomino polyominoes} based on the length of the perimeter of their bounding
     * boxes, but preferring polyominoes with the shortest short sides.
     */
    private class MinPerimeterComparatorWithShape implements Comparator<PlanarGrid> {
        @Override
        public int compare(final PlanarGrid o1, final PlanarGrid o2) {
            int width = o1.getWidth();
            int height = o1.getHeight();

            if (width < height) {
                width *= width;
            } else {
                height *= height;
            }

            int val1 = width + height;

            width = o2.getWidth();
            height = o2.getHeight();

            if (width < height) {
                width *= width;
            } else {
                height *= height;
            }

            int val2 = width + height;

            if (val1 < val2) {
                return -1;
            }
            if (val1 == val2) {
                return 0;
            }
            return 1;
        }
    }

    /**
     * Comparator for sorting {@link Polyomino polyominoes} based on the minimum number of their extensions.
     */
    private class MinNumOfExtensionsComparator implements Comparator<Polyomino> {
        @Override
        public int compare(final Polyomino o1, final Polyomino o2) {
            int numExt1 = o1.getPolyominoExtensions().size();
            int numExt2 = o2.getPolyominoExtensions().size();

            if (numExt1 < numExt2) {
                return -1;
            }
            if (numExt1 == numExt2) {
                return 0;
            }
            return 1;
        }
    }

    /**
     * Comparator for sorting {@link Polyomino polyominoes} based on the minimum number of their extension directions.
     */
    private class MinNumOfExtensionDirectionsComparator implements Comparator<Polyomino> {
        @Override
        public int compare(final Polyomino o1, final Polyomino o2) {
            Function<UniqueTriple<Direction, Integer, Integer>, Direction> detectDirections =
                    (UniqueTriple<Direction, Integer, Integer> polyExt) -> polyExt.getFirst();

            Set<Direction> dirSet =
                    o1.getPolyominoExtensions().parallelStream().map(detectDirections).collect(Collectors.toSet());
            int numDir1 = dirSet.size();

            dirSet = o2.getPolyominoExtensions().parallelStream().map(detectDirections).collect(Collectors.toSet());
            int numDir2 = dirSet.size();

            if (numDir1 < numDir2) {
                return -1;
            }
            if (numDir1 == numDir2) {
                return 0;
            }
            return 1;
        }
    }

    /**
     * Comparator partitioning {@link Polyomino polyominoes} into polyominoes that only have extensions going one
     * direction (greater) and all other cases (less).
     */
    private class SingleExtensionSideGreaterThanRestComparator implements Comparator<Polyomino> {
        @Override
        public int compare(final Polyomino o1, final Polyomino o2) {
            Function<UniqueTriple<Direction, Integer, Integer>, Direction> detectDirections =
                    (UniqueTriple<Direction, Integer, Integer> polyExt) -> polyExt.getFirst();

            Set<Direction> dirSet =
                    o1.getPolyominoExtensions().parallelStream().map(detectDirections).collect(Collectors.toSet());
            int numDir1 = dirSet.size();

            dirSet = o2.getPolyominoExtensions().parallelStream().map(detectDirections).collect(Collectors.toSet());
            int numDir2 = dirSet.size();

            // more and less than 1 Extensions are considered more minimal
            numDir1 = numDir1 == 1 ? 1 : 0;
            numDir2 = numDir2 == 1 ? 1 : 0;

            if (numDir1 < numDir2) {
                return -1;
            }
            if (numDir1 == numDir2) {
                return 0;
            }
            return 1;
        }
    }

    /**
     * Comparator partitioning {@link Polyomino polyominoes} into polyominoes that only have extensions going two
     * orthogonal directions (one horizontally, one vertically -- greater) and all other cases (less).
     */
    private class CornerCasesGreaterThanRestComparator implements Comparator<Polyomino> {
        @Override
        public int compare(final Polyomino o1, final Polyomino o2) {
            Function<UniqueTriple<Direction, Integer, Integer>, Direction> detectDirections =
                    (UniqueTriple<Direction, Integer, Integer> polyExt) -> polyExt.getFirst();

            Set<Direction> dirSet =
                    o1.getPolyominoExtensions().parallelStream().map(detectDirections).collect(Collectors.toSet());
            int numDir1 = dirSet.size();
            // more and less than 2 Extensions are considered more minimal
            numDir1 = numDir1 == 2 ? 1 : 0;
            if (numDir1 == 1) {
                if (dirSet.parallelStream().filter(x -> x.isHorizontal()).collect(Collectors.counting()) % 2 == 0) {
                    // Both extensions are horizontal or vertical and thus not a "corner" case
                    numDir1 = 0;
                }
            }

            dirSet = o2.getPolyominoExtensions().parallelStream().map(detectDirections).collect(Collectors.toSet());
            int numDir2 = dirSet.size();
            // more and less than 2 Extensions are considered more minimal
            numDir2 = numDir2 == 2 ? 1 : 0;
            if (numDir2 == 1) {
                if (dirSet.parallelStream().filter(x -> x.isHorizontal()).collect(Collectors.counting()) % 2 == 0) {
                    // Both extensions are horizontal or vertical and thus not a "corner" case
                    numDir2 = 0;
                }
            }

            if (numDir1 < numDir2) {
                return -1;
            }
            if (numDir1 == numDir2) {
                return 0;
            }
            return 1;
        }
    }

}
