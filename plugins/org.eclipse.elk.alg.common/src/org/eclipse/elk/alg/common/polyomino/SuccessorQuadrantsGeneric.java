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

import java.util.Set;
import java.util.function.BiFunction;
import java.util.function.Function;
import java.util.stream.Collectors;

import org.eclipse.elk.alg.common.polyomino.structures.Direction;
import org.eclipse.elk.alg.common.polyomino.structures.Polyomino;
import org.eclipse.elk.alg.common.utils.UniqueTriple;
import org.eclipse.elk.core.util.Pair;

/**
 * <p>
 * Implements one way of computing the successor of a point (x', y') in ℕ x ℕ by modifying an existing way, given as an
 * instance of {@link BiFunction}. The modification is as follows.
 * </p>
 * 
 * <p>
 * case: polyomino has no extensions, only NORTH and SOUTH extensions, only EAST and WEST extensions, or extensions in
 * all four directions: simply use the old function
 * </p>
 * <p>
 * case: polyomino has only NORTH extensions, or only WEST, NORTH and EAST extensions: only return values of the old
 * function, if the y-coordinate y of the result satisfies the constraint y<0, otherwise skip to the next candidate
 * until the constraint holds.
 * </p>
 * <p>
 * case: polyomino has only EAST extensions, or only NORTH, EAST and SOUTH extensions: only return values of the old
 * function, if the x-coordinate x of the result satisfies the constraint x>=0.
 * </p>
 * <p>
 * case: polyomino has only SOUTH extensions, or only EAST, SOUTH and WEST extensions: only return values of the old
 * function, if y>=0.
 * </p>
 * <p>
 * case: polyomino has only WEST extensions, or only SOUTH, WEST AND NORTH extensions: only return values of the old
 * function, if x<0.
 * </p>
 * <p>
 * case: polyomino has only NORTH and EAST extensions: only return values of the old function, if y<0, x>=0.
 * </p>
 * 
 * <p>
 * case: polyomino has only SOUTH and EAST extensions: only return values of the old function, if y>=0, x>=0.
 * </p>
 * 
 * <p>
 * case: polyomino has only SOUTH and WEST extensions: only return values of the old function, if y>=0, x<0.
 * </p>
 * 
 * <p>
 * case: polyomino has only NORTH and WEST extensions: only return values of the old function, if y<0, x<0.
 * </p>
 */
public class SuccessorQuadrantsGeneric
        implements BiFunction<Pair<Integer, Integer>, Polyomino, Pair<Integer, Integer>> {

    private Polyomino lastPoly;
    private boolean posX, posY, negX, negY;
    private BiFunction<Pair<Integer, Integer>, Polyomino, Pair<Integer, Integer>> costFun;

    /**
     * Constructor for a function computing the successor of a point (x', y') in ℕ x ℕ by modifying an existing way,
     * given as a {@link BiFunction}. Based on the directions the extensions of a polyomino are facing the number of
     * valid candidate positions are further restricted by only allowing certain quadrants of the underlying coordinate
     * system.
     */
    public SuccessorQuadrantsGeneric(
            final BiFunction<Pair<Integer, Integer>, Polyomino, Pair<Integer, Integer>> costFun) {
        this.costFun = costFun;
    }

    @Override
    public Pair<Integer, Integer> apply(final Pair<Integer, Integer> coords, final Polyomino poly) {
        if (!poly.equals(lastPoly)) {
            lastPoly = poly;
            Function<UniqueTriple<Direction, Integer, Integer>, Direction> detectDirections =
                    polyExt -> polyExt.getFirst();
            Set<Direction> dirSet =
                    poly.getPolyominoExtensions().parallelStream().map(detectDirections).collect(Collectors.toSet());
            posX = true;
            posY = true;
            negX = true;
            negY = true;

            boolean containsPos = dirSet.contains(Direction.NORTH);
            boolean containsNeg = dirSet.contains(Direction.SOUTH);
            // Beware of the non mathematical coordinate system (inverted y-axis)!
            if (containsPos && !containsNeg) {
                posY = false;
            }
            if (!containsPos && containsNeg) {
                negY = false;
            }

            containsPos = dirSet.contains(Direction.EAST);
            containsNeg = dirSet.contains(Direction.WEST);
            if (containsPos && !containsNeg) {
                negX = false;
            }
            if (!containsPos && containsNeg) {
                posX = false;
            }

        }

        Pair<Integer, Integer> nextCoords = costFun.apply(coords, poly);

        int newX = nextCoords.getFirst();
        int newY = nextCoords.getSecond();

        // (0,0) is part of the positive quadrant
        boolean invalid = false;
        if (newX < 0) {
            if (!negX) {
                invalid = true;
            }
        } else {
            if (!posX) {
                invalid = true;
            }
        }
        if (newY < 0) {
            if (!negY) {
                invalid = true;
            }
        } else {
            if (!posY) {
                invalid = true;
            }
        }

        if (invalid) {
            return apply(nextCoords, poly);
        } else {
            return nextCoords;
        }

    }

}
