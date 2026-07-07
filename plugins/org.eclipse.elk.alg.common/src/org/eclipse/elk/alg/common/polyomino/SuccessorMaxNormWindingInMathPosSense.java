/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.common.polyomino;

import java.util.function.BiFunction;

import org.eclipse.elk.alg.common.polyomino.structures.Polyomino;
import org.eclipse.elk.core.util.Pair;

/**
 * <p>
 * Implements one way of computing the successor of a point (x, y) in ℕ x ℕ based on enumerating points in increasing
 * order of the cost function max(|x|, |y|) given in Freivalds et al. (2002) (see {@link PolyominoCompactor} for full
 * bibliographic entry). As the cost function returns the same value for some distinct points, there are different ways
 * to enumerate these points. In this implementation the order of the positions forms a spiral winding in the
 * mathematical positive sense.
 * </p>
 * Example: a 3x3 grid with coordinates </br>
 * {@code (-1,-1) ( 0,-1) ( 1,-1)}</br>
 * {@code (-1, 0) ( 0, 0) ( 1, 0)}</br>
 * {@code (-1, 1) ( 0, 1) ( 1, 1)}</br>
 * will be enumerated in this order (marked by an X) </br>
 * {@code _ _ _}</br>
 * {@code _ X _}</br>
 * {@code _ _ _}</br>
 * </br>
 * {@code _ _ _}</br>
 * {@code X O _}</br>
 * {@code _ _ _}</br>
 * </br>
 * {@code X _ _}</br>
 * {@code O O _}</br>
 * {@code _ _ _}</br>
 * </br>
 * {@code O X _}</br>
 * {@code O O _}</br>
 * {@code _ _ _}</br>
 * </br>
 * {@code O O X}</br>
 * {@code O O _}</br>
 * {@code _ _ _}</br>
 * </br>
 * {@code O O O}</br>
 * {@code O O X}</br>
 * {@code _ _ _}</br>
 * </br>
 * {@code O O O}</br>
 * {@code O O O}</br>
 * {@code _ _ X}</br>
 * </br>
 * {@code O O O}</br>
 * {@code O O O}</br>
 * {@code _ X O}</br>
 * </br>
 * {@code O O O}</br>
 * {@code O O O}</br>
 * {@code X O O}</br>
 * </br>
 */
public class SuccessorMaxNormWindingInMathPosSense
        implements BiFunction<Pair<Integer, Integer>, Polyomino, Pair<Integer, Integer>> {

    @Override
    public Pair<Integer, Integer> apply(final Pair<Integer, Integer> coords, final Polyomino poly) {
        int x = coords.getFirst();
        int y = coords.getSecond();
        int cost = Math.max(Math.abs(x), Math.abs(y));
        if (x < cost && y == -cost) {
            return new Pair<Integer, Integer>(x + 1, y);
        }
        if (x == cost && y < cost) {
            return new Pair<Integer, Integer>(x, y + 1);
        }
        if (x >= -cost && y == cost) {
            return new Pair<Integer, Integer>(x - 1, y);
        }
        // if (x == -cost && y > -cost)
        return new Pair<Integer, Integer>(x, y - 1);
    }

}
