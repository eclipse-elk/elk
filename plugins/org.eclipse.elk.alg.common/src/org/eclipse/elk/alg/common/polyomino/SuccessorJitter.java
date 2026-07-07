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
 * to enumerate these points. In this implementation the order of the positions is supposed to look a little bit random
 * starting in the middle of each side clockwise and after each rotation stepping a step further to the corners where
 * the sides meet (thus the name "jittery").
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
 * {@code _ X _}</br>
 * {@code _ O _}</br>
 * {@code _ _ _}</br>
 * </br>
 * {@code _ O _}</br>
 * {@code _ O X}</br>
 * {@code _ _ _}</br>
 * </br>
 * {@code _ O _}</br>
 * {@code _ O O}</br>
 * {@code _ X _}</br>
 * </br>
 * {@code _ O _}</br>
 * {@code X O O}</br>
 * {@code _ O _}</br>
 * </br>
 * {@code _ O X}</br>
 * {@code O O O}</br>
 * {@code _ O _}</br>
 * </br>
 * {@code _ O O}</br>
 * {@code O O O}</br>
 * {@code _ O X}</br>
 * </br>
 * {@code _ O O}</br>
 * {@code O O O}</br>
 * {@code X O O}</br>
 * </br>
 * {@code X O O}</br>
 * {@code O O O}</br>
 * {@code O O O}</br>
 * </br>
 */
public class SuccessorJitter implements BiFunction<Pair<Integer, Integer>, Polyomino, Pair<Integer, Integer>> {

    @Override
    public Pair<Integer, Integer> apply(final Pair<Integer, Integer> coords, final Polyomino poly) {
        int x = coords.getFirst();
        int y = coords.getSecond();

        int newX = x;
        int newY = y;

        int cost = Math.max(Math.abs(x), Math.abs(y));

        if (x <= 0 && x == y) {
            newX = 0;
            newY = y - 1;
        } else {
            if (x == -cost && y != cost) {
                newX = y;
                newY = x;
                if (y >= 0) {
                    newX++;
                }
            } else {
                newX = -y;
                newY = x;
            }
        }
        return new Pair<Integer, Integer>(newX, newY);
    }

}
