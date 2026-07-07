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
 * Implements one way of computing the successor of a point (x, y) in ℕ x ℕ by combining two other traversal functions,
 * given as instances of {@link BiFunction}. {@code normalFun} is used on polyominoes without external extensions,
 * {@code externalFun} is used on polyominoes with external extensions.
 * </p>
 */
public class SuccessorCombination implements BiFunction<Pair<Integer, Integer>, Polyomino, Pair<Integer, Integer>> {

    private BiFunction<Pair<Integer, Integer>, Polyomino, Pair<Integer, Integer>> normalFun;
    private BiFunction<Pair<Integer, Integer>, Polyomino, Pair<Integer, Integer>> externalFun;

    /**
     * Constructor for a function computing the successor of a point (x, y) in ℕ x ℕ by combining two other traversal
     * functions.
     * 
     * @param normalFun
     *            traversal function used on polyominoes without external extensions
     * @param externalFun
     *            traversal function used on polyominoes with external extensions
     */
    public SuccessorCombination(final BiFunction<Pair<Integer, Integer>, Polyomino, Pair<Integer, Integer>> normalFun,
            final BiFunction<Pair<Integer, Integer>, Polyomino, Pair<Integer, Integer>> externalFun) {
        this.normalFun = normalFun;
        this.externalFun = externalFun;
    }

    @Override
    public Pair<Integer, Integer> apply(final Pair<Integer, Integer> coords, final Polyomino poly) {
        if (poly.getPolyominoExtensions().size() > 0) {
            return externalFun.apply(coords, poly);
        } else {
            return normalFun.apply(coords, poly);
        }
    }

}
