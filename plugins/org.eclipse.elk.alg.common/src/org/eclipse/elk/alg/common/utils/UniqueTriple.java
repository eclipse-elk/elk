/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.common.utils;

/**
 * Data structure for storing three (different) values. Each triple is unique, two different instantiations will be
 * unequal, even if filled with the same three objects.
 *
 * @param <F>
 *            type of first value
 * @param <S>
 *            type of second value
 * @param <T>
 *            type of third value
 */
public class UniqueTriple<F, S, T> {
    private F first;
    private S second;
    private T third;

    /**
     * Creates a Triple with initial values. These are read-only afterwards.
     * 
     * @param f
     *            First value
     * @param s
     *            Second value
     * @param t
     *            Third value
     */
    public UniqueTriple(final F f, final S s, final T t) {
        first = f;
        second = s;
        third = t;
    }

    /**
     * Gets the first value.
     * 
     * @return The first value
     */
    public F getFirst() {
        return first;
    }

    /**
     * Gets the second value.
     * 
     * @return The second value
     */
    public S getSecond() {
        return second;
    }

    /**
     * Gets the third value.
     * 
     * @return The third value
     */
    public T getThird() {
        return third;
    }
    
    @Override
    public String toString() {
        return "(" + first + ", " + second + ", " + third + ")";
    }

}
