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
package org.eclipse.elk.core.util;

import java.util.Objects;

/**
 * Data structure for storing three (different) values.
 *
 * @param <F>
 *            type of first value
 * @param <S>
 *            type of second value
 * @param <T>
 *            type of third value
 */
public class Triple<F, S, T> {
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
    public Triple(final F f, final S s, final T t) {
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

    /*
     * (non-Javadoc)
     * 
     * @see java.lang.Object#hashCode()
     */
    @Override
    public int hashCode() {
        return Objects.hash(first, second, third);
    }

    /*
     * (non-Javadoc)
     * 
     * @see java.lang.Object#equals(java.lang.Object)
     */
    @Override
    public boolean equals(final Object obj) {
        return Objects.equals(first, ((Triple<?, ?, ?>) obj).getFirst())
                && Objects.equals(second, ((Triple<?, ?, ?>) obj).getSecond())
                && Objects.equals(third, ((Triple<?, ?, ?>) obj).getThird());
    }
    
    @Override
    public String toString() {
        return "(" + first + ", " + second + ", " + third + ")";
    }

}
