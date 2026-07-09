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
 * Data structure for storing four (different) values.
 * 
 * @param <A>
 *            type of first value
 * @param <B>
 *            type of second value
 * @param <C>
 *            type of third value
 * @param <D>
 *            type of fourth value
 */
public class Quadruple<A, B, C, D> {
    private A first;
    private B second;
    private C third;
    private D fourth;

    /**
     * Creates a Quadruple with initial values. These are read-only afterwards.
     * 
     * @param a
     *            First value
     * @param b
     *            Second value
     * @param c
     *            Third value
     * @param d
     *            Fourth value
     */
    public Quadruple(final A a, final B b, final C c, final D d) {
        first = a;
        second = b;
        third = c;
        fourth = d;
    }

    /**
     * Gets the first value.
     * 
     * @return The first value
     */
    public A getFirst() {
        return first;
    }

    /**
     * Gets the second value.
     * 
     * @return The second value
     */
    public B getSecond() {
        return second;
    }

    /**
     * Gets the third value.
     * 
     * @return The third value
     */
    public C getThird() {
        return third;
    }

    /**
     * Gets the fourth value.
     * 
     * @return The fourth value
     */
    public D getFourth() {
        return fourth;
    }

    /*
     * (non-Javadoc)
     * 
     * @see java.lang.Object#hashCode()
     */
    @Override
    public int hashCode() {
        return Objects.hash(first, second, third, fourth);
    }

    /*
     * (non-Javadoc)
     * 
     * @see java.lang.Object#equals(java.lang.Object)
     */
    @Override
    public boolean equals(final Object obj) {
        return Objects.equals(first, ((Quadruple<?, ?, ?, ?>) obj).getFirst())
                && Objects.equals(second, ((Quadruple<?, ?, ?, ?>) obj).getSecond())
                && Objects.equals(third, ((Quadruple<?, ?, ?, ?>) obj).getThird())
                && Objects.equals(fourth, ((Quadruple<?, ?, ?, ?>) obj).getFourth());
    }

    @Override
    public String toString() {
        return "(" + first + ", " + second + ", " + third + ", " + fourth + ")";
    }
}
