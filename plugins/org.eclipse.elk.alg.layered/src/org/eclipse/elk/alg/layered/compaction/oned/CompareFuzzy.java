/*******************************************************************************
 * Copyright (c) 2016 Kiel University and others.
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
package org.eclipse.elk.alg.layered.compaction.oned;


/**
 * Internal Class for tolerance affected double comparisons.
 */
public final class CompareFuzzy {
    /** Default fuzzy comparison tolerance. */
    public static final double TOLERANCE = 0.0001;

    private CompareFuzzy() {
    }

    /** Tolerance-aware equality: {@code true} iff |a-b| &le; tolerance, with NaN handling. */
    public static boolean fuzzyEquals(final double a, final double b, final double tolerance) {
        return Math.abs(a - b) <= tolerance
                || Double.valueOf(a).equals(Double.valueOf(b));
    }

    /** Tolerance-aware comparison consistent with {@link #fuzzyEquals}. */
    public static int fuzzyCompare(final double a, final double b, final double tolerance) {
        if (fuzzyEquals(a, b, tolerance)) {
            return 0;
        }
        if (a < b) {
            return -1;
        }
        if (a > b) {
            return 1;
        }
        return Boolean.compare(Double.isNaN(a), Double.isNaN(b));
    }

    // SUPPRESS CHECKSTYLE NEXT 20 Javadoc
    public static boolean eq(final double d1, final double d2) {
        return fuzzyEquals(d1, d2, TOLERANCE);
    }

    public static boolean gt(final double d1, final double d2) {
        return fuzzyCompare(d1, d2, TOLERANCE) > 0;
    }

    public static boolean lt(final double d1, final double d2) {
        return fuzzyCompare(d1, d2, TOLERANCE) < 0;
    }

    public static boolean ge(final double d1, final double d2) {
        return fuzzyCompare(d1, d2, TOLERANCE) >= 0;
    }

    public static boolean le(final double d1, final double d2) {
        return fuzzyCompare(d1, d2, TOLERANCE) <= 0;
    }
}
