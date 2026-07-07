/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.common.compaction.oned;

import com.google.common.math.DoubleMath;

/**
 * Internal Class for tolerance affected double comparisons.
 */
public final class CompareFuzzy {
    /** Epsilon for double compares. */
    public static final double TOLERANCE = 0.0001;
    
    private CompareFuzzy() {
    }

    // SUPPRESS CHECKSTYLE NEXT 20 Javadoc
    public static boolean eq(final double d1, final double d2) {
        return DoubleMath.fuzzyEquals(d1, d2, TOLERANCE);
    }

    public static boolean gt(final double d1, final double d2) {
        return DoubleMath.fuzzyCompare(d1, d2, TOLERANCE) > 0;
    }

    public static boolean lt(final double d1, final double d2) {
        return DoubleMath.fuzzyCompare(d1, d2, TOLERANCE) < 0;
    }

    public static boolean ge(final double d1, final double d2) {
        return DoubleMath.fuzzyCompare(d1, d2, TOLERANCE) >= 0;
    }

    public static boolean le(final double d1, final double d2) {
        return DoubleMath.fuzzyCompare(d1, d2, TOLERANCE) <= 0;
    }
}
