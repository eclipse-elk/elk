/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.radial.options;

import org.eclipse.elk.alg.radial.sorting.IDSorter;
import org.eclipse.elk.alg.radial.sorting.IRadialSorter;
import org.eclipse.elk.alg.radial.sorting.PolarCoordinateSorter;
import org.eclipse.elk.graph.properties.AdvancedPropertyValue;

/**
 *
 */
public enum SortingStrategy {
    /** Do no sorting. */
    NONE,
    /** Sort by polar coordinates of parent. */
    @AdvancedPropertyValue
    POLAR_COORDINATE,
    /** Sort by given order id. */
    ID;

    /**
     * Instantiate the chosen wedge strategy.
     * 
     * @return A wedge compactor.
     */
    public IRadialSorter create() {
        switch (this) {
        case NONE:
            return null;
        case POLAR_COORDINATE:
            return new PolarCoordinateSorter();
        case ID:
            return new IDSorter();
        default:
            throw new IllegalArgumentException(
                    "No implementation is available for the layout option " + this.toString());
        }
    }
}