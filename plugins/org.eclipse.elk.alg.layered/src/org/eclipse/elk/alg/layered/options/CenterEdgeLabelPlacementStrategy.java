/*******************************************************************************
 * Copyright (c) 2016, 2017 Kiel University and others.
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
package org.eclipse.elk.alg.layered.options;

import org.eclipse.elk.graph.properties.ExperimentalPropertyValue;

/**
 * Strategies regarding how to place center edge labels. Implementation-wise, the placement strategy influences which
 * layer a center edge label dummy is moved into.
 */
public enum CenterEdgeLabelPlacementStrategy {
    
    /** Consider all layers a labeled edge crosses and place its labels in the one at the median spot. */
    MEDIAN_LAYER,
    
    /**
     * Places labels in the layer closest to the edge's tail. This is similar to simply defining tail labels, but is a
     * bit more space-efficient.
     */
    TAIL_LAYER,
    
    /**
     * Places labels in the layer closest to the edge's head. This is similar to simply defining tail labels, but is a
     * bit more space-efficient.
     */
    HEAD_LAYER,
    
    /** Places labels in layers that will likely result in the most compact diagram. */
    SPACE_EFFICIENT_LAYER,
    
    /** Place labels in the widest of all layers its edge crosses. */
    @ExperimentalPropertyValue
    WIDEST_LAYER,
    
    /** Places labels in a layer that we think will be closest to the edge's physical center later. */
    @ExperimentalPropertyValue
    CENTER_LAYER;
    
    
    /**
     * Checks whether this strategy is based on information about the width of layers.
     * 
     * @return {@code true} iff that is the case.
     */
    public boolean usesLabelSizeInformation() {
        return this == WIDEST_LAYER
                || this == CENTER_LAYER
                || this == SPACE_EFFICIENT_LAYER;
    }
    
}