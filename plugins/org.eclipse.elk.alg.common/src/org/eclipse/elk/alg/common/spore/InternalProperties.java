/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.common.spore;

import org.eclipse.elk.graph.properties.IProperty;
import org.eclipse.elk.graph.properties.Property;

/**
 * Container for property definitions for internal use.
 */
public final class InternalProperties {
    /** Fuzziness of fuzzy comparisons used to avoid non-deterministic behavior caused by double imprecision. */
    public static final double FUZZINESS = 0.0001;
    
    // CHECKSTYLEOFF VisibilityModifier
    /** Enables SVG debug output. */
    public static final IProperty<Boolean> DEBUG_SVG = new Property<Boolean>("debugSVG", false);
    
    /** This flag is set in {@link GrowTreePhase} if two nodes of an edge in the minimum spanning tree 
     *  were overlapping. 
     */
    public static final IProperty<Boolean> OVERLAPS_EXISTED = new Property<Boolean>("overlapsExisted", true);
    
    /**
     * hidden default constructor.
     */
    private InternalProperties() {
        
    }
}
