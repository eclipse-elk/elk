/*******************************************************************************
 * Copyright (c) 2011, 2015 Kiel University and others.
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
package org.eclipse.elk.alg.force.options;

import java.util.Random;

import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.graph.properties.IProperty;
import org.eclipse.elk.graph.properties.Property;

/**
 * Internal property definitions for the force layout algorithm.
 * 
 * @author msp
 */
public final class InternalProperties {
    
    /** the original object from which a graph element was created. */
    public static final IProperty<Object> ORIGIN = new Property<Object>("origin");
    
    /** random number generator for the algorithm. */
    public static final IProperty<Random> RANDOM = new Property<Random>("random");
    
    /** upper left corner of the graph's bounding box. */
    public static final IProperty<KVector> BB_UPLEFT = new Property<KVector>("boundingBox.upLeft");
    /** lower right corner of the graph's bounding box. */
    public static final IProperty<KVector> BB_LOWRIGHT = new Property<KVector>("boundingBox.lowRight");
    
    /**
     * Hidden default constructor.
     */
    private InternalProperties() {
    }

}
