/*******************************************************************************
 * Copyright (c) 2011, 2015 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.layered.options;

/**
 * Enumeration of layer constraint types. May be set on nodes to constrain in which layer
 * they may appear.
 *
 * @see org.eclipse.elk.alg.layered.intermediate.EdgeAndLayerConstraintEdgeReverser
 *        EdgeAndLayerConstraintEdgeReverser
 * @see org.eclipse.elk.alg.layered.intermediate.LayerConstraintProcessor
 *        LayerConstraintProcessor
 * @author msp
 */
public enum LayerConstraint {
    
    /** no constraint on the layering. */
    NONE,
    /** put into the first layer. */
    FIRST,
    /** put into a separate first layer; used internally. */
    FIRST_SEPARATE,
    /** put into the last layer. */
    LAST,
    /** put into a separate last layer; used internally. */
    LAST_SEPARATE;
    
}
