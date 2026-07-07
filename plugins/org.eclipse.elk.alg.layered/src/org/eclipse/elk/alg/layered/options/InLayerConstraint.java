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
 * Enumeration of in-layer constraint types. In-layer constraints divide a layer into three
 * parts: the normal part, a top part and a bottom part. This constraint can be set on nodes
 * to define in which part they may appear.
 * 
 * @see org.eclipse.elk.alg.layered.intermediate.InLayerConstraintProcessor
 *        InLayerConstraintProcessor
 * @author cds
 */
public enum InLayerConstraint {
    
    /** no constraint on in-layer placement. */
    NONE,
    /** float node to the top of the layer, along with other nodes posessing this constraint. */
    TOP,
    /** float node to the bottom of the layer, along with other nodes posessing this constraint. */
    BOTTOM;
    
}
