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
 * Enumeration of edge constraints. Edge constraints can be set on ports to constrain the
 * type of edges incident to that port.
 * 
 * @see org.eclipse.elk.alg.layered.intermediate.EdgeAndLayerConstraintEdgeReverser
 *        EdgeAndLayerConstraintEdgeReverser
 * @author cds
 */
public enum EdgeConstraint {
    
    /** no constraint on incident edges. */
    NONE,
    /** node may have only incoming edges. */
    INCOMING_ONLY,
    /** node may have only outgoing edges. */
    OUTGOING_ONLY;
    
}
