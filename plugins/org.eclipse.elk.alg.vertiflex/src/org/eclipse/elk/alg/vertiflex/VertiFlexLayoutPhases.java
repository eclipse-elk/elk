/*******************************************************************************
 * Copyright (c) 2023 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.alg.vertiflex;

/**
 * This algorithm's layout phases. Each phase has its own strategy factory.
 */
public enum VertiFlexLayoutPhases {
    
    /** Phase 1. Vertical placement of nodes.*/
    P1_NODE_Y_PLACEMENT,
    
    /** Phase 2. Horizontal placement of nodes.*/
    P2_NODE_RELATIVE_PLACEMENT,
    
    /** Phase 3. Computation of absolute coordinates.*/
    P3_NODE_ABSOLUTE_PLACEMENT,
    
    /** Phase 4. Edge routing.*/
    P4_EDGE_ROUTING;
    
}