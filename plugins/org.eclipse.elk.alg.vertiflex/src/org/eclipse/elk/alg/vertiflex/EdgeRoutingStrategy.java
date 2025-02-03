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
 * Strategies for routing edges in the tree layout.
 *
 */
public enum EdgeRoutingStrategy {
    
    /**
     * Straight lines between nodes.
     */
    STRAIGHT,
    
    /**
     * Allow one bend point in edges to enable more compact layouts while also maintaining the node model order.
     */
    BEND

}
