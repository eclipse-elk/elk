/*******************************************************************************
 * Copyright (c) 2020 sdo and others.
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



/**
 * Strategy to sort dummy nodes compared to nodes with no connection to the previous layer.
 * Dummy nodes nodes are not part of the input model and 
 */
public enum LongEdgeOrderingStrategy {
    /**
     * Dummy nodes are sorted over normal nodes.
     */
    DUMMY_NODE_OVER,
    /**
     * Dummy nodes are sorted under normal nodes.
     */
    DUMMY_NODE_UNDER,
    /**
     * Dummy nodes are equal to normal nodes.
     */
    EQUAL;
    
    /**
     * Returns a value for a comparator to implement the desired sorting strategy.
     * @return a model order node comparator value
     */
    public int returnValue() {
        switch (this) {
        case DUMMY_NODE_OVER:
            return Integer.MAX_VALUE;
        case DUMMY_NODE_UNDER:
            return Integer.MIN_VALUE;
        default:
            return 0;
        }
        
    }
}
