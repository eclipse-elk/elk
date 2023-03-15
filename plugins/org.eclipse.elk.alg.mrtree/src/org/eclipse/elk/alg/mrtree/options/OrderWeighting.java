/*******************************************************************************
 * Copyright (c) 2013 - 2022 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.alg.mrtree.options;

/**
 * Layout option for the choice of the weighting for the node order.
 *
 * @author sor
 * @author sgu
 * @author sdo
 */
public enum OrderWeighting {
    
    /** Doesn't change node ordering. 
     * When using this nodes will be in the order of the definition of their incoming edges */
    MODEL_ORDER,
    /** Chooses the number of descendants for the weighting. */
    DESCENDANTS,
    /** Chooses maximal number of descendants at the same level. */
    FAN,
    /** Chooses the nodes with the highest position constraint. 
     * Without set constraints this will behave just like NONE */
    CONSTRAINT;
    
}
