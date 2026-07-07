/*******************************************************************************
 * Copyright (c) 2022 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.mrtree.options;

/**
 * 
 * @author jnc, sdo
 */
public enum EdgeRoutingMode {
    /** Does nothing. */
    NONE,
    /** Routes edges from node middle to node middle. */
    MIDDLE_TO_MIDDLE,
    /** Similar to MiddleToMiddle but avoids edge node overlaps. */
    AVOID_OVERLAP
}
