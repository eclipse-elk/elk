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
 * Describes what interactive layout phases should take as reference point for comparison
 * of node positions.
 * 
 * @author cds
 */
public enum InteractiveReferencePoint {
    
    /** the node's center point. */
    CENTER,
    /** the node's top left corner. */
    TOP_LEFT;
    
}
