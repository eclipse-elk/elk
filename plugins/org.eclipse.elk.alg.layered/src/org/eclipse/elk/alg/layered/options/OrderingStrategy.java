/*******************************************************************************
 * Copyright (c) 2020 Kiel University and others.
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
 * Strategy to order nodes and ports before crossing minimization.
 */
public enum OrderingStrategy {
    /**
     * Nothing is ordered.
     */
    NONE,
    /**
     * Nodes and edges are ordered.
     */
    NODES_AND_EDGES,
    /**
     * The node ordering is only used as a secondary criterion. Edge order shall be preserved.
     */
    PREFER_EDGES,
    /**
     * Prefers node order.
     */
    PREFER_NODES;

}
