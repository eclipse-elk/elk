/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.common.spore;

/**
 * Interface for an overlap handling function used in the {@link ScanlineOverlapCheck}.
 */
@FunctionalInterface
public interface IOverlapHandler {
    /**
     * This is what is done with overlapping nodes.
     * 
     * @param n1 one node
     * @param n2 the other node
     */
    void handle(Node n1, Node n2);
}
