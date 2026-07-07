/*******************************************************************************
 * Copyright (c) 2023 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 * 
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later 
 *******************************************************************************/
package org.eclipse.elk.alg.radial.intermediate.rotation;

import org.eclipse.elk.graph.ElkNode;

/**
 * An interface for rotating the radial layout.
 *
 */
public interface IRadialRotator {
    
    /**
     * Rotate the graph.
     * 
     * @param graph
     *            The graph which is already radial.
     */
    void rotate(ElkNode graph);

}