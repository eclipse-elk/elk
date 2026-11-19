/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
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
package org.eclipse.elk.alg.disco.transform;

import org.eclipse.elk.alg.disco.graph.DCGraph;

/**
 * Interface for transforming any graph structure into a {@link DCGraph}, in order to use layout algorithms provided by
 * DisCo for compacting disconnected components of the given graph.
 *
 * @param <G>
 */
public interface IGraphTransformer<G> {

    /**
     * Transforms the given graph into a {@link DCGraph}.
     * 
     * @param graph
     * @return
     */
    DCGraph importGraph(G graph);

    /**
     * Applies the changes made to its associated {@link DCGraph} back to the original graph (presumably after a layout
     * algorithm has been used).
     */
    void applyLayout();
}
