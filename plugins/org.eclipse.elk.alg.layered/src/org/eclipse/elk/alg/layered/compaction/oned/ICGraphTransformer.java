/*******************************************************************************
 * Copyright (c) 2016 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.layered.compaction.oned;

/**
 * Interface for classes that import a {@link CGraph}.
 * 
 * @param <T> the type of data structure that is transformed into a {@link CGraph}
 */
public interface ICGraphTransformer<T> {
    /**
     * Transforms the input graph into a {@link CGraph} consisting of {@link CNode}s that may
     * be grouped in {@link CGroup}s.
     * 
     * @param inputGraph
     *          the graph to transform into a {@link CGraph}
     * @return a {@link CGraph}
     */
    CGraph transform(T inputGraph);
    
    /**
     * Updates the properties of the input graph and applies the compacted positions to the
     * elements of the input graph.
     */
    void applyLayout();
    
}
