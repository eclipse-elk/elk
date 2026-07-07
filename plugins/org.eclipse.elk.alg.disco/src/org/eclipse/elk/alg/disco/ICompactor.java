/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.disco;

import org.eclipse.elk.alg.disco.graph.DCComponent;
import org.eclipse.elk.alg.disco.graph.DCGraph;
import org.eclipse.elk.alg.disco.transform.IGraphTransformer;

/**
 * Interface for algorithms geared towards compacting the connected components of {@link DCGraph DCGraphs}.
 */
public interface ICompactor {

    /**
     * <p>
     * Compacts the connected components of a {@link DCGraph}. Sets an offset field of each {@link DCComponent}, i.e. a
     * vector relative to the absolute position of each component before compaction. Moreover sets a field of the
     * {@link DCGraph} containing the dimensions of a minimum bounding box able to depict the whole resulting graph.
     * </p>
     * 
     * <p>
     * The values of these fields are later going to be used in a class implementing {@link IGraphTransformer} to apply
     * the resulting layout to the original graph structure the {@link DCGraph} was created from.
     * </p>
     * 
     * @param graph
     *            DCGraph to compact
     */
    void compact(DCGraph graph);
}