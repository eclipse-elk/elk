/*******************************************************************************
 * Copyright (c) 2010, 2017 Kiel University and others.
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
package org.eclipse.elk.core.alg;

import org.eclipse.elk.core.util.IElkProgressMonitor;

/**
 * A layout processor processes a graph. Layout processors are the secondary components of layout algorithms, the
 * primary being {@link ILayoutPhase layout phases}. Layout processors are inserted before or after phases to do further
 * processing on a graph.
 * 
 * <p>
 * The {@link AlgorithmAssembler} class can be used to build algorithms by specifying phases and letting the assembler
 * worry about instantiating all required processors.
 * </p>
 *
 * @param <G>
 *            the type of graph the processor will process.
 * @see ILayoutPhase
 * @see AlgorithmAssembler
 */
public interface ILayoutProcessor<G> {

    /**
     * Performs the processor's work on the given graph.
     * 
     * @param graph
     *            the graph to process.
     * @param progressMonitor
     *            a progress monitor to indicate progress with.
     */
    void process(G graph, IElkProgressMonitor progressMonitor);

}
