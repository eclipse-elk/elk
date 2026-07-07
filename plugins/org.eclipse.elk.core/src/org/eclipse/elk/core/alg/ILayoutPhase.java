/*******************************************************************************
 * Copyright (c) 2010, 2017 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.core.alg;

/**
 * A layout phase processes a graph and may depend on {@link ILayoutProcessor layout processors} for further processing
 * of the graph.
 * 
 * <p>
 * Layout phases are the primary component of layout algorithms, providing their main structure. Each phase may have
 * different implementations, which are classes that implement this interface. Phases can require pre- or
 * post-processing by {@link ILayoutProcessor layout processors} executed in processing slots before, after, and between
 * the algorithm's phases.
 * </p>
 * 
 * <p>
 * The {@link AlgorithmAssembler} class can be used to build algorithms by specifying phases and letting the assembler
 * worry about instantiating all required processors.
 * </p>
 * 
 * @param <P>
 *            enumeration of all available phases. This is not an enumeration of all phase implementations.
 * @param <G>
 *            type of the graph the created phase will operate on.
 * @see AlgorithmAssembler
 */
public interface ILayoutPhase<P extends Enum<P>, G> extends ILayoutProcessor<G> {

    /**
     * Returns a layout processor configuration that specifies which {@link ILayoutProcessor layout processors} this
     * phase would require to be executed at which point in the algorithm to process the given graph.
     * 
     * @param graph
     *            the graph for which to return a configuration.
     * @return a configuration. Returning {@code null} is synonymous to returning an empty configuration (that is, no
     *         dependencies).
     */
    LayoutProcessorConfiguration<P, G> getLayoutProcessorConfiguration(G graph);

}
