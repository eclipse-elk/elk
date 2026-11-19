/*******************************************************************************
 * Copyright (c) 2015, 2017 Kiel University and others.
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

/**
 * Classes that implement this interface can create layout phases.
 * 
 * <p>
 * Algorithms that are structured into phases often have different implementations available for each phase. Thus, each
 * phase will provide a custom enumeration of all available implementations. That implementation is the prime candidate
 * for implementing this interface and instantiating the correct phase implementation subject to the enumeration
 * constant that was called upon to do so.
 * </p>
 * 
 * <p>
 * If there is only a single implementation of each phase, the enumeration that lists all phases can serve as the phase
 * factory.
 * </p>
 * 
 * @param <P>
 *            enumeration of all available phases. This is not an enumeration of all phase implementations.
 * @param <G>
 *            type of the graph the created phase will operate on.
 * @see ILayoutPhase
 */
public interface ILayoutPhaseFactory<P extends Enum<P>, G> extends ILayoutProcessorFactory<G> {

    /**
     * Returns an implementation of {@link ILayoutPhase}. The actual implementation returned depends on the actual type
     * that implements this method.
     * 
     * @return new layout processor.
     */
    ILayoutPhase<P, G> create();

}
