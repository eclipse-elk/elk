/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.common;

/**
 * Interface for a cost function used in the processing order phase.
 */
@FunctionalInterface
public interface ICostFunction {
    /** Cost of an edge.
     *
     * @param edge
     * @return cost associated with an edge created in the structure phase
     */
    double cost(TEdge edge);
}
