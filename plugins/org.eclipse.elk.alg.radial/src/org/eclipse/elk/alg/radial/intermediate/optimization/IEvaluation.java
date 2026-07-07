/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.radial.intermediate.optimization;

import org.eclipse.elk.graph.ElkNode;

/**
 * The interface for translation optimizations.
 */
public interface IEvaluation {

    /**
     * Evaluate the current graph status.
     * 
     * @param graph The graph.
     * @return a value representing the score of the evaluation
     */
    double evaluate(ElkNode graph);
}