/*******************************************************************************
 * Copyright (c) 2014, 2020 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.layered.options;

/**
 * Sets the variant of the greedy switch heuristic.
 * 
 * Note that there are two layout options using this enumeration:
 * {@link LayeredOptions#CROSSING_MINIMIZATION_GREEDY_SWITCH_TYPE} and
 * {@link LayeredOptions#CROSSING_MINIMIZATION_GREEDY_SWITCH_HIERARCHICAL_TYPE}.
 */
public enum GreedySwitchType {
    /** Only consider crossings to one side of the free layer. Calculate crossing matrix on demand. */
    ONE_SIDED,
    /** Consider crossings to both sides of the free layer. Calculate crossing matrix on demand. */
    TWO_SIDED,
    /** Do not use greedy switch. */
    OFF;
}
