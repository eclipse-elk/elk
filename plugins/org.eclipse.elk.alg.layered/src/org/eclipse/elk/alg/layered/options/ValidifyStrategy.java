/*******************************************************************************
 * Copyright (c) 2016 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.layered.options;

import org.eclipse.elk.alg.layered.intermediate.wrapping.ICutIndexCalculator;

/**
 * The cuts calculated by a {@link ICutIndexCalculator} are not necessarily valid (when using
 * {@link WrappingStrategy#SINGLE_EDGE}). The specified 'validify' takes care of transforming forbidden cuts into valid
 * ones.
 */
public enum ValidifyStrategy {
    
    /** Do not touch my cuts!. */
    NO,
    /** Just increase forbidden cuts until they are valid. */
    GREEDY,
    /** Be a bit smarter and check if the lastly valid cut is closer than the next valid cut. */
    LOOK_BACK,
    
}
