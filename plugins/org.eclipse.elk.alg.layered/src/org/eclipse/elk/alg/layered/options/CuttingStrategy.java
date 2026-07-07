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

/**
 * Specifies the strategy employed to calculate cut indexes during graph wrapping.
 * 
 * @see WrappingStrategy
 */
public enum CuttingStrategy {

    /** Aspect ratio-driven cut calculation heuristic. */
    ARD,
    /** Max scale-driven cut calculation heuristic. */
    MSD,
    /** Cuts are manually specified by a user via {@link LayeredOptions#WRAPPING_CUTTING_CUTS}. */
    MANUAL,
    
}
