/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
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
package org.eclipse.elk.alg.disco.options;

import org.eclipse.elk.alg.common.polyomino.PolyominoCompactor;

/**
 * <p>
 * Possible strategies for packing different connected components on a canvas in order to save space and enhance
 * readability of a graph.
 * </p>
 * 
 * <p>
 * POLYOMINO - First naive implementation of Freivalds et al. (2002) realized in {@link PolyominoCompactor}
 * </p>
 */
public enum CompactionStrategy {
    /**
     * There is only one strategy implemented, yet.
     */
    POLYOMINO
}
