/*******************************************************************************
 * Copyright (c) 2022 - 2023 Kiel University and others.
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
package org.eclipse.elk.alg.topdownpacking;

import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.graph.ElkNode;

/**
 * A node arranger has to provide a size prediction without computing a layout.
 */
public interface INodeArranger {
    
    /**
     * Computes the predicted required size of the graph and returns it without computing a full layout.
     * @param graph the input graph
     * @return the predicted size
     */
    KVector getPredictedSize(ElkNode graph);

}
