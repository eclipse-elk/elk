/*******************************************************************************
 * Copyright (c) 2024 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 * 
 * SPDX-License-Identifier: EPL-2.0 
 *******************************************************************************/
package org.eclipse.elk.core.options;

import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.graph.ElkNode;

/**
 * A topdown size approximator returns an estimated size of the graph drawing after performing layout using some 
 * heuristic.
 *
 */
public interface ITopdownSizeApproximator {

    /**
     * Returns an approximated required size for a given node.
     * @param node the node
     * @return the size as a vector
     */
    public KVector getSize(ElkNode node);
}
