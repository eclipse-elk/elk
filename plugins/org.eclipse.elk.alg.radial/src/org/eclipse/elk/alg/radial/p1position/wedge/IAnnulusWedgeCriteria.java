/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.radial.p1position.wedge;

import org.eclipse.elk.graph.ElkNode;

/**
 * An interface for defining how the annulus wedge is determined.
 */
public interface IAnnulusWedgeCriteria {

    /**
     * Calculate the space for a wedge.
     * 
     * @param node
     *            A node which shall be placed in a wedge.
     * @return The space a node requires for its wedge.
     */
    double calculateWedgeSpace(ElkNode node);

}