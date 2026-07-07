/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.radial.options;

import org.eclipse.elk.alg.radial.p1position.wedge.AnnulusWedgeByLeafs;
import org.eclipse.elk.alg.radial.p1position.wedge.AnnulusWedgeByNodeSpace;
import org.eclipse.elk.alg.radial.p1position.wedge.IAnnulusWedgeCriteria;

/**
 * The list of possible wedge criteria.
 */
public enum AnnulusWedgeCriteria {

    /**
     * Create the wedge by the number of leaves a child has.
     */
    LEAF_NUMBER,

    /**
     * Create the wedge by the size a node requires.
     */
    NODE_SIZE;

    /**
     * Instantiate the chosen wedge strategy.
     * 
     * @return A wedge compactor.
     */
    public IAnnulusWedgeCriteria create() {
        switch (this) {
        case LEAF_NUMBER:
            return new AnnulusWedgeByLeafs();
        case NODE_SIZE:
            return new AnnulusWedgeByNodeSpace();
        default:
            throw new IllegalArgumentException(
                    "No implementation is available for the layout option " + this.toString());
        }
    }
}