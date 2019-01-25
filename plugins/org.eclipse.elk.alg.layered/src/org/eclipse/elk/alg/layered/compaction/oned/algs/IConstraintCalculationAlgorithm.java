/*******************************************************************************
 * Copyright (c) 2016 Kiel University and others.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the Eclipse Public License v1.0
 * which accompanies this distribution, and is available at
 * http://www.eclipse.org/legal/epl-v10.html
 *
 * Contributors:
 *     Kiel University - initial API and implementation
 *******************************************************************************/
package org.eclipse.elk.alg.layered.compaction.oned.algs;

import org.eclipse.elk.alg.layered.compaction.oned.OneDimensionalCompactor;

/**
 * An algorithm that calculates separation constraints in one dimension that are induced by a
 * set of rectangles in the plane.
 */
@FunctionalInterface
public interface IConstraintCalculationAlgorithm {

    /**
     * @param compactor
     *            the instance of the surrounding {@link OneDimensionalCompactor}.
     */
    void calculateConstraints(OneDimensionalCompactor compactor);
    
}