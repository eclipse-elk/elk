/*******************************************************************************
 * Copyright (c) 2016 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.layered.compaction.components;

import org.eclipse.elk.core.math.ElkRectangle;
import org.eclipse.elk.core.options.Direction;

/**
 * Container class for edges that represent external edges, i.e. edges that connect to ports that
 * span hierarchy boundaries.
 * 
 * @param <E>
 *            the type of the element that is being represented by this extension.
 */
public interface IExternalExtension<E> {

    /**
     * @return for instance the underlying edge this extension represents.
     */
    E getRepresentative();

    /**
     * @return the rectangle that represents the external extension.
     */
    ElkRectangle getRepresentor();
    
    /**
     * @return an optional placeholder along the original diagram's boundary.
     */
    default ElkRectangle getPlaceholder() {
        return null;
    }
    
    /**
     * @return the rectangle to which this extension (conceptually) connects.
     */
    ElkRectangle getParent();

    /**
     * @return the direction into which this extension points.
     */
    Direction getDirection();

}
