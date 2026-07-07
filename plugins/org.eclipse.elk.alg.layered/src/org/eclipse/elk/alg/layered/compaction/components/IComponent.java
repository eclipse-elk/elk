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

import java.util.List;
import java.util.Set;

import org.eclipse.elk.core.math.ElkRectangle;
import org.eclipse.elk.core.options.PortSide;

/**
 * A component represents, for instance, a set of nodes that are fully connected. In other words,
 * between any pair of nodes of the component exists an undirected path.
 * 
 * @param <N>
 *            the type of the contained nodes
 * @param <E>
 *            the type of the contained edges
 */
public interface IComponent<N, E> {

    /**
     * @return a hull the nodes of this connected component. This, for instance, can be a
     *         rectilinear convex hull or a polyomino.
     */
    List<ElkRectangle> getHull();

    /**
     * @return a list with external extensions of this component.
     */
    List<IExternalExtension<E>> getExternalExtensions();
    
    /**
     * @return a set of the port sides of this component's external ports (if there are any,
     *         otherwise an empty set).
     */
    Set<PortSide> getExternalExtensionSides();

}
