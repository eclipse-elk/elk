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

/**
 * A set of {@link IComponent}s. 
 * 
 * @param <N>
 *            the type of the contained nodes
 * @param <E>
 *            the type of the contained edges
 */
public interface IConnectedComponents<N, E> extends Iterable<IComponent<N, E>> {

    /**
     * @return the components
     */
    List<IComponent<N, E>> getComponents();

    /**
     * @return whether any of the components contains external extensions
     */
    boolean isContainsExternalExtensions();

}
