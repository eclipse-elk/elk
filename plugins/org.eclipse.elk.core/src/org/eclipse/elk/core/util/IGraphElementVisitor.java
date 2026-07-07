/*******************************************************************************
 * Copyright (c) 2015 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.core.util;

import org.eclipse.elk.graph.ElkGraphElement;

/**
 * A graph element visitor is applied to all elements of a graph.
 */
@FunctionalInterface
public interface IGraphElementVisitor {

    /**
     * Visit the given graph element.
     */
    void visit(ElkGraphElement element);

}
