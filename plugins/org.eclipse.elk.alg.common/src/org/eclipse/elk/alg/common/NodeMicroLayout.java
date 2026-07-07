/*******************************************************************************
 * Copyright (c) 2020 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.common;

import org.eclipse.elk.alg.common.nodespacing.NodeDimensionCalculation;
import org.eclipse.elk.core.util.adapters.ElkGraphAdapters;
import org.eclipse.elk.core.util.adapters.GraphAdapters.GraphAdapter;
import org.eclipse.elk.graph.ElkNode;

/**
 * Utility class to execute what can be considered as "node micro layout". That is, automatically computing node
 * dimensions if requested by the user, positioning ports, positioning node and port labels, etc.
 * 
 * Currently, simply all possible micro layouts are executed. It may make sense to perform only a subset, e.g. only
 * place node labels. This may be implemented in the future.
 * 
 * Furthermore, all of {@link NodeDimensionCalculation}'s methods should be moved to this utility class at some point as
 * this class's naming makes it a bit clearer what actually happens.
 */
public final class NodeMicroLayout {

    private final GraphAdapter<?> adapter;

    private NodeMicroLayout(final GraphAdapter<?> adapter) {
        this.adapter = adapter;
    }

    /**
     * @return a new micro layout instance for the passed graph. Use {@link NodeMicroLayout#execute()} to finally
     *         perform the layout process.
     */
    public static NodeMicroLayout forGraph(final ElkNode elkGraph) {
        return forGraph(ElkGraphAdapters.adapt(elkGraph));
    }

    /**
     * @return a new micro layout instance for the passed graph. Use {@link NodeMicroLayout#execute()} to finally
     *         perform the layout process.
     */
    public static NodeMicroLayout forGraph(final GraphAdapter<?> adapter) {
        NodeMicroLayout builder = new NodeMicroLayout(adapter);
        return builder;
    }

    /**
     * Perform the actual layout.
     */
    public void execute() {
        NodeDimensionCalculation.sortPortLists(adapter);
        NodeDimensionCalculation.calculateLabelAndNodeSizes(adapter);
        NodeDimensionCalculation.calculateNodeMargins(adapter);
    }

}
