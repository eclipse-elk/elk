/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.core.comments;

import java.util.Collection;
import java.util.Objects;
import java.util.stream.Collectors;

import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.util.ElkGraphUtil;

/**
 * Implementation of {@link IDataProvider} for ELK graphs. Comments are nodes that have the
 * {@link CoreOptions#COMMENT_BOX COMMENT_BOX} property set. Possible attachment targets are the siblings of a comment
 * node.
 */
public class ElkGraphDataProvider implements IDataProvider<ElkNode, ElkNode> {

    private final ElkNode graph;
    
    /**
     * Creates a new instance that provides data from the given graph.
     * 
     * @param graph
     *            non-{@code null} graph to provide data for.
     * @throws IllegalArgumentException
     *             if {@code graph == null}.
     */
    public ElkGraphDataProvider(final ElkNode graph) {
        Objects.requireNonNull(graph, "Comment attachment must be run on a graph.");
        this.graph = graph;
    }
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////
    // IDataProvider
    
    @Override
    public Collection<ElkNode> provideComments() {
        // Simply return all nodes marked as comments
        return graph.getChildren().stream()
            .filter(node -> node.getProperty(CoreOptions.COMMENT_BOX))
            .collect(Collectors.toList());
    }

    @Override
    public Collection<ElkNode> provideTargets() {
        // Return all non-comment nodes
        return graph.getChildren().stream()
                .filter(node -> !node.getProperty(CoreOptions.COMMENT_BOX))
                .collect(Collectors.toList());
    }

    @Override
    public Collection<IDataProvider<ElkNode, ElkNode>> provideSubHierarchies() {
        return graph.getChildren().stream()
            .filter(node -> !node.getChildren().isEmpty())
            .map(node -> new ElkGraphDataProvider(node))
            .collect(Collectors.toList());
    }

    /**
     * Attaches a comment to an attachment target by inserting an edge that connects the two.
     */
    @Override
    public void attach(final ElkNode comment, final ElkNode target) {
        ElkGraphUtil.createSimpleEdge(comment, target);
    }

}
