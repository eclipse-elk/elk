/*******************************************************************************
 * Copyright (c) 2013, 2015 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.

 * This Source Code may also be made available under the following Secondary
 * Licenses when the conditions for such availability set forth in the Eclipse
 * Public License v. 2.0 are satisfied: GPL-3.0 which is available at
 * https://www.gnu.org/licenses/gpl-3.0-standalone.html.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.mrtree;

import org.eclipse.elk.alg.mrtree.graph.TGraph;

/**
 * Interface for importer classes for graph structure. Graph importers
 * should usually subclass
 * {@link org.eclipse.elk.alg.mrtree.ElkGraphImporter KGraphImporter}
 * instead of implementing this interface directly.
 * 
 * <p>Graph importers are encouraged to set the {@link MrTreeOptions#GRAPH_PROPERTIES}
 * property on imported graphs.</p>
 *
 * @param <T> the type of graph that this importer can transform into a layered graph.
 * @author sor
 * @author sgu
 * @author msp
 */
public interface IGraphImporter<T> {
    
    /**
     * Create a t-graph from the given graph.
     * 
     * @param graph the graph to turn into a t-graph
     * @return a t-graph, or {@code null} if the input was not recognized
     */
    TGraph importGraph(T graph);
    
    /**
     * Apply the computed layout of a t-graph to the original graph.
     * 
     * @param tGraph the graph for which layout is applied
     */
    void applyLayout(TGraph tGraph);

}
