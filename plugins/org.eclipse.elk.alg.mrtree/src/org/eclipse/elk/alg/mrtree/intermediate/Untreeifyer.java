/*******************************************************************************
 * Copyright (c) 2013 - 2022 University and others.
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
package org.eclipse.elk.alg.mrtree.intermediate;

import java.util.List;

import org.eclipse.elk.alg.mrtree.graph.TEdge;
import org.eclipse.elk.alg.mrtree.graph.TGraph;
import org.eclipse.elk.alg.mrtree.options.InternalProperties;
import org.eclipse.elk.core.alg.ILayoutProcessor;
import org.eclipse.elk.core.util.IElkProgressMonitor;

/**
 * This processor should run as a post-processor after the first phase. It should run through 
 * the list of edges that once destroyed the tree property and simply insert that edges directly after
 * layouting the graph.
 * 
 * @author sor
 * @author sgu
 */
public class Untreeifyer implements ILayoutProcessor<TGraph> {

    /**
     * {@inheritDoc}
     */
    public void process(final TGraph tGraph, final IElkProgressMonitor progressMonitor) {
        progressMonitor.begin("Untreeify", 1);
        
        List<TEdge> edges = tGraph.getProperty(InternalProperties.REMOVABLE_EDGES);

        for (TEdge tEdge : edges) {
            tEdge.getSource().getOutgoingEdges().add(tEdge);
            tEdge.getTarget().getIncomingEdges().add(tEdge);
        }
        
        progressMonitor.done();
    }
}
