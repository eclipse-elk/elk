/*******************************************************************************
 * Copyright (c) 2022 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.mrtree.intermediate;

import org.eclipse.elk.alg.mrtree.graph.TGraph;
import org.eclipse.elk.alg.mrtree.options.InternalProperties;
import org.eclipse.elk.core.alg.ILayoutProcessor;
import org.eclipse.elk.core.util.IElkProgressMonitor;

/**
 * Sets the graphs x/y max/min properties.
 * 
 * @author jnc, sdo
 */
public class GraphBoundsProcessor implements ILayoutProcessor<TGraph> {
    /**
     * {@inheritDoc}
     */
    public void process(final TGraph tGraph, final IElkProgressMonitor progressMonitor) {
        progressMonitor.begin("Process graph bounds", 1);
        
        tGraph.setProperty(InternalProperties.GRAPH_XMIN, 
                tGraph.getNodes().stream().mapToDouble(x -> x.getPosition().x).min().getAsDouble());
        tGraph.setProperty(InternalProperties.GRAPH_YMIN, 
                tGraph.getNodes().stream().mapToDouble(x -> x.getPosition().y).min().getAsDouble());
        tGraph.setProperty(InternalProperties.GRAPH_XMAX, 
                tGraph.getNodes().stream().mapToDouble(x -> x.getPosition().x + x.getSize().x).max().getAsDouble());
        tGraph.setProperty(InternalProperties.GRAPH_YMAX, 
                tGraph.getNodes().stream().mapToDouble(x -> x.getPosition().y + x.getSize().y).max().getAsDouble());
        
        progressMonitor.done();
    }
}
