/*******************************************************************************
 * Copyright (c) 2013 - 2022 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.mrtree;

import java.util.List;

import org.eclipse.elk.alg.common.NodeMicroLayout;
import org.eclipse.elk.alg.mrtree.graph.TGraph;
import org.eclipse.elk.alg.mrtree.options.MrTreeOptions;
import org.eclipse.elk.core.AbstractLayoutProvider;
import org.eclipse.elk.core.util.IElkProgressMonitor;
import org.eclipse.elk.graph.ElkNode;

/**
 * Layout provider to connect the tree layouter to the Eclipse based layout services and orchestrate
 * the pre layout processing.
 * 
 * @author sor
 * @author sgu
 * @author sdo
 */
public class TreeLayoutProvider extends AbstractLayoutProvider {

    // /////////////////////////////////////////////////////////////////////////////
    // Variables

    /** the layout algorithm used for this layout. */
    private MrTree klayTree = new MrTree();
    /** connected components processor. */
    private ComponentsProcessor componentsProcessor = new ComponentsProcessor();
    
    private final float defaultWork = 0.1f;

    // /////////////////////////////////////////////////////////////////////////////
    // Regular Layout

    /**
     * {@inheritDoc}
     */
    @Override
    public void layout(final ElkNode layoutGraph, final IElkProgressMonitor progressMonitor) {
        
        // If requested, compute nodes's dimensions, place node labels, ports, port labels, etc.
        if (!layoutGraph.getProperty(MrTreeOptions.OMIT_NODE_MICRO_LAYOUT)) {
            NodeMicroLayout.forGraph(layoutGraph)
                           .execute();
        }
        // build tGraph
        IElkProgressMonitor pm = progressMonitor.subTask(defaultWork);
        pm.begin("build tGraph", 1);
        IGraphImporter<ElkNode> graphImporter = new ElkGraphImporter();
        TGraph tGraph = graphImporter.importGraph(layoutGraph);
        pm.done();

        // split the input graph into components
        pm = progressMonitor.subTask(defaultWork);
        pm.begin("Split graph", 1);
        List<TGraph> components = componentsProcessor.split(tGraph);
        pm.done();

        // perform the actual layout on the components
        for (TGraph comp : components) {
            klayTree.doLayout(comp, progressMonitor.subTask((1.0f - defaultWork * 4) / components.size()));
        }

        // pack the components back into one graph
        pm = progressMonitor.subTask(defaultWork);
        pm.begin("Pack components", 1);
        tGraph = componentsProcessor.pack(components);
        pm.done();

        // apply the layout results to the original graph
        pm = progressMonitor.subTask(defaultWork);
        pm.begin("Apply layout results", 1);
        graphImporter.applyLayout(tGraph);
        pm.done();
    }

}
