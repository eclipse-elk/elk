/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.disco;

import org.eclipse.elk.alg.disco.graph.DCGraph;
import org.eclipse.elk.alg.disco.options.DisCoMetaDataProvider;
import org.eclipse.elk.alg.disco.options.DisCoOptions;
import org.eclipse.elk.alg.disco.transform.ElkGraphTransformer;
import org.eclipse.elk.core.AbstractLayoutProvider;
import org.eclipse.elk.core.data.LayoutAlgorithmData;
import org.eclipse.elk.core.data.LayoutMetaDataService;
import org.eclipse.elk.core.util.IElkProgressMonitor;
import org.eclipse.elk.graph.ElkNode;

/**
 * Layout provider for the layout algorithms for positioning disconnected components relative to each other.
 */
public class DisCoLayoutProvider extends AbstractLayoutProvider {

    ///////////////////////////////////////////////////////////////////////////////
    // Variables

    /** simplified version of the input graph consisting of polygons approximating the original graph's shape. */
    private DCGraph result;

    /** minimum spacing between each component. */
    private double componentSpacing;

    ///////////////////////////////////////////////////////////////////////////////
    // Layout

    /*
     * (non-Javadoc)
     */
    @Override
    public void layout(final ElkNode layoutGraph, final IElkProgressMonitor progressMonitor) {
        progressMonitor.begin("Connected Components Compaction", 1);

        componentSpacing = layoutGraph.getProperty(DisCoOptions.SPACING_COMPONENT_COMPONENT);

        // If desired, apply a layout algorithm to the connected components themselves
        if (layoutGraph.hasProperty(DisCoOptions.COMPONENT_COMPACTION_COMPONENT_LAYOUT_ALGORITHM)) {
            String requestedAlgorithm =
                    layoutGraph.getProperty(DisCoOptions.COMPONENT_COMPACTION_COMPONENT_LAYOUT_ALGORITHM);
            LayoutAlgorithmData lad = LayoutMetaDataService.getInstance().getAlgorithmDataBySuffix(requestedAlgorithm);
            if (lad != null) {
                AbstractLayoutProvider layoutProvider = lad.getInstancePool().fetch();
                layoutProvider.layout(layoutGraph, progressMonitor.subTask(1));
            }
        }
        
        // 1.) Transform given KGraph into DCGraph
        ElkGraphTransformer transformer = new ElkGraphTransformer(componentSpacing);
        result = transformer.importGraph(layoutGraph);

        // 2.) Choose strategy and compact the DCGraph (only polyomino compaction at the moment)
        switch (layoutGraph.getProperty(DisCoMetaDataProvider.COMPONENT_COMPACTION_STRATEGY)) {

        // The implementation based on Freivalds et al. (2002), uses lower resolution polyominoes of the DCGraph's
        // components
        case POLYOMINO:
            new DisCoPolyominoCompactor().compact(result);
            // Only for debugging purposes, see org.eclipse.elk.alg.disco.debug.views.DisCoGraphRenderer
            layoutGraph.setProperty(DisCoOptions.DEBUG_DISCO_POLYS, result.getProperty(DisCoOptions.DEBUG_DISCO_POLYS));
            break;
        default:
            System.out.println("DisCo: no compaction strategy used for connected components.");
        }

        // 3.) Apply new layout to input graph
        transformer.applyLayout();

        // Only for debugging purposes, see org.eclipse.elk.alg.disco.debug.views.DisCoGraphRenderer
        layoutGraph.setProperty(DisCoOptions.DEBUG_DISCO_GRAPH, result);
        
        progressMonitor.done();
    }

}
