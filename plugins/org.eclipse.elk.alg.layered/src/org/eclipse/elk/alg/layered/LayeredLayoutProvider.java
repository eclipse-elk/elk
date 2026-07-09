/*******************************************************************************
 * Copyright (c) 2010, 2015 Kiel University and others.
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
package org.eclipse.elk.alg.layered;

import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.transform.ElkGraphTransformer;
import org.eclipse.elk.alg.layered.graph.transform.IGraphTransformer;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.core.AbstractLayoutProvider;
import org.eclipse.elk.core.options.HierarchyHandling;
import org.eclipse.elk.core.testing.IWhiteBoxTestable;
import org.eclipse.elk.core.testing.TestController;
import org.eclipse.elk.core.util.IElkProgressMonitor;
import org.eclipse.elk.graph.ElkNode;

/**
 * Layout provider to connect the layered layouter to the Eclipse based layout services.
 * 
 * @see ElkLayered
 * 
 * @author msp
 * @author cds
 */
public final class LayeredLayoutProvider extends AbstractLayoutProvider implements IWhiteBoxTestable {

    ///////////////////////////////////////////////////////////////////////////////
    // Variables

    /** the layout algorithm used for regular layout runs. */
    private final ElkLayered elkLayered = new ElkLayered();


    ///////////////////////////////////////////////////////////////////////////////
    // Regular Layout

    @Override
    public void layout(final ElkNode elkgraph, final IElkProgressMonitor progressMonitor) {
        // Import the graph (layeredGraph won't be null since the ElkGraphTransformer always returns an LGraph
        // instance, even though the IGraphTransformer interface would allow null as a return value)
        IGraphTransformer<ElkNode> graphTransformer = new ElkGraphTransformer();
        LGraph layeredGraph = graphTransformer.importGraph(elkgraph);

        // Check if hierarchy handling for a compound graph is requested
        if (elkgraph.getProperty(LayeredOptions.HIERARCHY_HANDLING) == HierarchyHandling.INCLUDE_CHILDREN) {
            // Layout for all hierarchy levels is requested
            elkLayered.doCompoundLayout(layeredGraph, progressMonitor);
        } else {
            // Only the top-level graph is processed
            elkLayered.doLayout(layeredGraph, progressMonitor);
        }
        
        if (!progressMonitor.isCanceled()) {
            // Apply the layout results to the original graph
            graphTransformer.applyLayout(layeredGraph);
        }
    }


    ///////////////////////////////////////////////////////////////////////////////
    // Layout Testing
    
    /**
     * Import the given KGraph and return a test execution state prepared for a test run with the
     * resulting {@link LGraph}. The layout test run methods can immediately be called on the
     * returned object.
     * 
     * <p><strong>Note:</strong> This method does not apply the layout back to the original KGraph!</p>
     * 
     * @param elkgraph the KGraph to be used for the layout test run.
     * @return an initialized test execution state
     */
    public ElkLayered.TestExecutionState startLayoutTest(final ElkNode elkgraph) {
        // Import the graph (layeredGraph won't be null since the KGraphImporter always returns an
        // LGraph instance, even though the IGraphImporter interface would allow null as a return
        // value)
        IGraphTransformer<ElkNode> graphImporter = new ElkGraphTransformer();
        LGraph layeredGraph = graphImporter.importGraph(elkgraph);
        
        // Prepare a layout test and return the test execution state
        return elkLayered.prepareLayoutTest(layeredGraph);
    }
    
    /**
     * Return the layered layout algorithm.
     * 
     * @return the layout algorithm
     */
    public ElkLayered getLayoutAlgorithm() {
        return elkLayered;
    }

    @Override
    public void setTestController(final TestController controller) {
        elkLayered.setTestController(controller);
    }

}
