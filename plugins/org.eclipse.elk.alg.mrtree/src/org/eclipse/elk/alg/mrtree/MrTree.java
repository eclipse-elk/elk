/*******************************************************************************
 * Copyright (c) 2013 - 2022 Kiel University and others.
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

import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.List;

import org.eclipse.elk.alg.mrtree.graph.TGraph;
import org.eclipse.elk.alg.mrtree.options.MrTreeOptions;
import org.eclipse.elk.core.alg.AlgorithmAssembler;
import org.eclipse.elk.core.alg.ILayoutProcessor;
import org.eclipse.elk.core.util.BasicProgressMonitor;
import org.eclipse.elk.core.util.IElkProgressMonitor;

/**
 * Implements a layout algorithm for trees. The tree layouter uses the algorithm from
 * "A Node-Positioning Algorithm for General Trees, John Q.Walker II" to layout trees. To do this it
 * uses four phases plus a pre-processing to build a corresponding data structure. The first phase
 * "treeifying" transforms the given graph into a tree if necessary. To do this, edges which destroy
 * the tree property will be removed and stored, so that they can be reinserted during a post
 * processing. In the second phase "orderNodes" the nodes of each level are separated into leaves
 * and inner nodes. And then whitespace in the level is filled with leaves. The third phase
 * "NodePlacer" uses the algorithm first mentioned from John Q.Walker II to compute the actual
 * position of the nodes. The last phase routeEdges sets the positions for the edges corresponding
 * to the positions of the nodes.
 * 
 * Each phase uses intermediate processors for small computations on the graph. The corresponding
 * processors are defined in each phase. Some are defined multiple times, but they are invoked only
 * once.
 * 
 * The processors can determine roots of components and fan outs or level neighbors of nodes.
 * 
 * @author sor
 * @author sgu
 * @author sdo
 */
public final class MrTree {

    // /////////////////////////////////////////////////////////////////////////////
    // Variables
    
    /** The thing that will assemble the processors that will constitute our algorithm. */
    private AlgorithmAssembler<TreeLayoutPhases, TGraph> algorithmAssembler =
            AlgorithmAssembler.<TreeLayoutPhases, TGraph>create(TreeLayoutPhases.class);
    
    /** Our algorithm. */
    private List<ILayoutProcessor<TGraph>> algorithm;

    // /////////////////////////////////////////////////////////////////////////////
    // Regular Layout

    /**
     * Does a layout on the given graph.
     * 
     * @param tgraph
     *            the graph to layout.
     * @param progressMonitor
     *            a progress monitor to show progress information in.
     * @return tree graph with layout applied.
     */
    public TGraph doLayout(final TGraph tgraph, final IElkProgressMonitor progressMonitor) {
        progressMonitor.begin("Tree layout", 1);

        // elkjs-exclude-start
        if (tgraph.getProperty(MrTreeOptions.DEBUG_MODE)) {
            progressMonitor.log("MrTree! called at "
                    + LocalDateTime.now().format(DateTimeFormatter.ofPattern("HH:mm:ss")));
        }
        // elkjs-exclude-end
        
        // set up the phases and processors depending on user options
        updateModules(tgraph);

        // do layout for each component
        layout(tgraph, progressMonitor.subTask(1.0f));

        progressMonitor.done();

        return tgraph;
    }

    // /////////////////////////////////////////////////////////////////////////////
    // Options and Modules Management

    /**
     * Update the modules depending on user options.
     * 
     * @param graph
     *            the graph to be laid out.
     */
    private void updateModules(final TGraph graph) {
        algorithmAssembler.reset();
        
        // Configure phases
        algorithmAssembler.setPhase(TreeLayoutPhases.P1_TREEIFICATION, TreeLayoutPhases.P1_TREEIFICATION);
        algorithmAssembler.setPhase(TreeLayoutPhases.P2_NODE_ORDERING, TreeLayoutPhases.P2_NODE_ORDERING);
        algorithmAssembler.setPhase(TreeLayoutPhases.P3_NODE_PLACEMENT, TreeLayoutPhases.P3_NODE_PLACEMENT);
        algorithmAssembler.setPhase(TreeLayoutPhases.P4_EDGE_ROUTING, TreeLayoutPhases.P4_EDGE_ROUTING);
        
        // Assemble the algorithm
        algorithm = algorithmAssembler.build(graph);
    }

    // TODO add this for intermediate processing configuration for block direction
    // /**
    // * Returns an intermediate processing configuration with processors not tied to specific
    // phases.
    // *
    // * @param graph
    // * the tree graph to be processed. The configuration may vary depending on certain
    // * properties of the graph.
    // * @return intermediate processing configuration. May be {@code null}.
    // */
    // private IntermediateProcessingConfiguration getIntermediateProcessingConfiguration(
    // final TGraph graph) {
    //
    // // get graph properties
    //
    // // Basic configuration
    // IntermediateProcessingConfiguration configuration = new
    // IntermediateProcessingConfiguration();
    //
    // // graph transformations for unusual layout directions
    // switch (graph.getProperty(LayoutOptions.DIRECTION)) {
    // case LEFT:
    // // ADD LEFT_DIR_POSTPROCESSOR to IntermediateProcessingConfiguration
    // break;
    // case RIGHT:
    // // ADD RIGHT_DIR_POSTPROCESSOR to IntermediateProcessingConfiguration
    // break;
    // case DOWN:
    // // ADD DOWN_DIR_POSTPROCESSOR to IntermediateProcessingConfiguration
    // break;
    // case UP:
    // // ADD UP_DIR_POSTPROCESSOR to IntermediateProcessingConfiguration
    // break;
    // default:
    // // This is either RIGHT or UNDEFINED, which is just mapped to RIGHT. Either way, we
    // // don't
    // // need any processors here
    // break;
    // }
    //
    // // Additional dependencies
    //
    // return configuration;
    // }

    // /////////////////////////////////////////////////////////////////////////////
    // Layout

    /**
     * Perform the phases of the tree layouter.
     * 
     * @param graph
     *            the graph that is to be laid out
     * @param themonitor
     *            a progress monitor, or {@code null}
     */
    private void layout(final TGraph graph, final IElkProgressMonitor themonitor) {
        IElkProgressMonitor monitor = themonitor;
        if (monitor == null) {
            monitor = new BasicProgressMonitor();
        }
        monitor.begin("Layout", algorithm.size());

        // elkjs-exclude-start
        if (graph.getProperty(MrTreeOptions.DEBUG_MODE)) {
            // Debug Mode!
            // Prints the algorithm configuration and outputs the whole graph to a file
            // before each slot execution
            themonitor.log("ELK MrTree uses the following " + algorithm.size() + " modules:");
            for (int i = 0; i < algorithm.size(); i++) {
                themonitor.log("   Slot " + i + ": " + algorithm.get(i).getClass().getName());
            }
        }
        // elkjs-exclude-end
        // invoke each layout processor
        for (ILayoutProcessor<TGraph> processor : algorithm) {
            if (monitor.isCanceled()) {
                return;
            }
            processor.process(graph, monitor.subTask(1));
        }
        monitor.done();
    }
}
