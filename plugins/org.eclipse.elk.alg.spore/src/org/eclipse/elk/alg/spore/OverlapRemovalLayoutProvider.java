/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.spore;

import java.util.List;
import java.util.Set;

import org.eclipse.elk.alg.common.TEdge;
import org.eclipse.elk.alg.common.spore.IOverlapHandler;
import org.eclipse.elk.alg.common.spore.InternalProperties;
import org.eclipse.elk.alg.common.spore.ScanlineOverlapCheck;
import org.eclipse.elk.alg.common.utils.SVGImage;
import org.eclipse.elk.alg.spore.graph.Graph;
import org.eclipse.elk.alg.spore.options.SporeOverlapRemovalOptions;
import org.eclipse.elk.alg.spore.options.OverlapRemovalStrategy;
import org.eclipse.elk.alg.spore.options.RootSelection;
import org.eclipse.elk.alg.spore.options.SporeCompactionOptions;
import org.eclipse.elk.alg.spore.options.SpanningTreeCostFunction;
import org.eclipse.elk.alg.spore.options.StructureExtractionStrategy;
import org.eclipse.elk.alg.spore.options.TreeConstructionStrategy;
import org.eclipse.elk.core.AbstractLayoutProvider;
import org.eclipse.elk.core.alg.AlgorithmAssembler;
import org.eclipse.elk.core.alg.ILayoutProcessor;
import org.eclipse.elk.core.data.LayoutAlgorithmData;
import org.eclipse.elk.core.data.LayoutMetaDataService;
import org.eclipse.elk.core.util.ElkUtil;
import org.eclipse.elk.core.util.IElkProgressMonitor;
import org.eclipse.elk.graph.ElkNode;

import com.google.common.collect.Sets;

/**
 * This layout provider performs a node overlap removal procedure proposed by Nachmanson et al.
 * in the paper "Node Overlap Removal by Growing a Tree" (2016).
 * This implementation of said approach lacks the scanline overlap check and therefore only removes overlap
 * along edges of the minimum spanning tree.
 */
public class OverlapRemovalLayoutProvider extends AbstractLayoutProvider {
    
    private AlgorithmAssembler<SPOrEPhases, Graph> algorithmAssembler =
            AlgorithmAssembler.<SPOrEPhases, Graph>create(SPOrEPhases.class);
    private List<ILayoutProcessor<Graph>> algorithm;

    @Override
    public void layout(final ElkNode layoutGraph, final IElkProgressMonitor progressMonitor) {
        // If desired, apply a layout algorithm
        if (layoutGraph.hasProperty(SporeCompactionOptions.UNDERLYING_LAYOUT_ALGORITHM)) {
            String requestedAlgorithm = layoutGraph.getProperty(SporeOverlapRemovalOptions.UNDERLYING_LAYOUT_ALGORITHM);
            LayoutAlgorithmData lad = LayoutMetaDataService.getInstance().getAlgorithmDataBySuffix(requestedAlgorithm);
            if (lad != null) {
                AbstractLayoutProvider layoutProvider = lad.getInstancePool().fetch();
                layoutProvider.layout(layoutGraph, progressMonitor.subTask(1));
            }
        }
        
        // set algorithm properties
        layoutGraph.setProperty(SporeCompactionOptions.PROCESSING_ORDER_ROOT_SELECTION, RootSelection.CENTER_NODE);
        layoutGraph.setProperty(SporeCompactionOptions.PROCESSING_ORDER_SPANNING_TREE_COST_FUNCTION, 
                SpanningTreeCostFunction.INVERTED_OVERLAP);
        layoutGraph.setProperty(SporeCompactionOptions.PROCESSING_ORDER_TREE_CONSTRUCTION, 
                TreeConstructionStrategy.MINIMUM_SPANNING_TREE);
        int maxIterations = layoutGraph.getProperty(SporeOverlapRemovalOptions.OVERLAP_REMOVAL_MAX_ITERATIONS);
        
        progressMonitor.begin("Overlap removal", 1);
        
        // initialize debug output
        String debugOutputFile = null;
        if (layoutGraph.getProperty(SporeOverlapRemovalOptions.DEBUG_MODE)) {
            debugOutputFile = ElkUtil.debugFolderPath("spore") + "45scanlineOverlaps";
        }
        SVGImage svg = new SVGImage(debugOutputFile);
        
        // set overlap handler and import ElkGraph
        Set<TEdge> overlapEdges = Sets.newHashSet();
        IOverlapHandler overlapHandler = (n1, n2) -> overlapEdges.add(new TEdge(n1.originalVertex, n2.originalVertex));
        IGraphImporter<ElkNode> graphImporter = new ElkGraphImporter();
        Graph graph = graphImporter.importGraph(layoutGraph);
            
        boolean overlapsExisted = true;
        int iteration = 0;
        
        // repeat overlap removal
        while (iteration < maxIterations && overlapsExisted) {
            
            // scanline overlap check
            if (layoutGraph.getProperty(SporeOverlapRemovalOptions.OVERLAP_REMOVAL_RUN_SCANLINE)) {
                overlapEdges.clear();
                new ScanlineOverlapCheck(overlapHandler, svg).sweep(graph.vertices);
                if (overlapEdges.isEmpty()) {
                    break; // don't bother if nothing overlaps
                }
                graph.tEdges = overlapEdges;
            }
            
            // assembling and executing the algorithm
            algorithmAssembler.reset();
            algorithmAssembler.setPhase(SPOrEPhases.P1_STRUCTURE, 
                    StructureExtractionStrategy.DELAUNAY_TRIANGULATION);
            algorithmAssembler.setPhase(SPOrEPhases.P2_PROCESSING_ORDER, 
                    graph.treeConstructionStrategy);
            algorithmAssembler.setPhase(SPOrEPhases.P3_EXECUTION, 
                    OverlapRemovalStrategy.GROW_TREE);
            algorithm = algorithmAssembler.build(graph);
            for (ILayoutProcessor<Graph> processor : algorithm) {
                processor.process(graph, progressMonitor.subTask(1));
            }
            
            // update node positions 
            graphImporter.updateGraph(graph);
            
            overlapsExisted = graph.getProperty(InternalProperties.OVERLAPS_EXISTED);
            iteration++;
        }

        // apply node positions to ElkGraph
        graphImporter.applyPositions(graph);

        progressMonitor.done();
    }

}
