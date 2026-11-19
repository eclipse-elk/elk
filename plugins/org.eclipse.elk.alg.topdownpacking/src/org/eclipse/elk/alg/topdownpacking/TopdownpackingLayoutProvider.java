/*******************************************************************************
 * Copyright (c) 2022 - 2023 Kiel University and others.
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
package org.eclipse.elk.alg.topdownpacking;

import java.util.List;

import org.eclipse.elk.alg.topdownpacking.options.TopdownpackingOptions;
import org.eclipse.elk.core.AbstractLayoutProvider;
import org.eclipse.elk.core.ITopdownLayoutProvider;
import org.eclipse.elk.core.alg.AlgorithmAssembler;
import org.eclipse.elk.core.alg.ILayoutProcessor;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.util.IElkProgressMonitor;
import org.eclipse.elk.graph.ElkNode;

/**
 * A simple box packing algorithm that places nodes as evenly sized rectangles. This algorithm uses fixed sizes and 
 * therefore requires the option 'Topdown Layout' to be set to true to result in a correct layout. Using this algorithm
 * makes parallel nodes possible by allowing a node know its child layout size before the child layout is computed.
 * 
 */
public class TopdownpackingLayoutProvider extends AbstractLayoutProvider implements ITopdownLayoutProvider {

    private final AlgorithmAssembler<TopdownPackingPhases, GridElkNode> algorithmAssembler =
            AlgorithmAssembler.<TopdownPackingPhases, GridElkNode>create(TopdownPackingPhases.class);
    
    @Override
    public void layout(ElkNode layoutGraph, IElkProgressMonitor progressMonitor) {
        
        // convert graph to GridElkNode
        GridElkNode wrappedGraph = new GridElkNode(layoutGraph);
        
        // assemble algorithm
        List<ILayoutProcessor<GridElkNode>> algorithm = assembleAlgorithm(wrappedGraph);
        
        
        progressMonitor.begin("Topdown Packing", algorithm.size());
        
        for (ILayoutProcessor<GridElkNode> processor : algorithm) {
            processor.process(wrappedGraph, progressMonitor.subTask(1));
        }
        
        progressMonitor.done();
        
    }
    
    public List<ILayoutProcessor<GridElkNode>> assembleAlgorithm(GridElkNode graph) {
        algorithmAssembler.reset();
        
        algorithmAssembler.setPhase(TopdownPackingPhases.P1_NODE_ARRANGEMENT, 
                graph.getProperty(TopdownpackingOptions.NODE_ARRANGEMENT_STRATEGY));
        algorithmAssembler.setPhase(TopdownPackingPhases.P2_WHITESPACE_ELIMINATION, 
                graph.getProperty(TopdownpackingOptions.WHITESPACE_ELIMINATION_STRATEGY));
        
        return algorithmAssembler.build(graph);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public KVector getPredictedGraphSize(ElkNode graph) {
        // FIXME: enforce that all node placement strategies implement INodePlacer
        INodeArranger nodePlacer = (INodeArranger) graph.getProperty(TopdownpackingOptions.NODE_ARRANGEMENT_STRATEGY).create();
        return nodePlacer.getPredictedSize(graph);
    }
}
