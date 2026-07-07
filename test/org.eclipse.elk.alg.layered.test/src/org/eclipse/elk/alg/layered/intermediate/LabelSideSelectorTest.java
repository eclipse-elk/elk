/*******************************************************************************
 * Copyright (c) 2013, 2020 Kiel University and others.
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
package org.eclipse.elk.alg.layered.intermediate;

import static org.junit.Assert.assertTrue;

import java.util.Deque;
import java.util.LinkedList;
import java.util.List;
import java.util.Random;

import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.test.framework.LayoutTestRunner;
import org.eclipse.elk.alg.test.framework.annotations.Algorithm;
import org.eclipse.elk.alg.test.framework.annotations.Configurator;
import org.eclipse.elk.alg.test.framework.annotations.DefaultConfiguration;
import org.eclipse.elk.alg.test.framework.annotations.GraphResourceProvider;
import org.eclipse.elk.alg.test.framework.annotations.TestAfterProcessor;
import org.eclipse.elk.alg.test.framework.io.AbstractResourcePath;
import org.eclipse.elk.alg.test.framework.io.FileExtensionFilter;
import org.eclipse.elk.alg.test.framework.io.ModelResourcePath;
import org.eclipse.elk.core.options.EdgeLabelPlacement;
import org.eclipse.elk.core.options.LabelSide;
import org.eclipse.elk.graph.ElkEdge;
import org.eclipse.elk.graph.ElkLabel;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.util.ElkGraphUtil;
import org.junit.runner.RunWith;

import com.google.common.collect.Lists;
import com.google.common.collect.Streams;

@RunWith(LayoutTestRunner.class)
@Algorithm(LayeredOptions.ALGORITHM_ID)
@DefaultConfiguration()
public class LabelSideSelectorTest {

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Sources

    @GraphResourceProvider
    public List<AbstractResourcePath> testGraphs() {
        return Lists.newArrayList(
                new ModelResourcePath("realworld/ptolemy/**/").withFilter(new FileExtensionFilter("elkg")),
                new ModelResourcePath("tickets/layered/599_unimportantPortIndex.elkt"));
    }
    

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Configuration
    
    /** Fixed random seed to keep tests reproducible. */
    private final int RANDOM_SEED = 0;
    
    @Configurator
    public void configurator(final ElkNode graph) {
        Deque<ElkNode> nodeQueue = new LinkedList<>();
        nodeQueue.add(graph);
        
        Random rand = new Random(RANDOM_SEED);
        
        while (!nodeQueue.isEmpty()) {
            ElkNode node = nodeQueue.poll();
            
            for (ElkEdge edge : node.getContainedEdges()) {
                if (edge.getLabels().isEmpty()) {
                    ElkGraphUtil.createLabel(edge.getIdentifier(), edge);
                }
                
                ElkLabel label = edge.getLabels().get(0);
                
                double r = rand.nextDouble();
                if (r < 0.3) {
                    label.setProperty(LayeredOptions.EDGE_LABELS_PLACEMENT, EdgeLabelPlacement.HEAD);
                } else if (r < 0.7) {
                    label.setProperty(LayeredOptions.EDGE_LABELS_PLACEMENT, EdgeLabelPlacement.CENTER);
                } else {
                    label.setProperty(LayeredOptions.EDGE_LABELS_PLACEMENT, EdgeLabelPlacement.TAIL);
                }
            }
        }
    }


    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Tests

    /**
     * All labels on ports and edges have an assigned {@link LabelSide}.
     */
    @TestAfterProcessor(LabelSideSelector.class)
    public void testRemovedNodes(final Object graph) {
        LGraph lGraph = (LGraph) graph;
        
        lGraph.getLayers().stream()
            .flatMap(layer -> layer.getNodes().stream())
            .flatMap(node -> Streams.stream(node.getOutgoingEdges()))
            .flatMap(edge -> edge.getLabels().stream())
            .forEach(label -> assertTrue(label.getProperty(InternalProperties.LABEL_SIDE) != LabelSide.UNKNOWN));
    }

}
