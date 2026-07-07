/*******************************************************************************
 * Copyright (c) 2018, 2019 Kiel University and others.
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
package org.eclipse.elk.alg.test.tests;

import static org.junit.Assert.assertTrue;

import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.layered.options.LayeringStrategy;
import org.eclipse.elk.alg.layered.p2layers.NetworkSimplexLayerer;
import org.eclipse.elk.alg.test.framework.LayoutTestRunner;
import org.eclipse.elk.alg.test.framework.annotations.Algorithm;
import org.eclipse.elk.alg.test.framework.annotations.Configurator;
import org.eclipse.elk.alg.test.framework.annotations.FailIfNotExecuted;
import org.eclipse.elk.alg.test.framework.annotations.GraphProvider;
import org.eclipse.elk.alg.test.framework.annotations.TestAfterProcessor;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.options.HierarchyHandling;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.util.ElkGraphUtil;
import org.junit.runner.RunWith;

/**
 * A few simple whitebox test examples.
 */
@RunWith(LayoutTestRunner.class)
@Algorithm(LayeredOptions.ALGORITHM_ID)
public class WhiteBoxTest {

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Sources

    /**
     * Returns a concrete graph to test with.
     */
    @GraphProvider
    public ElkNode compoundGraph() {
        // Create a graph with a compound node in it
        ElkNode topLevelGraph = ElkGraphUtil.createGraph();
        ElkNode compoundChild = ElkGraphUtil.createNode(topLevelGraph);
        ElkGraphUtil.createNode(compoundChild);

        return topLevelGraph;
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Configuration

    /**
     * Ensures ELK Layered will use the network simplex layerer.
     */
    @Configurator
    public void configure(final ElkNode graph) {
        graph.setProperty(CoreOptions.HIERARCHY_HANDLING, HierarchyHandling.INCLUDE_CHILDREN);
        graph.setProperty(LayeredOptions.LAYERING_STRATEGY, LayeringStrategy.NETWORK_SIMPLEX);
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Tests

    /**
     * Checks that there are no layerless nodes left.
     */
    @TestAfterProcessor(NetworkSimplexLayerer.class)
    @FailIfNotExecuted
    public void testNoLayerlessNodes(final Object graph) {
        assertTrue("There are layerless nodes left!", ((LGraph) graph).getLayerlessNodes().isEmpty());
    }

    /**
     * Checks that the layering is proper.
     */
    @TestAfterProcessor(NetworkSimplexLayerer.class)
    @FailIfNotExecuted
    public void testProperLayering(final Object graph) {
        LGraph lGraph = (LGraph) graph;
        for (Layer layer : lGraph) {
            int sourceLayerIndex = layer.getIndex();

            for (LNode lnode : layer) {
                for (LEdge ledge : lnode.getOutgoingEdges()) {
                    int targetLayerIndex = ledge.getTarget().getNode().getLayer().getIndex();
                    assertTrue("Edge points leftwards!", sourceLayerIndex <= targetLayerIndex);
                }
            }
        }
    }

}
