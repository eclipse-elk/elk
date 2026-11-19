/*******************************************************************************
 * Copyright (c) 2026 Kiel University and others.
 *
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 ******************************************************************************/
package org.eclipse.elk.alg.layered.intermediate.wrapping;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;

import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.intermediate.greedyswitch.TestGraphCreator;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.layered.options.WrappingStrategy;
import org.junit.Test;

/** Tests that wrapping heuristics avoid cuts next to edge label dummies. */
public class GraphWrappingCutIndexTest extends TestGraphCreator {

    @Test
    public void testArdAvoidsEdgeLabelCut() {
        assertHeuristicAvoidsEdgeLabelCut(new ARDCutIndexHeuristic(), WrappingStrategy.SINGLE_EDGE, "ARD");
    }

    @Test
    public void testArdAvoidsEdgeLabelCutForMultiEdgeWrapping() {
        assertHeuristicAvoidsEdgeLabelCut(new ARDCutIndexHeuristic(), WrappingStrategy.MULTI_EDGE, "ARD");
    }

    @Test
    public void testMsdAvoidsEdgeLabelCut() {
        assertHeuristicAvoidsEdgeLabelCut(new MSDCutIndexHeuristic(), WrappingStrategy.SINGLE_EDGE, "MSD");
    }

    @Test
    public void testMsdAvoidsEdgeLabelCutForMultiEdgeWrapping() {
        assertHeuristicAvoidsEdgeLabelCut(new MSDCutIndexHeuristic(), WrappingStrategy.MULTI_EDGE, "MSD");
    }

    private void assertHeuristicAvoidsEdgeLabelCut(final ICutIndexCalculator heuristic,
            final WrappingStrategy strategy, final String heuristicName) {
        LGraph graph = graphWithLabel(8, strategy);
        GraphStats stats = new GraphStats(graph);

        List<Integer> cuts = heuristic.getCutIndexes(graph, stats);

        assertTrue(cuts.size() > 0);
        for (int cut : cuts) {
            assertTrue(heuristicName + " returned forbidden cut " + cut + " for " + strategy,
                    stats.isCutAllowedForWrapping(cut));
        }
    }

    @Test
    public void testCutsAdjacentToLabelAreForbidden() {
        LGraph graph = graphWithLabel(5, WrappingStrategy.SINGLE_EDGE);
        GraphStats stats = new GraphStats(graph);

        assertFalse(stats.isCutAllowedForWrapping(2));
        assertFalse(stats.isCutAllowedForWrapping(3));
    }

    @Test
    public void testMultiEdgeCutsAdjacentToLabelAreForbidden() {
        LGraph graph = graphWithLabel(5, WrappingStrategy.MULTI_EDGE);
        GraphStats stats = new GraphStats(graph);

        assertFalse(stats.isCutAllowedForWrapping(2));
        assertFalse(stats.isCutAllowedForWrapping(3));
    }

    private LGraph graphWithLabel(final int layerCount, final WrappingStrategy strategy) {
        graph = new LGraph();
        Layer[] layers = makeLayers(layerCount);
        for (Layer layer : layers) {
            LNode node = addNodeToLayer(layer);
            node.getSize().x = 10;
            node.getSize().y = 10;
        }
        LNode label = addNodeToLayer(layers[2]);
        label.setType(LNode.NodeType.LABEL);
        label.getSize().x = 2;
        label.getSize().y = 2;
        graph.setProperty(LayeredOptions.ASPECT_RATIO, 0.1);
        graph.setProperty(LayeredOptions.WRAPPING_STRATEGY, strategy);
        return graph;
    }
}
