/*******************************************************************************
 * Copyright (c) 2026 Kiel University and others.
 *
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 ******************************************************************************/
package org.eclipse.elk.alg.layered.intermediate;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;

import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.intermediate.greedyswitch.TestGraphCreator;
import org.eclipse.elk.alg.layered.intermediate.wrapping.ARDCutIndexHeuristic;
import org.eclipse.elk.alg.layered.intermediate.wrapping.GraphStats;
import org.eclipse.elk.alg.layered.intermediate.wrapping.MSDCutIndexHeuristic;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.junit.Test;

/** Tests that wrapping heuristics avoid cuts next to edge label dummies. */
public class GraphWrappingCutIndexTest extends TestGraphCreator {

    @Test
    public void testArdAvoidsEdgeLabelCut() {
        LGraph graph = graphWithLabel(8);
        GraphStats stats = new GraphStats(graph);

        List<Integer> cuts = new ARDCutIndexHeuristic().getCutIndexes(graph, stats);

        assertTrue(cuts.size() > 0);
        for (int cut : cuts) {
            assertTrue("ARD returned forbidden cut " + cut, stats.isCutAllowed(cut));
        }
    }

    @Test
    public void testMsdAvoidsEdgeLabelCut() {
        LGraph graph = graphWithLabel(8);
        GraphStats stats = new GraphStats(graph);

        List<Integer> cuts = new MSDCutIndexHeuristic().getCutIndexes(graph, stats);

        assertTrue(cuts.size() > 0);
        for (int cut : cuts) {
            assertTrue("MSD returned forbidden cut " + cut, stats.isCutAllowed(cut));
        }
    }

    @Test
    public void testCutsAdjacentToLabelAreForbidden() {
        LGraph graph = graphWithLabel(5);
        GraphStats stats = new GraphStats(graph);

        assertFalse(stats.isCutAllowed(2));
        assertFalse(stats.isCutAllowed(3));
    }

    private LGraph graphWithLabel(final int layerCount) {
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
        return graph;
    }
}
