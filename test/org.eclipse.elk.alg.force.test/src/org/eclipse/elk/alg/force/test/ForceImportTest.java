/*******************************************************************************
 * Copyright (c) 2020 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 * 
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later 
 *******************************************************************************/
package org.eclipse.elk.alg.force.test;

import static org.hamcrest.Matchers.is;
import static org.junit.Assert.assertThat;

import java.util.List;

import org.eclipse.elk.alg.force.ComponentsProcessor;
import org.eclipse.elk.alg.force.ElkGraphImporter;
import org.eclipse.elk.alg.force.graph.FGraph;
import org.eclipse.elk.alg.force.options.ForceOptions;
import org.eclipse.elk.graph.ElkEdge;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.util.ElkGraphUtil;
import org.junit.Test;

/**
 * Tests related to importing an ElkGraph into an {@link FGraph}.
 */
public class ForceImportTest {

    // CHECKSTYLEOFF MagicNumber

    /**
     * Test importing a simple graph.
     */
    @Test
    public void testImport() {
        FGraph fGraph = createSimpleGraph();
        checkSimpleGraph(fGraph);
    }

    /**
     * Ensure that connected components are split when desired.
     * 
     * Also verifies that the split graphs have the correct number of nodes, edges, and labels.
     */
    @Test
    public void testSeparateConnectedComponents() {
        ComponentsProcessor cp = new ComponentsProcessor();
        FGraph graph = createTwoComponentsGraph();
        graph.setProperty(ForceOptions.SEPARATE_CONNECTED_COMPONENTS, true);
        List<FGraph> graphs = cp.split(graph);
        assertThat(graphs.size(), is(2));
        graphs.forEach(g -> checkSimpleGraph(g));
    }

    /**
     * Ensure that a single {@link FGraph} is imported if connected components shall not be split.
     */
    @Test
    public void testDoNotSeparateConnectedComponents() {
        ComponentsProcessor cp = new ComponentsProcessor();
        FGraph graph = createTwoComponentsGraph();
        graph.setProperty(ForceOptions.SEPARATE_CONNECTED_COMPONENTS, false);
        List<FGraph> graphs = cp.split(graph);
        assertThat(graphs.size(), is(1));

        FGraph fGraph = graphs.get(0);
        assertThat(fGraph.getNodes().size(), is(6));
        assertThat(fGraph.getEdges().size(), is(4));
        assertThat(fGraph.getLabels().size(), is(4));
    }

    private void checkSimpleGraph(final FGraph fGraph) {
        assertThat(fGraph.getNodes().size(), is(3));
        assertThat(fGraph.getEdges().size(), is(2));
        assertThat(fGraph.getLabels().size(), is(2));
    }

    private ElkNode createElkGraph() {
        ElkNode elkGraph = ElkGraphUtil.createGraph();
        ElkNode n1 = ElkGraphUtil.createNode(elkGraph);
        ElkNode n2 = ElkGraphUtil.createNode(elkGraph);
        ElkNode n3 = ElkGraphUtil.createNode(elkGraph);
        ElkEdge edge = ElkGraphUtil.createSimpleEdge(n1, n2);
        ElkGraphUtil.createLabel("test", edge);
        ElkEdge edge2 = ElkGraphUtil.createSimpleEdge(n1, n3);
        ElkGraphUtil.createLabel("test2", edge2);
        return elkGraph;
    }

    private FGraph createSimpleGraph() {
        ElkGraphImporter importer = new ElkGraphImporter();
        return importer.importGraph(createElkGraph());
    }

    private FGraph createTwoComponentsGraph() {
        ElkNode g1 = createElkGraph();
        ElkNode g2 = createElkGraph();
        ElkNode graph = ElkGraphUtil.createGraph();
        graph.getChildren().addAll(g1.getChildren());
        graph.getChildren().addAll(g2.getChildren());
        ElkGraphImporter importer = new ElkGraphImporter();
        return importer.importGraph(graph);
    }
}
