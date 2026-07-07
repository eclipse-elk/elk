/*******************************************************************************
 * Copyright (c) 2020 Kiel University and others.
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 * 
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.common.nodespacing;

import static org.hamcrest.Matchers.closeTo;
import static org.hamcrest.Matchers.greaterThan;
import static org.hamcrest.Matchers.lessThan;
import static org.junit.Assert.assertThat;

import java.util.List;

import org.eclipse.elk.alg.common.NodeMicroLayout;
import org.eclipse.elk.alg.force.options.ForceOptions;
import org.eclipse.elk.alg.force.options.StressOptions;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.mrtree.options.MrTreeOptions;
import org.eclipse.elk.alg.radial.options.RadialOptions;
import org.eclipse.elk.alg.rectpacking.options.RectPackingOptions;
import org.eclipse.elk.alg.test.framework.LayoutTestRunner;
import org.eclipse.elk.alg.test.framework.annotations.Algorithm;
import org.eclipse.elk.alg.test.framework.annotations.Configurator;
import org.eclipse.elk.alg.test.framework.annotations.DefaultConfiguration;
import org.eclipse.elk.alg.test.framework.annotations.GraphResourceProvider;
import org.eclipse.elk.alg.test.framework.io.AbstractResourcePath;
import org.eclipse.elk.alg.test.framework.io.ModelResourcePath;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.graph.ElkLabel;
import org.eclipse.elk.graph.ElkNode;
import org.junit.Test;
import org.junit.runner.RunWith;

import com.google.common.collect.Lists;

/**
 * Tests the successful invocation of {@link NodeMicroLayout} for algorithms that support it.
 */
@RunWith(LayoutTestRunner.class)
@Algorithm(ForceOptions.ALGORITHM_ID)
@Algorithm(LayeredOptions.ALGORITHM_ID)
@Algorithm(MrTreeOptions.ALGORITHM_ID)
@Algorithm(RadialOptions.ALGORITHM_ID)
@Algorithm(RectPackingOptions.ALGORITHM_ID)
@Algorithm(StressOptions.ALGORITHM_ID)
@DefaultConfiguration(nodes = true, ports = true)
public class NodeMicroLayoutTest {

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Sources

    @GraphResourceProvider
    public List<AbstractResourcePath> testGraphs() {
        return Lists.newArrayList(new ModelResourcePath("tests/core/label_placement/node_labels/variants.elkt"));
    }

    @Configurator
    public void omitNodeMicroLayout(final ElkNode graph) {
        if (graph.getProperty(CoreOptions.ALGORITHM) == LayeredOptions.ALGORITHM_ID) {
            // Layered does not provide an option to omit the node micro layout.
            return;
        }
        graph.setProperty(CoreOptions.OMIT_NODE_MICRO_LAYOUT, true);
    }

    @Configurator
    public void applyNodeMicroLayout(final ElkNode graph) {
        graph.setProperty(CoreOptions.OMIT_NODE_MICRO_LAYOUT, null);
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Tests

    @Test
    public void testMinimumSize(final ElkNode graph) {
        final ElkNode node = graph.getChildren().get(0);
        if (!omitMicroLayout(graph)) {
            assertThat(node.getWidth(), closeTo(101, 0.01));
            assertThat(node.getHeight(), closeTo(102, 0.01));
        } else {
            assertThat(node.getWidth(), closeTo(100, 0.01));
            assertThat(node.getHeight(), closeTo(100, 0.01));
        }
    }

    @Test
    public void testNodeLabelPositions(final ElkNode graph) {
        final ElkNode node = graph.getChildren().get(0);
        if (!omitMicroLayout(graph)) {
            verifyNodeLabelsHaveBeenPositioned(node);
        } else {
            verifyNodeLabelsHaveNotBeenPositioned(node);
        }
    }

    private void verifyNodeLabelsHaveBeenPositioned(final ElkNode node) {
        final double beforeCenterThreshold = 40;
        final double afterCenterThreshold = 60;

        // INSIDE V_TOP H_CENTER
        ElkLabel a = getLabel(node, "A");
        assertThat(a.getY(), greaterThan(0.0));

        // INSIDE V_BOTTOM H_RIGHT
        ElkLabel b = getLabel(node, "B");
        assertThat(b.getX(), greaterThan(afterCenterThreshold));
        assertThat(b.getX(), lessThan(node.getWidth()));
        assertThat(b.getY(), greaterThan(afterCenterThreshold));
        assertThat(b.getY(), lessThan(node.getHeight()));

        // OUTSIDE V_BOTTOM H_CENTER
        ElkLabel c = getLabel(node, "C");
        assertThat(c.getX(), greaterThan(beforeCenterThreshold));
        assertThat(c.getX(), lessThan(afterCenterThreshold));
        assertThat(c.getY(), greaterThan(node.getHeight()));

        // OUTSIDE V_CENTER H_LEFT
        ElkLabel d = getLabel(node, "D");
        assertThat(d.getX(), lessThan(0.0));
        assertThat(d.getY(), greaterThan(beforeCenterThreshold));
        assertThat(d.getY(), lessThan(afterCenterThreshold));
    }

    private boolean omitMicroLayout(final ElkNode graph) {
        return graph.hasProperty(CoreOptions.OMIT_NODE_MICRO_LAYOUT)
                && graph.getProperty(CoreOptions.OMIT_NODE_MICRO_LAYOUT);
    }

    private void verifyNodeLabelsHaveNotBeenPositioned(final ElkNode node) {
        node.getLabels().forEach(label -> assertThat(label.getX(), closeTo(0.0, 0.01)));
    }

    private ElkLabel getLabel(final ElkNode node, final String text) {
        return node.getLabels().stream().filter(label -> label.getText().equals(text)).findFirst().get();
    }

}
