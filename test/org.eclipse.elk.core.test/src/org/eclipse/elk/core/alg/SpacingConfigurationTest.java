/*******************************************************************************
 * Copyright (c) 2020 Kiel University and others.
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
package org.eclipse.elk.core.alg;

import static org.junit.Assert.assertEquals;

import java.util.Arrays;
import java.util.Collection;

import org.eclipse.elk.alg.test.PlainJavaInitialization;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.util.ElkSpacings;
import org.eclipse.elk.core.util.ElkSpacings.ElkCoreSpacingsBuilder;
import org.eclipse.elk.core.util.ElkSpacings.AbstractSpacingsBuilder;
import org.eclipse.elk.core.util.ElkUtil;
import org.eclipse.elk.core.util.IGraphElementVisitor;
import org.eclipse.elk.graph.ElkGraphElement;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.properties.IProperty;
import org.eclipse.elk.graph.util.ElkGraphUtil;
import org.junit.Test;
import org.junit.experimental.runners.Enclosed;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameter;
import org.junit.runners.Parameterized.Parameters;

import com.google.common.collect.Streams;

/**
 * Test the utility class {@link ElkSpacings}, which allows to conveniently configure spacings for the overall graph.
 */
@RunWith(Enclosed.class)
public class SpacingConfigurationTest {

    static {
        PlainJavaInitialization.initializePlainJavaLayout();
    }

    private static double DOUBLE_EQ_EPSILON = 10e-5;
    
    public static ElkNode createAndConfigureTestGraph(final IGraphElementVisitor configurator) {
        ElkNode simpleGraph = createSimpleGraph();
        ElkUtil.applyVisitors(simpleGraph, configurator);
        return simpleGraph;
    }

    public static void checkOptionValue(ElkNode graph, Class<? extends ElkGraphElement> clazz,
            IProperty<Double> option, final double value) {
        // It's important to also check the root node here.
        Streams.stream(ElkGraphUtil.propertiesSkippingIteratorFor(graph, true))
                .filter(e -> clazz.isInstance(e))
                .map(ElkGraphElement.class::cast).forEach(e -> {
                    assertEquals(option.getId(), value, e.getProperty(option), DOUBLE_EQ_EPSILON);
                });
    }

    public static class Basics {
        
        @Test
        public void testCreateVisitor() {
            ElkNode node = ElkGraphUtil.createNode(null);
            IGraphElementVisitor visitor = ElkSpacings.withBaseValue(3.0).toVisitor();
            visitor.visit(node);
            assertEquals(3.0, node.getProperty(ElkCoreSpacingsBuilder.BASE_SPACING_OPTION), DOUBLE_EQ_EPSILON);
        }
        
        @Test
        public void testApplyTo() {
            ElkNode node = ElkGraphUtil.createNode(null);
            ElkSpacings.withBaseValue(3.0).apply(node);
            assertEquals(3.0, node.getProperty(ElkCoreSpacingsBuilder.BASE_SPACING_OPTION), DOUBLE_EQ_EPSILON);
        }

        @Test(expected = IllegalArgumentException.class)
        public void testSpacingWithNegativeFactor() {
            createAndConfigureTestGraph(
                    ElkSpacings.withBaseValue(33.0)
                               .withFactor(CoreOptions.SPACING_EDGE_EDGE, -Double.MIN_VALUE)
                               .toVisitor());
        }
        
        @Test(expected = IllegalArgumentException.class)
        public void testSpacingWithNegativeValue() {
            createAndConfigureTestGraph(
                    ElkSpacings.withBaseValue(33.0)
                               .withValue(CoreOptions.SPACING_EDGE_EDGE, -Double.MIN_VALUE)
                               .toVisitor());
        }
    }
    
    /**
     * Tests the basic functionality of {@link AbstractSpacingsBuilder#withFactor(IProperty, double)} and
     * {@link AbstractSpacingsBuilder#withValue(IProperty, double)}.
     */
    @RunWith(Parameterized.class)
    public static class WithFactorAndValueTest {

        @Parameters(name = "{0}")
        public static Collection<Object[]> data() {
            return Arrays.asList(new Object[][] { 
                { CoreOptions.SPACING_COMPONENT_COMPONENT },
                { CoreOptions.SPACING_EDGE_EDGE },
                { CoreOptions.SPACING_EDGE_LABEL },
                { CoreOptions.SPACING_EDGE_NODE },
                { CoreOptions.SPACING_LABEL_LABEL },
                { CoreOptions.SPACING_LABEL_NODE },
                { CoreOptions.SPACING_LABEL_PORT_HORIZONTAL },
                { CoreOptions.SPACING_LABEL_PORT_VERTICAL },
                { CoreOptions.SPACING_NODE_SELF_LOOP }, 
                { CoreOptions.SPACING_PORT_PORT } 
            });
        }
        
        @Parameter(0)
        public IProperty<Double> optionToTest;
        
        @Test
        public void testDefaultFactors() {
            ElkCoreSpacingsBuilder builder = ElkSpacings.withBaseValue(33.0);
            double value = builder.getFactors().get(optionToTest);
            assertEquals(optionToTest.getDefault() / ElkCoreSpacingsBuilder.BASE_SPACING_OPTION.getDefault(), value,
                    DOUBLE_EQ_EPSILON);
        }
        
        @Test
        public void testSpacingWithFactor() {
            ElkNode graph =
                    createAndConfigureTestGraph(ElkSpacings.withBaseValue(33.0)
                                                           .withFactor(optionToTest, 2.0)
                                                           .toVisitor());
            checkOptionValue(graph, ElkNode.class, optionToTest, 66.0);
        }

        @Test
        public void testSpacingWithValue() {
            ElkNode graph =
                    createAndConfigureTestGraph(ElkSpacings.withBaseValue(33.0)
                                                           .withValue(optionToTest, 24.0)
                                                           .toVisitor());
            checkOptionValue(graph, ElkNode.class, optionToTest, 24.0);
        }
        
        @Test
        public void testOverwriteIfRequested() {
            ElkNode graph = createSimpleGraph();
            final ElkNode firstChild = graph.getChildren().get(0);
            final ElkNode secondChild = graph.getChildren().get(1);
            firstChild.setProperty(optionToTest, 22.0);
            ElkUtil.applyVisitors(graph, ElkSpacings.withBaseValue(1)
                                                    .withValue(optionToTest, 3)
                                                    .withOverwrite(true)
                                                    .toVisitor());
            assertEquals("Should overwrite", 3.0, firstChild.getProperty(optionToTest), DOUBLE_EQ_EPSILON);
            assertEquals(3.0, secondChild.getProperty(optionToTest), DOUBLE_EQ_EPSILON);
        }
        
        @Test
        public void testDontOverwriteIfNotRequested() {
            ElkNode graph = createSimpleGraph();
            final ElkNode firstChild = graph.getChildren().get(0);
            final ElkNode secondChild = graph.getChildren().get(1);
            firstChild.setProperty(optionToTest, 22.0);
            ElkUtil.applyVisitors(graph, ElkSpacings.withBaseValue(1)
                                                    .withValue(optionToTest, 3)
                                                    .withOverwrite(false)
                                                    .toVisitor());
            assertEquals("Should not overwrite", 22.0, firstChild.getProperty(optionToTest), DOUBLE_EQ_EPSILON);
            assertEquals(3.0, secondChild.getProperty(optionToTest), DOUBLE_EQ_EPSILON);
        }
    }
    
    /**
     * Tests that trying to configure factors and values for non-spacing layout options is rejected.
     */
    @RunWith(Parameterized.class)
    public static class InvalidOptionConfigurationTest {

        @Parameters(name = "{0}")
        public static Collection<Object[]> data() {
            return Arrays.asList(
                    new Object[][] { { ElkCoreSpacingsBuilder.BASE_SPACING_OPTION }, { CoreOptions.ASPECT_RATIO } });
        }
        
        @Parameter(0)
        public IProperty<Double> optionToTest;
        
        @Test(expected = IllegalArgumentException.class)
        public void testSpacingWithFactor() {
            createAndConfigureTestGraph(
                    ElkSpacings.withBaseValue(33.0)
                               .withFactor(optionToTest, 2.0)
                               .toVisitor());
        }
        
        @Test(expected = IllegalArgumentException.class)
        public void testSpacingWithValue() {
            createAndConfigureTestGraph(
                    ElkSpacings.withBaseValue(33.0)
                               .withValue(optionToTest, 24.0)
                               .toVisitor());
        }
        
    }
    

    /**
     * Create a simple test graph. The graph has at least two nodes and an edge, the nodes have predefined sizes with
     * fixed size constraint, but neither nodes nor edges have predefined coordinates, i.e. they are all set to (0,0).
     * The graph size is initially 0 as well.
     */
    private static ElkNode createSimpleGraph() {
        ElkNode graph = ElkGraphUtil.createGraph();

        ElkNode node1 = ElkGraphUtil.createNode(graph);
        node1.setWidth(30);
        node1.setHeight(30);

        ElkNode node2 = ElkGraphUtil.createNode(graph);
        node2.setWidth(30);
        node2.setHeight(30);

        ElkGraphUtil.createSimpleEdge(node1, node2);

        return graph;
    }

}
