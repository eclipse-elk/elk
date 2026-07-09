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
package org.eclipse.elk.alg.layered;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.Arrays;
import java.util.Collection;

import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.layered.options.LayeredSpacings;
import org.eclipse.elk.alg.layered.options.LayeredSpacings.LayeredSpacingsBuilder;
import org.eclipse.elk.alg.test.PlainJavaInitialization;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.util.BasicProgressMonitor;
import org.eclipse.elk.core.util.ElkSpacings;
import org.eclipse.elk.core.util.ElkUtil;
import org.eclipse.elk.core.util.IGraphElementVisitor;
import org.eclipse.elk.graph.ElkEdge;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.properties.IProperty;
import org.eclipse.elk.graph.properties.Property;
import org.eclipse.elk.graph.util.ElkGraphUtil;
import org.junit.Test;
import org.junit.experimental.runners.Enclosed;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameter;
import org.junit.runners.Parameterized.Parameters;

/**
 * Tests adjustments by {@link LayeredSpacings} to {@link ElkSpacings}. For instance, whether layered's potentially
 * overwritten default values are considered properly.
 */
@RunWith(Enclosed.class)
public class LayeredSpacingTest {

    static {
        PlainJavaInitialization.initializePlainJavaLayout();
    }

    private static double DOUBLE_EQ_EPSILON = 10e-5;
    
    private static ElkNode createAndConfigureTestGraph(final IGraphElementVisitor configurator) {
        ElkNode simpleGraph = createSimpleGraph();
        ElkUtil.applyVisitors(simpleGraph, configurator);
        return simpleGraph;
    }

    @RunWith(Parameterized.class)
    public static class WithFactorAndValueTest {

        @Parameters(name="{0}")
        public static Collection<Object[]> data(){
            return Arrays.asList(new Object[][] {
                {LayeredOptions.SPACING_EDGE_EDGE},
                {LayeredOptions.SPACING_EDGE_LABEL},
                {LayeredOptions.SPACING_EDGE_NODE},
                {LayeredOptions.SPACING_LABEL_LABEL},
                {LayeredOptions.SPACING_LABEL_NODE},
                {LayeredOptions.SPACING_LABEL_PORT_HORIZONTAL},
                {LayeredOptions.SPACING_LABEL_PORT_VERTICAL},
                {LayeredOptions.SPACING_NODE_SELF_LOOP},
                {LayeredOptions.SPACING_PORT_PORT},
                //
                {LayeredOptions.SPACING_EDGE_EDGE_BETWEEN_LAYERS},
                {LayeredOptions.SPACING_EDGE_NODE_BETWEEN_LAYERS},
                {LayeredOptions.SPACING_NODE_NODE_BETWEEN_LAYERS}
            });
        }
        
        @Parameter(0)
        public IProperty<Double> optionToTest;
        
        /**
         * Check that in the (potentially overwritten) default values of the layered options and not the core options
         * are used.
         */
        @Test
        public void testDefaultFactorsProperlyOverwritten() {
            LayeredSpacingsBuilder builder = LayeredSpacings.withBaseValue(33.0);
            double value = builder.getFactors().get(optionToTest);
            assertEquals(optionToTest.getDefault() / LayeredSpacingsBuilder.BASE_SPACING_OPTION.getDefault(), value,
                    DOUBLE_EQ_EPSILON);
        }

    }
    
    /**
     * Tests that the spacing values have actually properly been regarded in the final drawing. 
     * 
     * TODO only {@link LayeredOptions#SPACING_NODE_NODE_BETWEEN_LAYERS} is tested so far - add tests for the others!
     */
    public static class AfterLayoutTest {
        
        @Test
        public void testSpacingNodeNodeBetweenLayers() {
            ElkNode graph = createAndConfigureTestGraph(
                    LayeredSpacings.withBaseValue(33.0)
                                   .withFactor(LayeredOptions.SPACING_NODE_NODE_BETWEEN_LAYERS, 2.0)
                                   .toVisitor());
            
            LayeredLayoutProvider layoutProvider = new LayeredLayoutProvider();
            layoutProvider.layout(graph, new BasicProgressMonitor());
            
            for (ElkEdge edge : graph.getContainedEdges()) {
                final ElkNode source = ElkGraphUtil.getSourceNode(edge);
                final double spacingNodeNodeBetweenLayers =
                        ElkGraphUtil.getTargetNode(edge).getX() - (source.getX() + source.getWidth());
                assertEquals(66, spacingNodeNodeBetweenLayers, DOUBLE_EQ_EPSILON);
            }
        }

    }
    
    /**
     * As described in #104 there is a conceptual problem with spacings default values. Code that is implemented
     * independent of the layout algorithm (such as node label placement code) is not aware of the default values that,
     * for instance, the layered algorithm specifies. It can only access the layout options via the {@link CoreOptions}.
     * 
     * However, as soon as a spacing configuration is applied to a graph, spacing values should have been set
     * explicitly. Thus, the problem with the default values must not occur in that case.
     */
    public static class SpacingOverrides {
        

        final static IProperty<Double> PROPERTY_WITH_DIFFERENT_DEFAULT = new Property<Double>(
                LayeredOptions.SPACING_EDGE_NODE, LayeredOptions.SPACING_EDGE_NODE.getDefault() + 3.0);

        @Test 
        public void testExtractsCorrectDefaultForLayeredProperty() {
            ElkNode node = createSimpleGraph();
            final double snn = node.getProperty(LayeredOptions.SPACING_EDGE_NODE);
            assertEquals(LayeredOptions.SPACING_EDGE_NODE.getDefault(), snn, DOUBLE_EQ_EPSILON);
            // Note that after the 'getProperty' call the default value is actually _set_ on the node
        }
     
        @Test 
        public void testExtractsCorrectDefaultForInheritingProperty() {
            ElkNode node = createSimpleGraph();
            final double sne = node.getProperty(PROPERTY_WITH_DIFFERENT_DEFAULT);
            assertEquals(LayeredOptions.SPACING_EDGE_NODE.getDefault() + 3.0, sne, DOUBLE_EQ_EPSILON);
            // Note that after the 'getProperty' call the default value is actually _set_ on the node
        }
        
        @Test
        public void testExtractsCorrectDefaultValueAfterSpacingConfiguration() {
            ElkNode node = createSimpleGraph();
            LayeredSpacingsBuilder builder = LayeredSpacings.withBaseValue(35.0);
            builder.apply(node);
            
            
            assertTrue(node.hasProperty(LayeredOptions.SPACING_EDGE_NODE));
            final double expected = 35.0 * builder.getFactors().get(LayeredOptions.SPACING_EDGE_NODE);
            assertEquals(expected, node.getProperty(PROPERTY_WITH_DIFFERENT_DEFAULT), DOUBLE_EQ_EPSILON);
            assertEquals(expected, node.getProperty(LayeredOptions.SPACING_EDGE_NODE), DOUBLE_EQ_EPSILON);
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
