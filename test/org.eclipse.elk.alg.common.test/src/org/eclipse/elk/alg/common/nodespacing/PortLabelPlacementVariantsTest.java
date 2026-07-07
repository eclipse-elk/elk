/*******************************************************************************
 * Copyright (c) 2020 Kiel University and others.
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
package org.eclipse.elk.alg.common.nodespacing;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertThat;
import static org.junit.Assert.assertTrue;

import java.util.Collections;
import java.util.List;

import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.test.framework.LayoutTestRunner;
import org.eclipse.elk.alg.test.framework.annotations.Algorithm;
import org.eclipse.elk.alg.test.framework.annotations.ConfiguratorProvider;
import org.eclipse.elk.alg.test.framework.annotations.DefaultConfiguration;
import org.eclipse.elk.alg.test.framework.annotations.GraphResourceProvider;
import org.eclipse.elk.alg.test.framework.io.AbstractResourcePath;
import org.eclipse.elk.alg.test.framework.io.ModelResourcePath;
import org.eclipse.elk.core.LayoutConfigurator;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.options.PortLabelPlacement;
import org.eclipse.elk.core.options.PortSide;
import org.eclipse.elk.graph.ElkLabel;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.ElkPort;
import org.eclipse.emf.common.util.ECollections;
import org.hamcrest.number.OrderingComparison;
import org.junit.Test;
import org.junit.runner.RunWith;

import com.google.common.collect.Lists;

/**
 * Tests the various variants of port label placement that are configurable using {@link PortLabelPlacement}.
 */
@RunWith(LayoutTestRunner.class)
@Algorithm(LayeredOptions.ALGORITHM_ID)
@DefaultConfiguration(nodes = true, ports = true)
public class PortLabelPlacementVariantsTest {

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Sources

    @GraphResourceProvider
    public List<AbstractResourcePath> testGraphs() {
        return Lists.newArrayList(
                new ModelResourcePath("tests/core/label_placement/port_labels/variants.elkt"));
    }
    
    @ConfiguratorProvider
    public LayoutConfigurator westSidePorts() {
        return configuratorFor(PortSide.WEST);
    }

    // TODO activate if behavior for NORTH/SOUTH is consistent
    // @ConfiguratorProvider
    // public LayoutConfigurator northSidePort() {
    // return configuratorFor(PortSide.NORTH);
    // }
    
    @ConfiguratorProvider
    public LayoutConfigurator eastSidePorts() {
        return configuratorFor(PortSide.EAST);
    }
    
    // TODO activate if behavior for NORTH/SOUTH is consistent
    // @ConfiguratorProvider
    // public LayoutConfigurator southSidePorts() {
    // return configuratorFor(PortSide.SOUTH);
    // }

    private LayoutConfigurator configuratorFor(final PortSide ps) {
        LayoutConfigurator config = new LayoutConfigurator();
        // By default this would override already set options
        config.configure(ElkPort.class).setProperty(LayeredOptions.PORT_SIDE, ps);
        return config;
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Tests
    
    /* - - -  OUTSIDE - - -  */
    
    @Test
    public void testOutsideTwoDefault(final ElkNode graph) {
        List<ElkPort> ports = getPortsOfNode(graph, "outside_two_default");
        assertEquals("Unexpected test graph", 2, ports.size());
        assertBelowOrRight(ports.get(0));
        assertAboveOrLeft(ports.get(1));
    }
    
    @Test
    public void testOutsideTwoSameSide(final ElkNode graph) {
        List<ElkPort> ports = getPortsOfNode(graph, "outside_two_same_side");
        assertEquals("Unexpected test graph", 2, ports.size());
        assertBelowOrRight(ports.get(0));
        assertBelowOrRight(ports.get(1));
    }
    
    @Test
    public void testOutsideTwoNextToPort(final ElkNode graph) {
        List<ElkPort> ports = getPortsOfNode(graph, "outside_two_next_to_port");
        assertEquals("Unexpected test graph", 2, ports.size());
        assertCentered(ports.get(0));
        assertCentered(ports.get(1));
    }
    
    @Test
    public void testOutsideThreeDefault(final ElkNode graph) {
        List<ElkPort> ports = getPortsOfNode(graph, "outside_three_default");
        assertEquals("Unexpected test graph", 3, ports.size());
        assertBelowOrRight(ports.get(0));
        assertBelowOrRight(ports.get(1));
        assertBelowOrRight(ports.get(2));
    }
    
    @Test
    public void testOutsideThreeSpaceEfficient(final ElkNode graph) {
        List<ElkPort> ports = getPortsOfNode(graph, "outside_three_space_efficient");
        assertEquals("Unexpected test graph", 3, ports.size());
        assertBelowOrRight(ports.get(0));
        assertBelowOrRight(ports.get(1));
        assertAboveOrLeft(ports.get(2));
    }
    
    /* - - -  INSIDE - - -  */
    
    @Test
    public void testInsideTwoDefault(final ElkNode graph) {
        List<ElkPort> ports = getPortsOfNode(graph, "inside_two_default");
        assertEquals("Unexpected test graph", 2, ports.size());
        assertCentered(ports.get(0));
        assertCentered(ports.get(1));
    }

    @Test
    public void testInsideTwoDefaultHierarchical(final ElkNode graph) {
        List<ElkPort> ports = getPortsOfNode(graph, "inside_two_default_hierarchical");
        assertEquals("Unexpected test graph", 2, ports.size());
        assertBelowOrRight(ports.get(0));
        assertBelowOrRight(ports.get(1));
    }

    @Test
    public void testInsideTwoNextToPortHierarchical(final ElkNode graph) {
        List<ElkPort> ports = getPortsOfNode(graph, "inside_two_next_to_port_hierarchical");
        assertEquals("Unexpected test graph", 2, ports.size());
        assertCentered(ports.get(0));
        assertCentered(ports.get(1));
    }

    @Test
    public void testInsideTwoWithOneEdge(final ElkNode graph) {
        List<ElkPort> ports = getPortsOfNode(graph, "inside_two_with_one_edge");
        assertEquals("Unexpected test graph", 2, ports.size());
        assertBelowOrRight(ports.get(0));
        assertBelowOrRight(ports.get(1));
    }
    
    @Test
    public void testInsideTwoWithOneEdgeNextToPort(final ElkNode graph) {
        List<ElkPort> ports = getPortsOfNode(graph, "inside_two_with_one_edge_next_to_port");
        assertEquals("Unexpected test graph", 2, ports.size());
        // Unfortunately this test is not symmetric
        ElkPort portWithEdge = ports.get(0).getOutgoingEdges().isEmpty() ? ports.get(1) : ports.get(0);
        ElkPort theOther = ports.get(0).getOutgoingEdges().isEmpty() ? ports.get(0) : ports.get(1);
        assertBelowOrRight(portWithEdge);
        assertCentered(theOther);
    }
    
    @Test
    public void testInsideThreeWithOneEdge(final ElkNode graph) {
        List<ElkPort> ports = getPortsOfNode(graph, "inside_three_with_one_edge");
        assertEquals("Unexpected test graph", 3, ports.size());
        assertBelowOrRight(ports.get(0));
        assertBelowOrRight(ports.get(1));
        assertBelowOrRight(ports.get(2));
    }
    
    @Test
    public void testInsideThreeWithOneEdgeNextToPort(final ElkNode graph) {
        List<ElkPort> ports = getPortsOfNode(graph, "inside_three_with_one_edge_next_to_port");
        assertEquals("Unexpected test graph", 3, ports.size());
        assertCentered(ports.get(0));
        assertBelowOrRight(ports.get(1));
        assertCentered(ports.get(2));
    }
    
    /* - - - Checks and utility - - - */
    
    private void assertBelowOrRight(final ElkPort port) {
        final ElkLabel label = getLabel(port);
        // Note that port label positions are relative to the port
        double portPosition = 0;
        if (PortSide.SIDES_EAST_WEST.contains(port.getProperty(CoreOptions.PORT_SIDE))) {
            portPosition = label.getY();
        } else if (PortSide.SIDES_NORTH_SOUTH.contains(port.getProperty(CoreOptions.PORT_SIDE))) {
            portPosition = label.getX();
        } else {
            assertTrue(false);
        }
        assertThat(portPosition, OrderingComparison.greaterThan(0.0));
    }

    private void assertAboveOrLeft(final ElkPort port) {
        final ElkLabel label = getLabel(port);
        // Note that port label positions are relative to the port
        double portPosition = 0;
        if (PortSide.SIDES_EAST_WEST.contains(port.getProperty(CoreOptions.PORT_SIDE))) {
            portPosition = label.getY();
        } else if (PortSide.SIDES_NORTH_SOUTH.contains(port.getProperty(CoreOptions.PORT_SIDE))) {
            portPosition = label.getX();
        } else {
            assertTrue(false);
        }
        assertThat(portPosition, OrderingComparison.lessThan(0.0));
    }
    
    private void assertCentered(final ElkPort port) {
        final ElkLabel label = getLabel(port);
        if (PortSide.SIDES_EAST_WEST.contains(port.getProperty(CoreOptions.PORT_SIDE))) {
            final double portCenter = port.getHeight() / 2.0;
            final double labelCenter = label.getY() + label.getHeight() / 2.0;
            assertThat(labelCenter, OrderingComparison.comparesEqualTo(portCenter));
        } else if (PortSide.SIDES_NORTH_SOUTH.contains(port.getProperty(CoreOptions.PORT_SIDE))) {
            final double portCenter = port.getWidth() / 2.0;
            final double labelCenter = label.getX() + label.getWidth() / 2.0;
            assertThat(labelCenter, OrderingComparison.comparesEqualTo(portCenter));
        } else {
            assertTrue(false);
        }
    }

    private ElkLabel getLabel(final ElkPort port) {
        return port.getLabels().get(0);
    }

    private List<ElkPort> getPortsOfNode(final ElkNode parent, final String nodeId) {
        return parent.getChildren().stream()
            .filter(n -> n.getIdentifier().equals(nodeId))
            .map(n -> n.getPorts())
            .findFirst()
            // Sort the ports based on their index (which by definition is specified clock-wise).
            // For west and east sort them top-down, and for north and south sort them left-to-right.
            .map(l -> {
                ECollections.sort(l, (p1, p2) -> {
                    if (PortSide.SIDES_SOUTH_WEST.contains(p1.getProperty(CoreOptions.PORT_SIDE))) {
                        return p1.getProperty(LayeredOptions.PORT_INDEX)
                                - p2.getProperty(LayeredOptions.PORT_INDEX);
                    } else {
                        return p2.getProperty(LayeredOptions.PORT_INDEX)
                                - p1.getProperty(LayeredOptions.PORT_INDEX);
                    }
                });
                return l;
            })
            .map(l -> (List<ElkPort>) l) // only for the next line to work
            .orElse(Collections.emptyList());
    }
}
