/*******************************************************************************
 * Copyright (c) 2012, 2020 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.layered.components;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.Set;

import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.core.options.PortSide;
import org.junit.Test;

/**
 * Tests the {@link ComponentGroup} class.
 */
public class ComponentGroupTest {
    /**
     * Tests some valid combinations of layered graphs in a single component group.
     */
    @Test
    public void testValidConstraints() {
        // Add two N, E, W, S, and C components each to a single compound group
        ComponentGroup group = new ComponentGroup();
        assertTrue(group.add(generateGraph(PortSide.SIDES_NORTH)));
        assertTrue(group.add(generateGraph(PortSide.SIDES_NORTH)));
        assertTrue(group.add(generateGraph(PortSide.SIDES_SOUTH)));
        assertTrue(group.add(generateGraph(PortSide.SIDES_SOUTH)));
        assertTrue(group.add(generateGraph(PortSide.SIDES_WEST)));
        assertTrue(group.add(generateGraph(PortSide.SIDES_WEST)));
        assertTrue(group.add(generateGraph(PortSide.SIDES_EAST)));
        assertTrue(group.add(generateGraph(PortSide.SIDES_EAST)));
        assertTrue(group.add(generateGraph(PortSide.SIDES_NONE)));
        assertTrue(group.add(generateGraph(PortSide.SIDES_NONE)));
        
        // Add a horizontal and two E, W, and C components each to a single compound group
        group = new ComponentGroup();
        assertTrue(group.add(generateGraph(PortSide.SIDES_EAST_WEST)));
        assertTrue(group.add(generateGraph(PortSide.SIDES_EAST_WEST)));
        assertTrue(group.add(generateGraph(PortSide.SIDES_WEST)));
        assertTrue(group.add(generateGraph(PortSide.SIDES_WEST)));
        assertTrue(group.add(generateGraph(PortSide.SIDES_EAST)));
        assertTrue(group.add(generateGraph(PortSide.SIDES_EAST)));
        assertTrue(group.add(generateGraph(PortSide.SIDES_NONE)));
        assertTrue(group.add(generateGraph(PortSide.SIDES_NONE)));
        
        // Add a vertical and two N, S, and C components each to a single compound group
        group = new ComponentGroup();
        assertTrue(group.add(generateGraph(PortSide.SIDES_NORTH_SOUTH)));
        assertTrue(group.add(generateGraph(PortSide.SIDES_NORTH_SOUTH)));
        assertTrue(group.add(generateGraph(PortSide.SIDES_NORTH)));
        assertTrue(group.add(generateGraph(PortSide.SIDES_NORTH)));
        assertTrue(group.add(generateGraph(PortSide.SIDES_SOUTH)));
        assertTrue(group.add(generateGraph(PortSide.SIDES_SOUTH)));
        assertTrue(group.add(generateGraph(PortSide.SIDES_NONE)));
        assertTrue(group.add(generateGraph(PortSide.SIDES_NONE)));
        
        // A connected component connected to all four port sides allows corner components
        group = new ComponentGroup();
        assertTrue(group.add(generateGraph(PortSide.SIDES_NORTH_EAST_SOUTH_WEST)));
        assertTrue(group.add(generateGraph(PortSide.SIDES_NORTH_WEST)));
        assertTrue(group.add(generateGraph(PortSide.SIDES_NORTH_EAST)));
        assertTrue(group.add(generateGraph(PortSide.SIDES_SOUTH_WEST)));
        assertTrue(group.add(generateGraph(PortSide.SIDES_EAST_SOUTH)));
    }
    
    /**
     * Tests the invalid combinations of layered graphs in a single component group.
     */
    @Test
    public void testInvalidConstraints() {
        // Add a horizontal and a vertical component to a single component group
        ComponentGroup group = new ComponentGroup();
        assertTrue(group.add(generateGraph(PortSide.SIDES_EAST_WEST)));
        assertFalse(group.add(generateGraph(PortSide.SIDES_NORTH_SOUTH)));
        assertFalse(group.add(generateGraph(PortSide.SIDES_NORTH_SOUTH_WEST)));
        assertFalse(group.add(generateGraph(PortSide.SIDES_NORTH_EAST_SOUTH)));
        
        // Do the reverse
        group = new ComponentGroup();
        assertTrue(group.add(generateGraph(PortSide.SIDES_NORTH_SOUTH)));
        assertFalse(group.add(generateGraph(PortSide.SIDES_EAST_WEST)));
        assertFalse(group.add(generateGraph(PortSide.SIDES_NORTH_EAST_WEST)));
        assertFalse(group.add(generateGraph(PortSide.SIDES_EAST_SOUTH_WEST)));
        
        // Add graphs to each corner
        group = new ComponentGroup();
        assertTrue(group.add(generateGraph(PortSide.SIDES_NORTH_WEST)));
        assertTrue(group.add(generateGraph(PortSide.SIDES_NORTH_EAST)));
        assertTrue(group.add(generateGraph(PortSide.SIDES_SOUTH_WEST)));
        assertTrue(group.add(generateGraph(PortSide.SIDES_EAST_SOUTH)));
        assertFalse(group.add(generateGraph(PortSide.SIDES_NORTH_EAST_SOUTH)));
        assertFalse(group.add(generateGraph(PortSide.SIDES_NORTH_EAST_WEST)));
        assertFalse(group.add(generateGraph(PortSide.SIDES_EAST_SOUTH_WEST)));
        assertFalse(group.add(generateGraph(PortSide.SIDES_NORTH_SOUTH_WEST)));
        
        // A connected component connected to all four port sides collides with everything
        // except for corners
        group = new ComponentGroup();
        assertTrue(group.add(generateGraph(PortSide.SIDES_NORTH_EAST_SOUTH_WEST)));
        assertFalse(group.add(generateGraph(PortSide.SIDES_NONE)));
        assertFalse(group.add(generateGraph(PortSide.SIDES_NORTH_SOUTH)));
        assertFalse(group.add(generateGraph(PortSide.SIDES_EAST_WEST)));
        assertFalse(group.add(generateGraph(PortSide.SIDES_NORTH_EAST_SOUTH)));
        assertFalse(group.add(generateGraph(PortSide.SIDES_NORTH_EAST_WEST)));
        assertFalse(group.add(generateGraph(PortSide.SIDES_EAST_SOUTH_WEST)));
        assertFalse(group.add(generateGraph(PortSide.SIDES_NORTH_SOUTH_WEST)));
        assertFalse(group.add(generateGraph(PortSide.SIDES_NORTH_EAST_SOUTH_WEST)));
    }
    
    
    /**
     * Creates a new layered graph with its external port side connections property set to the
     * given set of connections.
     * 
     * @param connections the set of connections.
     * @return the layered graph.
     */
    private static LGraph generateGraph(final Set<PortSide> connections) {
        LGraph graph = new LGraph();
        graph.setProperty(InternalProperties.EXT_PORT_CONNECTIONS, connections);
        
        return graph;
    }
    
}
