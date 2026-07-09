/*******************************************************************************
 * Copyright (c) 2019 Kiel University and others.
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
package org.eclipse.elk.alg.layered.intermediate.loops;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.util.HashSet;
import java.util.Set;
import java.util.stream.Collectors;

import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.junit.Test;

/**
 * Tests whether the data model created by {@link SelfLoopHolder} is correct.
 */
public class SelfLoopHolderTest {

    @Test
    public void testNonSelfLoopGraph() {
        LGraph lGraph = SelfLoopTestGraphCreator.basicGraphWithoutSelfLoops();

        for (LNode lNode : lGraph.getLayerlessNodes()) {
            assertFalse(SelfLoopHolder.needsSelfLoopProcessing(lNode));
            
            // This should throw an AssertionError
            try {
                SelfLoopHolder.install(lNode);
                fail();
            } catch (AssertionError e) {
                // If this happens, everything's good
            }
        }
    }

    @Test
    public void testBasicSelfLoops() {
        LGraph lGraph = SelfLoopTestGraphCreator.basicGraphWithSelfLoops();

        // The second node is the one with self loops
        LNode lNode = lGraph.getLayerlessNodes().get(1);
        assertTrue(SelfLoopHolder.needsSelfLoopProcessing(lNode));

        SelfLoopHolder slHolder = SelfLoopHolder.install(lNode);
        assertNotNull(slHolder);
        assertEquals(slHolder, lNode.getProperty(InternalProperties.SELF_LOOP_HOLDER));
        
        // All ports have incident self loops
        assertEquals(lNode.getPorts().size(), slHolder.getSLPortMap().size());
        
        // There should be half as many hyper loops
        assertEquals(lNode.getPorts().size() / 2, slHolder.getSLHyperLoops().size());
    }
    
    @Test
    public void testDifferentSelfHyperLoops() {
        LGraph lGraph = SelfLoopTestGraphCreator.basicGraphWithoutSelfLoops();
        
        LNode lNode = lGraph.getLayerlessNodes().get(1);
        
        // Add further ports
        SelfLoopTestGraphCreator.ports(lNode, 2, 2, 2, 1);
        
        // Connect the first half of the ports with self loops
        Set<LEdge> firstLoops = new HashSet<>();
        for (int i = 1; i < lNode.getPorts().size() / 2; i++) {
            LEdge lEdge = SelfLoopTestGraphCreator.edge(lNode.getPorts().get(i - 1), lNode.getPorts().get(i));
            firstLoops.add(lEdge);
        }
        
        // Connect the other half
        Set<LEdge> secondLoops = new HashSet<>();
        for (int i = lNode.getPorts().size() / 2 + 1; i < lNode.getPorts().size(); i++) {
            LEdge lEdge = SelfLoopTestGraphCreator.edge(lNode.getPorts().get(i - 1), lNode.getPorts().get(i));
            secondLoops.add(lEdge);
        }
        
        // Install self loop holder and check that we have two self hyper loops
        SelfLoopHolder slHolder = SelfLoopHolder.install(lNode);
        assertEquals(2, slHolder.getSLHyperLoops().size());
        
        // Test the edges of the two hyper loops
        Set<LEdge> firstHLEdges = slHolder.getSLHyperLoops().get(0).getSLEdges().stream()
                .map(slEdge -> slEdge.getLEdge())
                .collect(Collectors.toSet());
        Set<LEdge> secondHLEdges = slHolder.getSLHyperLoops().get(1).getSLEdges().stream()
                .map(slEdge -> slEdge.getLEdge())
                .collect(Collectors.toSet());
        
        // They can match in either way
        boolean match = firstHLEdges.equals(firstLoops) && secondHLEdges.equals(secondLoops);
        match |= firstHLEdges.equals(secondLoops) && secondHLEdges.equals(firstLoops);
        
        assertTrue(match);
    }

}
