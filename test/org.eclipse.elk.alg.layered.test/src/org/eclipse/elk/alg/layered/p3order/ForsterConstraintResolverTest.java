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
package org.eclipse.elk.alg.layered.p3order;

import static org.junit.Assert.assertFalse;

import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.test.framework.LayoutTestRunner;
import org.eclipse.elk.alg.test.framework.annotations.Algorithm;
import org.eclipse.elk.alg.test.framework.annotations.DefaultConfiguration;
import org.eclipse.elk.alg.test.framework.annotations.GraphResourceProvider;
import org.eclipse.elk.alg.test.framework.annotations.TestAfterProcessor;
import org.eclipse.elk.alg.test.framework.io.AbstractResourcePath;
import org.eclipse.elk.alg.test.framework.io.ModelResourcePath;
import org.junit.runner.RunWith;

import com.google.common.collect.Lists;

/**
 * Tests whether in-layer constraints are respected. This test runs after the crossing minimizer which makes use of the
 * {@link ForsterConstraintResolver}.
 */
@RunWith(LayoutTestRunner.class)
@Algorithm(LayeredOptions.ALGORITHM_ID)
@DefaultConfiguration()
public class ForsterConstraintResolverTest {

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Sources

    @GraphResourceProvider
    public List<AbstractResourcePath> testGraphs() {
        return Lists.newArrayList(
                new ModelResourcePath("tickets/layered/528_nodeEdgeOverlap.elkt"),
                new ModelResourcePath("tickets/layered/528_violatedNodeOrder.elkt"));
    }


    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Tests

    /**
     * Test whether the order of nodes respects each node's {@link InternalProperties#IN_LAYER_SUCCESSOR_CONSTRAINTS}.
     */
    @TestAfterProcessor(LayerSweepCrossingMinimizer.class)
    public void testSuccessorConstraints(final Object graph) {
        for (Layer layer : (LGraph) graph) {
            // We iterate over the layer's nodes and remember the one's we've already seen. If a successor constraint
            // includes a node we've already seen, we fail
            Set<LNode> encounteredNodes = new HashSet<>();
            for (LNode node : layer) {
                if (node.hasProperty(InternalProperties.IN_LAYER_SUCCESSOR_CONSTRAINTS)) {
                    boolean violation = node.getProperty(InternalProperties.IN_LAYER_SUCCESSOR_CONSTRAINTS).stream()
                            .anyMatch(n -> encounteredNodes.contains(n));
                    assertFalse("Successor constraints are not respected for " + node, violation);
                }
                
                encounteredNodes.add(node);
            }
        }
    }
    
    /**
     * Nodes may be part of layout units, which must be adhered to even after in-layer constraints are resolved.
     */
    @TestAfterProcessor(LayerSweepCrossingMinimizer.class)
    public void testNonOverlappingLayoutUnits(final Object graph) {
        for (Layer layer : (LGraph) graph) {
            // We iterate over the layer's nodes and keep track of two things: the layout units we've encountered so
            // far, and the current layout unit. Layout units are identified through a representing node.
            Set<LNode> encounteredLayoutUnits = new HashSet<>();
            LNode currentLayoutUnit = null;
            
            for (LNode node : layer) {
                LNode newLayoutUnit = node.getProperty(InternalProperties.IN_LAYER_LAYOUT_UNIT);
                
                if (newLayoutUnit != null && newLayoutUnit != currentLayoutUnit) {
                    // We have a node which does not belong to the current layout unit. That's okay if the new layout
                    // unit has not been encountered before (if it has, nodes of different layout units overlap)
                    assertFalse("Node " + newLayoutUnit + " is overlapping with the current layout unit",
                            encounteredLayoutUnits.contains(newLayoutUnit));
                    
                    if (currentLayoutUnit != null) {
                        encounteredLayoutUnits.add(currentLayoutUnit);
                    }
                    
                    currentLayoutUnit = newLayoutUnit;
                }
            }
        }
    }

}
