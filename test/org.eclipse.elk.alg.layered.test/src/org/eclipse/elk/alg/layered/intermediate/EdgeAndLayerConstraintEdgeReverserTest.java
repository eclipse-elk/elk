/*******************************************************************************
 * Copyright (c) 2013, 2020 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.layered.intermediate;

import static org.junit.Assert.assertTrue;

import java.util.List;

import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.options.LayerConstraint;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.test.framework.LayoutTestRunner;
import org.eclipse.elk.alg.test.framework.annotations.Algorithm;
import org.eclipse.elk.alg.test.framework.annotations.GraphResourceProvider;
import org.eclipse.elk.alg.test.framework.annotations.TestAfterProcessor;
import org.eclipse.elk.alg.test.framework.io.AbstractResourcePath;
import org.eclipse.elk.alg.test.framework.io.ModelResourcePath;
import org.junit.runner.RunWith;

import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;

/**
 * Basic tests for the {@link EdgeAndLayerConstraintEdgeReverser}.
 */
@RunWith(LayoutTestRunner.class)
@Algorithm(LayeredOptions.ALGORITHM_ID)
public class EdgeAndLayerConstraintEdgeReverserTest {

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Sources

    @GraphResourceProvider
    public List<AbstractResourcePath> testGraphs() {
        return Lists.newArrayList(new ModelResourcePath("tests/layered/**/"));
    }


    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Tests
    
    @TestAfterProcessor(EdgeAndLayerConstraintEdgeReverser.class)
    public void testLayerConstraints(Object graph) {
        for (LNode node : ((LGraph) graph).getLayerlessNodes()) {
            LayerConstraint constr = node.getProperty(LayeredOptions.LAYERING_LAYER_CONSTRAINT);
            if (constr == LayerConstraint.FIRST) {
                // either no incoming edges or only from FIRST_SEPARATE nodes
                for (LEdge inc : node.getIncomingEdges()) {
                    LNode src = inc.getSource().getNode();
                    assertTrue(src.getProperty(LayeredOptions.LAYERING_LAYER_CONSTRAINT)
                            == LayerConstraint.FIRST_SEPARATE);
                }
                
            } else if (constr == LayerConstraint.LAST) {
                // either no outgoing edges or only to LAST_SEPARATE nodes
                for (LEdge inc : node.getOutgoingEdges()) {
                    LNode src = inc.getTarget().getNode();
                    assertTrue(src.getProperty(LayeredOptions.LAYERING_LAYER_CONSTRAINT)
                            == LayerConstraint.LAST_SEPARATE);
                }
                
            } else if (constr == LayerConstraint.FIRST_SEPARATE) {
                assertTrue(Iterables.isEmpty(node.getIncomingEdges()));
                
            } else if (constr == LayerConstraint.LAST_SEPARATE) {
                assertTrue(Iterables.isEmpty(node.getOutgoingEdges()));
            }
        }
    }

}
