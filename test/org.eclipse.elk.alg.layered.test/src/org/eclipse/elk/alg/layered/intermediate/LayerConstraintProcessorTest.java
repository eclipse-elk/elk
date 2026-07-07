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

import java.util.Collection;
import java.util.List;

import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.options.LayerConstraint;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.test.framework.LayoutTestRunner;
import org.eclipse.elk.alg.test.framework.annotations.Algorithm;
import org.eclipse.elk.alg.test.framework.annotations.DefaultConfiguration;
import org.eclipse.elk.alg.test.framework.annotations.GraphResourceProvider;
import org.eclipse.elk.alg.test.framework.annotations.TestAfterProcessor;
import org.eclipse.elk.alg.test.framework.io.AbstractResourcePath;
import org.eclipse.elk.alg.test.framework.io.ModelResourcePath;
import org.eclipse.elk.core.options.LabelSide;
import org.junit.runner.RunWith;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;

@RunWith(LayoutTestRunner.class)
@Algorithm(LayeredOptions.ALGORITHM_ID)
@DefaultConfiguration()
public class LayerConstraintProcessorTest {

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Sources

    @GraphResourceProvider
    public List<AbstractResourcePath> testGraphs() {
        return Lists.newArrayList(new ModelResourcePath("tests/layered/layer_constraints/**/"));
    }


    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Tests

    /**
     * If the first layer contains {@link LayerConstraint#FIRST_SEPARATE} no other type is allowed there. All
     * {@link LayerConstraint#FIRST} constraints must be in the next layer. Similar for
     * {@link LayerConstraint#LAST_SEPARATE} and {@link LayerConstraint#LAST}.
     */
    @TestAfterProcessor(LayerConstraintPostprocessor.class)
    public void testLayerConstraints(final Object graph) {
        LGraph lGraph = (LGraph) graph;
        
        Multimap<LayerConstraint, Integer> constraintMap = HashMultimap.create();

        int nextLayerId = 0;
        for (Layer layer : lGraph) {
            layer.id = nextLayerId++;
            for (LNode node : layer) {
                constraintMap.put(node.getProperty(LayeredOptions.LAYERING_LAYER_CONSTRAINT), layer.id);
            }
        }

        // if FIRST_SEPARATE nodes are available
        if (constraintMap.containsKey(LayerConstraint.FIRST_SEPARATE)) {
            // FIRST_SEPARATE must only occur in the first layer
            Collection<Integer> firstSep = constraintMap.get(LayerConstraint.FIRST_SEPARATE);
            assertTrue(firstSep.size() == 1 && firstSep.contains(0));

            // first layer must not contain a NONE node
            assertTrue(!constraintMap.get(LayerConstraint.NONE).contains(0));

            // all FIRST nodes have to be at position 1
            Collection<Integer> first = constraintMap.get(LayerConstraint.FIRST);
            if (!first.isEmpty()) {
                assertTrue(first.size() == 1 && first.contains(1));
            }
        } else {
            // all FIRST nodes must be at position 0
            Collection<Integer> first = constraintMap.get(LayerConstraint.FIRST);
            if (!first.isEmpty()) {
                assertTrue(first.size() == 1 && first.contains(0));
            }
        }

        int lastLayer = lGraph.getLayers().size() - 1;
        
        // if LAST_SEPARATE nodes are available
        if (constraintMap.containsKey(LayerConstraint.LAST_SEPARATE)) {
            // LAST_SEPARATE must only occur in the last layer
            Collection<Integer> lastSep = constraintMap.get(LayerConstraint.LAST_SEPARATE);
            assertTrue(lastSep.size() == 1 && lastSep.contains(lastLayer));

            // last layer must not contain a NONE node
            assertTrue(!constraintMap.get(LayerConstraint.NONE).contains(lastLayer));

            // all LAST nodes have to be at position lastLayer - 1
            Collection<Integer> last = constraintMap.get(LayerConstraint.LAST);
            if (!last.isEmpty()) {
                assertTrue(last.size() == 1 && last.contains(lastLayer - 1));
            }
        } else {
            // all LAST nodes must be at position lastLayer
            Collection<Integer> last = constraintMap.get(LayerConstraint.LAST);
            if (!last.isEmpty()) {
                assertTrue(last.size() == 1 && last.contains(lastLayer));
            }
        }
    }
    
    /**
     * All labels on ports and edges have an assigned {@link LabelSide}.
     */
    @TestAfterProcessor(LayerConstraintPostprocessor.class)
    public void noEmptyLayers(Object graph) {
        LGraph lGraph = (LGraph) graph;
        
        lGraph.getLayers().stream()
            .forEach(layer -> assertTrue(!layer.getNodes().isEmpty()));
    }
    
}
