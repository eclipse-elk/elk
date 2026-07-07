/*******************************************************************************
 * Copyright (c) 2013, 2020 Kiel University and others.
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
package org.eclipse.elk.alg.layered.intermediate;

import static org.junit.Assert.fail;

import java.util.List;

import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.test.framework.LayoutTestRunner;
import org.eclipse.elk.alg.test.framework.annotations.Algorithm;
import org.eclipse.elk.alg.test.framework.annotations.GraphResourceProvider;
import org.eclipse.elk.alg.test.framework.annotations.TestAfterProcessor;
import org.eclipse.elk.alg.test.framework.io.AbstractResourcePath;
import org.eclipse.elk.alg.test.framework.io.ModelResourcePath;
import org.eclipse.elk.core.options.PortSide;
import org.junit.runner.RunWith;

import com.google.common.collect.Lists;

/**
 * Basic tests for the {@link NorthSouthPortPostprocessor}.
 */
@RunWith(LayoutTestRunner.class)
@Algorithm(LayeredOptions.ALGORITHM_ID)
public class NorthSouthPortPreprocessorTest {

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Sources

    @GraphResourceProvider
    public List<AbstractResourcePath> testGraphs() {
        return Lists.newArrayList(new ModelResourcePath("tests/layered/north_south_ports/**"));
    }


    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Tests
    
    @TestAfterProcessor(NorthSouthPortPreprocessor.class)
    public void testIsolatedNorthSouthPorts(Object graph) {
        ((LGraph) graph).getLayers().stream()
            .flatMap(layer -> layer.getNodes().stream())
            .filter(node -> node.getProperty(LayeredOptions.PORT_CONSTRAINTS).isSideFixed())
            .forEach(node -> failOnNorthSouthPortConnections(node));
    }
    
    private void failOnNorthSouthPortConnections(LNode lNode) {
        boolean fail = lNode.getPorts().stream()
                .filter(lPort -> PortSide.SIDES_NORTH_SOUTH.contains(lPort.getSide()))
                .anyMatch(lPort -> lPort.getConnectedEdges().iterator().hasNext());
        
        if (fail) {
            fail();
        }
    }

}
