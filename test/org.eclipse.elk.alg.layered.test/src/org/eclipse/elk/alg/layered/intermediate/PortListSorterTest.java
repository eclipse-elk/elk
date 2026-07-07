/*******************************************************************************
 * Copyright (c) 2010, 2019 Kiel University and others.
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

import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.util.Deque;
import java.util.LinkedList;
import java.util.List;

import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.test.framework.LayoutTestRunner;
import org.eclipse.elk.alg.test.framework.annotations.Algorithm;
import org.eclipse.elk.alg.test.framework.annotations.Configurator;
import org.eclipse.elk.alg.test.framework.annotations.DefaultConfiguration;
import org.eclipse.elk.alg.test.framework.annotations.GraphResourceProvider;
import org.eclipse.elk.alg.test.framework.annotations.TestAfterProcessor;
import org.eclipse.elk.alg.test.framework.io.AbstractResourcePath;
import org.eclipse.elk.alg.test.framework.io.FileExtensionFilter;
import org.eclipse.elk.alg.test.framework.io.ModelResourcePath;
import org.eclipse.elk.core.options.PortConstraints;
import org.eclipse.elk.core.options.PortSide;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.ElkPort;
import org.junit.runner.RunWith;

import com.google.common.collect.Lists;

/**
 * Basic tests for the {@link PortListSorter}.
 * 
 * The ports on each side are either ordered according to their port index or,
 * if no index is set, depending on the position of the port. However, the indices for each side are
 * distinct and do not necessarily have to be gapless.
 * 
 * E.g.
 * 
 * <pre>
 *      0 1 2
 *   |--------|0
 *   |        |
 *   |________|7
 *      9  2
 * </pre>
 * 
 * 
 * The test assures that only {@link PortConstraints#FIXED_ORDER} is set to nodes, and assigns an
 * index to all ports. The test runs until the actual port positions are set, as after the
 * {@link PortListSorter} only the order within the nodes port list is adapted, not the actual
 * coordinates.
 */
@RunWith(LayoutTestRunner.class)
@Algorithm(LayeredOptions.ALGORITHM_ID)
@DefaultConfiguration()
public class PortListSorterTest {

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Sources

    @GraphResourceProvider
    public List<AbstractResourcePath> testGraphs() {
        return Lists.newArrayList(
                new ModelResourcePath("realworld/ptolemy/**/").withFilter(new FileExtensionFilter("elkg")));
    }
    

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Configuration
    
    @Configurator
    public void configurePorts(final ElkNode graph) {
        Deque<ElkNode> nodeQueue = new LinkedList<>(graph.getChildren());
        
        while (!nodeQueue.isEmpty()) {
            ElkNode node = nodeQueue.poll();
            nodeQueue.addAll(node.getChildren());
            
            if (node.getProperty(LayeredOptions.PORT_CONSTRAINTS) != PortConstraints.FIXED_ORDER) {
                node.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_ORDER);
                
                int index = 0;
                for (ElkPort port : node.getPorts()) {
                    port.setProperty(LayeredOptions.PORT_INDEX, index++);
                }
            }
        }
    }


    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Tests
    
    //LongEdgeJoiner
    
    @TestAfterProcessor(PortListSorter.class)
    public void testProperPortSideOrder(final Object graph) {
        LGraph lGraph = (LGraph) graph;
        
        for (Layer layer : lGraph) {
            for (LNode node : layer) {
                PortSide lastSide = null;
                for (LPort port : node.getPorts()) {
                    // upon a side change check that the next side is valid
                    if (lastSide != port.getSide()) {
                        if (lastSide == null) {
                            assertTrue(port.getSide() == PortSide.NORTH
                                    || port.getSide() == PortSide.EAST
                                    || port.getSide() == PortSide.SOUTH
                                    || port.getSide() == PortSide.WEST);
                        } else if (lastSide == PortSide.NORTH) {
                            assertTrue(port.getSide() == PortSide.EAST
                                    || port.getSide() == PortSide.SOUTH
                                    || port.getSide() == PortSide.WEST);
                        } else if (lastSide == PortSide.EAST) {
                            assertTrue(port.getSide() == PortSide.SOUTH
                                    || port.getSide() == PortSide.WEST);
                        } else if (lastSide == PortSide.SOUTH) {
                            assertTrue(port.getSide() == PortSide.WEST);
                        } else if (lastSide == PortSide.WEST) {
                            fail();
                        }
                        
                        lastSide = port.getSide();
                    }
                }
            }
        }
    }
    
    /**
     * Waits until port coordinates have been assigned to check whether they reflect the port order.
     */
    @TestAfterProcessor(LongEdgeJoiner.class)
    public void testMonotonicalCoordinates(final Object graph) {
        LGraph lGraph = (LGraph) graph;
        
        for (Layer layer : lGraph) {
            for (LNode node : layer) {
                Double lastPos = null;
                PortSide lastSide = null;

                for (LPort port : node.getPorts()) {
                    // upon side change reset the reference position
                    if (lastSide == null || lastSide != port.getSide()) {
                        if (port.getSide() == PortSide.NORTH || port.getSide() == PortSide.SOUTH) {
                            lastPos = port.getPosition().x;
                        } else if (port.getSide() == PortSide.EAST || port.getSide() == PortSide.WEST) {
                            lastPos = port.getPosition().y;
                        }
                        
                        lastSide = port.getSide();
                        continue;
                    }

                    // check the positions
                    if (port.getSide() == PortSide.NORTH) {
                        assertTrue(lastPos <= port.getPosition().x);
                    } else if (port.getSide() == PortSide.EAST) {
                        assertTrue(lastPos <= port.getPosition().y);
                    } else if (port.getSide() == PortSide.SOUTH) {
                        assertTrue(lastPos >= port.getPosition().x);
                    } else if (port.getSide() == PortSide.WEST) {
                        assertTrue(lastPos >= port.getPosition().y);
                    }
                }
            }
        }
    }
    
}
