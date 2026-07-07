/*******************************************************************************
 * Copyright (c) 2016 Kiel University and others.
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
package org.eclipse.elk.shared.test;

import org.eclipse.elk.alg.layered.LayeredLayoutProvider;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.test.GraphTestUtils;
import org.eclipse.elk.core.AbstractLayoutProvider;
import org.eclipse.elk.core.data.LayoutMetaDataService;
import org.eclipse.elk.core.options.Direction;
import org.eclipse.elk.core.options.PortConstraints;
import org.eclipse.elk.core.options.PortSide;
import org.eclipse.elk.core.options.SizeConstraint;
import org.eclipse.elk.core.util.BasicProgressMonitor;
import org.eclipse.elk.core.util.IElkProgressMonitor;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.ElkPort;
import org.junit.Assert;
import org.junit.BeforeClass;
import org.junit.Test;

/**
 * Test and demonstration of 'Using Layout Algorithms Directly'.
 */
public class DirectLayoutTest {
    
    @BeforeClass
    public static void setup() {
        LayoutMetaDataService.initElkReflect();
    }

    /**
     * Test a plain Java layout using the ELK Layered algorithm.
     */
    @Test
    public void testELKLayered() {
        // create a KGraph for layout
        ElkNode parentNode = GraphTestUtils.createSimpleGraph();

        // add layout options to the elements of the graph
        addLayeredOptions(parentNode);

        // create a progress monitor
        IElkProgressMonitor progressMonitor = new BasicProgressMonitor();

        // create the layout provider
        AbstractLayoutProvider layoutProvider = new LayeredLayoutProvider();
        layoutProvider.initialize(null);

        // perform layout on the created graph
        layoutProvider.layout(parentNode, progressMonitor);

        // execute a trivial junit test
        for (ElkNode child : parentNode.getChildren()) {
            if (child.getX() <= 0 || child.getY() <= 0) {
                Assert.fail("Not all node coordinates have been set properly.");
            }
        }
        
        // output layout information
        printLayoutInfo(parentNode, progressMonitor);
    }

    /**
     * Adds layout options to the elements of the given parent node.
     * 
     * @param parentNode
     *            parent node representing a graph
     */
    private static void addLayeredOptions(final ElkNode parentNode) {
        // set layout direction to horizontal
        parentNode.setProperty(LayeredOptions.DIRECTION, Direction.RIGHT);

        // add options for the child nodes
        for (ElkNode childNode : parentNode.getChildren()) {
            // CHECKSTYLEOFF MagicNumber
            // set some width and height for the child
            childNode.setWidth(30.0f);
            childNode.setHeight(30.0f);
            // set fixed size for the child
            childNode.setProperty(LayeredOptions.NODE_SIZE_CONSTRAINTS, SizeConstraint.fixed());
            // set port constraints to fixed port positions
            childNode.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_POS);

            // add options for the ports
            int i = 0;
            for (ElkPort port : childNode.getPorts()) {
                i++;
                // set position and side
                port.setY(i * 30.0 / (childNode.getPorts().size() + 1));
                if (childNode.getLabels().get(0).getText().equals("node1")) {
                    port.setX(30.0);
                    port.setProperty(LayeredOptions.PORT_SIDE, PortSide.EAST);
                } else {
                    port.setX(0.0);
                    port.setProperty(LayeredOptions.PORT_SIDE, PortSide.WEST);
                }
            }
        }
    }

    /**
     * Outputs layout information on the console.
     * 
     * @param parentNode
     *            parent node representing a graph
     * @param progressMonitor
     *            progress monitor for the layout run
     */
    private void printLayoutInfo(final ElkNode parentNode, final IElkProgressMonitor progressMonitor) {
        // print execution time of the algorithm run
        System.out.println("Execution time: " + progressMonitor.getExecutionTime() * 1000 + " ms");

        // print position of each node
        for (ElkNode childNode : parentNode.getChildren()) {
            System.out.println(childNode.getLabels().get(0).getText() + ": x = " + childNode.getX() 
                               + ", y = " + childNode.getY());
        }
    }
}
