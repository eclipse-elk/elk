/*******************************************************************************
 * Copyright (c) 2016 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.shared.test;

import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.test.GraphTestUtils;
import org.eclipse.elk.core.RecursiveGraphLayoutEngine;
import org.eclipse.elk.core.data.LayoutMetaDataService;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.util.BasicProgressMonitor;
import org.eclipse.elk.core.util.IElkProgressMonitor;
import org.eclipse.elk.graph.ElkNode;
import org.junit.Test;

import com.google.common.collect.Iterators;
import com.google.common.collect.UnmodifiableIterator;

/**
 * Test and demonstration of 'Plain Java Layout'. 
 */
public class PlainLayoutTest {

    /**
     * Test a plain Java layout on a hierarchical graph using the ELK Layered algorithm.
     */
    @Test
    public void testPlainLayout() {
        // create a hierarchical KGraph for layout
        ElkNode parentNode = GraphTestUtils.createHierarchicalGraph();

        // configure every hierarchical node to use ELK Layered (which would also be the default) 
        getAllElkNodes(parentNode).forEachRemaining(node -> {
            if (!node.getChildren().isEmpty()) {
                node.setProperty(CoreOptions.ALGORITHM, LayeredOptions.ALGORITHM_ID);
            }
        });

        // create a progress monitor
        IElkProgressMonitor progressMonitor = new BasicProgressMonitor();

        // initialize the meta data service with ELK Layered's meta data
        LayoutMetaDataService service = LayoutMetaDataService.getInstance();
        service.registerLayoutMetaDataProviders(new LayeredOptions());

        // instantiate a recursive graph layout engine and execute layout
        RecursiveGraphLayoutEngine layoutEngine = new RecursiveGraphLayoutEngine();
        layoutEngine.layout(parentNode, progressMonitor);

        // output layout information
        printLayoutInfo(parentNode, progressMonitor);
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
        // SUPPRESS CHECKSTYLE NEXT MagicNumber
        System.out.println("Execution time: " + progressMonitor.getExecutionTime() * 1000 + " ms");

        // print position of each node
        getAllElkNodes(parentNode).forEachRemaining(node -> {
            System.out.println(node.getLabels().get(0).getText() + ": x = " + node.getX() + ", y = "
                    + node.getY());
        });
    }
    
    private UnmodifiableIterator<ElkNode> getAllElkNodes(final ElkNode parentNode) {
        return Iterators.filter(parentNode.eAllContents(), ElkNode.class);
    }
}
