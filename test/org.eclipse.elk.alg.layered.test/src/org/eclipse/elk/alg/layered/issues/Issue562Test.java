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
package org.eclipse.elk.alg.layered.issues;

import static org.junit.Assert.fail;

import org.eclipse.elk.core.LayoutConfigurator;
import org.eclipse.elk.core.RecursiveGraphLayoutEngine;
import org.eclipse.elk.core.UnsupportedConfigurationException;
import org.eclipse.elk.core.data.LayoutAlgorithmResolver;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.util.ElkUtil;
import org.eclipse.elk.core.util.NullElkProgressMonitor;
import org.eclipse.elk.graph.ElkEdge;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.ElkPort;
import org.eclipse.elk.graph.util.ElkGraphUtil;
import org.junit.Test;

/**
 * Test for issue 562.
 */
public class Issue562Test {
    
    @Test
    public void test() {
        // Create the basic graph structure
        ElkNode graph = ElkGraphUtil.createGraph();
        ElkNode n1 = ElkGraphUtil.createNode(graph);
        ElkPort p1 = ElkGraphUtil.createPort(n1);
        ElkPort p2 = ElkGraphUtil.createPort(n1);
        ElkEdge e = ElkGraphUtil.createSimpleEdge(p1, p2);
        
        // Create a layout configurator to apply options, just as it would happen in the DiagramLayoutEngine
        LayoutConfigurator config = new LayoutConfigurator();
        
        config.configure(n1).setProperty(CoreOptions.INSIDE_SELF_LOOPS_ACTIVATE, true);
        config.configure(e).setProperty(CoreOptions.INSIDE_SELF_LOOPS_YO, true);
        
        // Apply the configurator and a layout algorithm resolver, just like the DiagramLayoutEngine does
        ElkUtil.applyVisitors(graph, config, new LayoutAlgorithmResolver());
        
        // Apply layout. We don't want an UnsupportedConfigurationException to be thrown
        try {
            new RecursiveGraphLayoutEngine().layout(graph, new NullElkProgressMonitor());
        } catch (UnsupportedConfigurationException ex) {
            fail(ex.toString());
        }
    }

}
