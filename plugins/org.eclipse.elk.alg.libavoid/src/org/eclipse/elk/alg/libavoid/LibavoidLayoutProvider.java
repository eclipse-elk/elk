/*******************************************************************************
 * Copyright (c) 2013, 2022 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.alg.libavoid;

import org.eclipse.elk.alg.common.nodespacing.NodeDimensionCalculation;
import org.eclipse.elk.alg.libavoid.server.LibavoidServer;
import org.eclipse.elk.alg.libavoid.server.LibavoidServerPool;
import org.eclipse.elk.core.AbstractLayoutProvider;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.options.Direction;
import org.eclipse.elk.core.options.PortSide;
import org.eclipse.elk.core.util.ElkUtil;
import org.eclipse.elk.core.util.IElkProgressMonitor;
import org.eclipse.elk.core.util.adapters.ElkGraphAdapters;
import org.eclipse.elk.core.util.adapters.ElkGraphAdapters.ElkGraphAdapter;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.ElkPort;

/**
 * A layout provider for KIML that performs layout using the Libavoid connector routing library. See
 * http://www.adaptagrams.org/documentation/ for further information on the library.
 * 
 * @author uru
 */
public class LibavoidLayoutProvider extends AbstractLayoutProvider {

    private LibavoidServerCommunicator comm = new LibavoidServerCommunicator();

    /**
     * {@inheritDoc}
     */
    @Override
    public void layout(final ElkNode parentNode, final IElkProgressMonitor progressMonitor) {

        // Prepare the graph
    	prepareGraph(parentNode);
        ElkGraphAdapter adapter = ElkGraphAdapters.adapt(parentNode);
        NodeDimensionCalculation.sortPortLists(adapter);
        NodeDimensionCalculation.calculateNodeMargins(adapter);
        
        // create an Libavoid server process instance or use an existing one
        LibavoidServer lvServer = LibavoidServerPool.INSTANCE.fetch();
        // send a layout request to the server process and apply the layout
        comm.requestLayout(parentNode, progressMonitor, lvServer);
        // if everything worked well, release the used process instance
        LibavoidServerPool.INSTANCE.release(lvServer);

    }
    
    private void prepareGraph(final ElkNode parentNode) {
    	Direction direction = parentNode.getProperty(CoreOptions.DIRECTION);
    	for (ElkNode node : parentNode.getChildren()) {
    		for (ElkPort port : node.getPorts()) {
    			
    			// Compute missing port sides
    			PortSide side = port.getProperty(CoreOptions.PORT_SIDE);
    			if (side == PortSide.UNDEFINED) {
    	        	side = ElkUtil.calcPortSide(port, direction);
    	        	port.setProperty(CoreOptions.PORT_SIDE, side);
    	        }
    			
    		}
    	}
    }

}
