/*******************************************************************************
 * Copyright (c) 2016 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.core.util.selection;

import java.util.Collections;
import java.util.Iterator;

import org.eclipse.elk.graph.ElkConnectableShape;
import org.eclipse.elk.graph.ElkEdge;
import org.eclipse.elk.graph.ElkGraphElement;
import org.eclipse.elk.graph.ElkPort;

import com.google.common.collect.Iterators;
import com.google.common.collect.Sets;

/**
 * The default {@link SelectionIterator} for usage in
 * {@link org.eclipse.elk.core.util.ElkUtil#getConnectedElements(ElkEdge, boolean) ElkUtil.getConnectedElements(...)}.
 * The iterator follows all edges transitively connected to the initial edge. Optionally the iterator can also include
 * the corresponding ports in the selection.
 * 
 * <p>
 * MIGRATE This iterator currently does not support hyperedges.
 * </p>
 */
public class DefaultSelectionIterator extends SelectionIterator {

    private static final long serialVersionUID = 2945508835051993888L;

    private boolean addPorts;

    private boolean followEdgeDirection;

    /**
     * Creates a new iterator which can optionally include ports in its selection and can be
     * configured to either iterate towards the tail or the head of the selected {@link ElkEdge}.
     * 
     * @param edge
     *            the edge to start with
     * @param addPorts
     *            flag to determine whether ports should be included in the selection
     * @param followEdgeDirection
     *            flag whether the iterator should iterate towards the head or the tail of the edge
     */
    public DefaultSelectionIterator(final ElkEdge edge, final boolean addPorts, final boolean followEdgeDirection) {
        super(edge);
        this.addPorts = addPorts;
        this.followEdgeDirection = followEdgeDirection;
    }

    @Override
    protected Iterator<? extends ElkGraphElement> getChildren(final Object object) {
        // Ensure that the visited set is properly initialized
        if (visited == null) {
            visited = Sets.newHashSet();
        }
        
        if (object instanceof ElkEdge) {
            ElkEdge edge = (ElkEdge) object;
            
            // Get the 
            ElkConnectableShape connectedShape =
                    followEdgeDirection ? edge.getTargets().get(0) : edge.getSources().get(0);
            
            if (!(connectedShape instanceof ElkPort) || visited.contains((ElkPort) connectedShape)) {
                // return an empty iterator if no target/source port is configured
                // or if the target/source port has been visited already, in order
                // to break infinite loops
                return Collections.emptyIterator();
            }
            
            ElkPort port = (ElkPort) connectedShape;
            visited.add(port);

            // for each object visited by this iterator check all the edges connected to
            // 'port' and visit those edges satisfying the criterion stated above
            Iterator<ElkEdge> resultEdges =
                    (followEdgeDirection ? port.getOutgoingEdges() : port.getIncomingEdges()).iterator();

            // If the port should be added to the selection, add it to the result set
            if (addPorts) {
                Iterator<ElkGraphElement> portIterator =
                        Iterators.singletonIterator((ElkGraphElement) port);
                return Iterators.concat(portIterator, resultEdges);
            } else {
                return resultEdges;
            }
        } else {
            return Collections.emptyIterator();
        }
    }
}
