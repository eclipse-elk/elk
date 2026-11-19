/*******************************************************************************
 * Copyright (c) 2023 Kiel University and others.
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
package org.eclipse.elk.alg.radial.intermediate;

import org.eclipse.elk.alg.radial.InternalProperties;
import org.eclipse.elk.alg.radial.options.RadialOptions;
import org.eclipse.elk.core.alg.ILayoutProcessor;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.util.IElkProgressMonitor;
import org.eclipse.elk.graph.ElkEdge;
import org.eclipse.elk.graph.ElkNode;

/**
 * Calculates the angles of outgoing edges so that they can be used as an input by subsequent child 
 * layouts. Only makes sense when used in a top-down layout.
 *
 */
public class EdgeAngleCalculator implements ILayoutProcessor<ElkNode> {

    /**
     * For each of edges connected to the root node we calculate its angle and store that information on the 
     * connected target node. This node can then later use that information as basis to align its own layout
     * to the incoming edge. Because this sets an option on child nodes, this is only useful when laying the
     * graph out in a top-down manner (or possibly in multiple layout runs). 
     */
    @Override
    public void process(ElkNode graph, IElkProgressMonitor progressMonitor) {

        ElkNode root = graph.getProperty(InternalProperties.ROOT_NODE);
        for (ElkEdge edge : root.getOutgoingEdges()) {
            
            KVector start = new KVector(edge.getSections().get(0).getStartX(), edge.getSections().get(0).getStartY());
            KVector end = new KVector(edge.getSections().get(0).getEndX(), edge.getSections().get(0).getEndY());
            
            KVector edgeVector = KVector.diff(end, start);
            double angle = Math.atan2(edgeVector.y, edgeVector.x);
            
            edge.getTargets().get(0).setProperty(RadialOptions.ROTATION_TARGET_ANGLE, angle);
        }

    }

}
