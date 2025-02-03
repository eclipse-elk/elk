/*******************************************************************************
 * Copyright (c) 2023 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 * 
 * SPDX-License-Identifier: EPL-2.0 
 *******************************************************************************/
package org.eclipse.elk.alg.vertiflex.p4edgerouting;

import org.eclipse.elk.alg.vertiflex.InternalProperties;
import org.eclipse.elk.alg.vertiflex.VertiFlexLayoutPhases;
import org.eclipse.elk.core.alg.ILayoutPhase;
import org.eclipse.elk.core.alg.LayoutProcessorConfiguration;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.util.IElkProgressMonitor;
import org.eclipse.elk.graph.ElkEdge;
import org.eclipse.elk.graph.ElkEdgeSection;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.util.ElkGraphUtil;

/**
 * Implements edge routing that includes bendpoints when straight edges would cross nodes.
 *
 */
public class BendEdgeRouter implements ILayoutPhase<VertiFlexLayoutPhases, ElkNode> {

    @Override
    public void process(final ElkNode graph, final IElkProgressMonitor progressMonitor) {
        
        progressMonitor.begin("BendEdgeRouter", 1);
        
        if (!graph.getChildren().isEmpty()) {
            ElkNode parent = graph.getProperty(InternalProperties.ROOT_NODE);
            
            routeEdges(parent);
        }
        
        progressMonitor.done();

    }

    @Override
    public LayoutProcessorConfiguration<VertiFlexLayoutPhases, ElkNode> getLayoutProcessorConfiguration(
            final ElkNode graph) {
        return null;
    }
    
    /** Route the edges with bendpoints. */
    private void routeEdges(final ElkNode node) {
        for (ElkEdge edge : ElkGraphUtil.allOutgoingEdges(node)) {
            ElkNode target = ElkGraphUtil.connectableShapeToNode(edge.getTargets().get(0));
            ElkEdgeSection section = ElkGraphUtil.firstEdgeSection(edge, true, true);
            
            
            double startX = node.getX() + node.getWidth() / 2;
            double startY = node.getY() + node.getHeight();
            double endX = target.getX() + target.getWidth() / 2;
            double endY = target.getY();
            
            section.setStartLocation(startX, startY);
            section.setEndLocation(endX, endY);
            
            double bendheight = target.getProperty(InternalProperties.EDGE_BEND_HEIGHT);
            double epsilon = 0.0001;
            // if the node is low place a bendpoint above it
            if (Math.abs(bendheight 
                    - (endY - target.getParent().getProperty(CoreOptions.SPACING_NODE_NODE) / 2)) > epsilon) {
                ElkGraphUtil.createBendPoint(section, endX, bendheight);
            }
            routeEdges(target);
        }
    }

}
