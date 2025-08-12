/*******************************************************************************
 * Copyright (c) 2023 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 * 
 * SPDX-License-Identifier: EPL-2.0 
 *******************************************************************************/
package org.eclipse.elk.alg.vertiflex.p1yplacement;

import org.eclipse.elk.alg.vertiflex.InternalProperties;
import org.eclipse.elk.alg.vertiflex.VertiFlexLayoutPhases;
import org.eclipse.elk.alg.vertiflex.options.VertiFlexOptions;
import org.eclipse.elk.core.alg.ILayoutPhase;
import org.eclipse.elk.core.alg.LayoutProcessorConfiguration;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.util.IElkProgressMonitor;
import org.eclipse.elk.graph.ElkNode;

/**
 * Node placer to position nodes vertically. Nodes that have a vertical constraint are placed according to that 
 * constraint and other nodes are positioned automatically according to their position in the tree.
 *
 */
public class NodeYPlacer implements ILayoutPhase<VertiFlexLayoutPhases, ElkNode> {

    private double layerDistance;
    private double nodeNodeSpacing;
    
    @Override
    public void process(final ElkNode graph, final IElkProgressMonitor progressMonitor) {

        progressMonitor.begin("YPlacer", 1);
        
        layerDistance = graph.getProperty(VertiFlexOptions.LAYER_DISTANCE);
        nodeNodeSpacing = graph.getProperty(CoreOptions.SPACING_NODE_NODE);

        if (!graph.getChildren().isEmpty()) {
            ElkNode parent = graph.getProperty(InternalProperties.ROOT_NODE);
            setYLevels(parent, 0.0);
        }
        
        progressMonitor.done();
    }
    
    /**
     * A Method to set the absolute Y coordinates of the nodes.
     * Uses a default distance between the nodes.
     * @param node and it's children are getting updated Y-coords.
     * @param minHeight: Node gets minimum this height.
     */
    private void setYLevels(final ElkNode node, double minHeight) {
        if (node.hasProperty(VertiFlexOptions.VERTICAL_CONSTRAINT)) {
            minHeight = node.getProperty(VertiFlexOptions.VERTICAL_CONSTRAINT);
        }
        node.setY(minHeight);
        double newMinHeight = minHeight + layerDistance + node.getHeight() 
            + Math.max(node.getProperty(CoreOptions.MARGINS).bottom, nodeNodeSpacing);
        for (int i = 0; i < node.getOutgoingEdges().size(); i++) {
            ElkNode child = (ElkNode) node.getOutgoingEdges().get(i).getTargets().get(0);
            setYLevels(child, newMinHeight);
        }
    }

    @Override
    public LayoutProcessorConfiguration<VertiFlexLayoutPhases, ElkNode> 
                getLayoutProcessorConfiguration(final ElkNode graph) {
        return null;
    }
    
}
