/*******************************************************************************
 * Copyright (c) 2023 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 * 
 * SPDX-License-Identifier: EPL-2.0 
 *******************************************************************************/
package org.eclipse.elk.alg.yconstree.p4absolute;

import org.eclipse.elk.alg.yconstree.InternalProperties;
import org.eclipse.elk.alg.yconstree.YconstreeLayoutPhases;
import org.eclipse.elk.core.alg.ILayoutPhase;
import org.eclipse.elk.core.alg.LayoutProcessorConfiguration;
import org.eclipse.elk.core.util.IElkProgressMonitor;
import org.eclipse.elk.graph.ElkNode;

/**
 * Computes absolute x coordinates from the previously computed relative coordinates.
 *
 */
public class AbsoluteXPlacer implements ILayoutPhase<YconstreeLayoutPhases, ElkNode> {
    
    private IElkProgressMonitor myProgressMonitor;
    
    @Override
    public void process(final ElkNode graph, final IElkProgressMonitor progressMonitor) {
        myProgressMonitor = progressMonitor;
        myProgressMonitor.begin("AbsolutPlacer", 1);
        
        try {
            if (!graph.getChildren().isEmpty()) {
                ElkNode parent = graph.getProperty(InternalProperties.ROOT_NODE);
                
                // first, move the root
                parent.setX(parent.getX() - findMinimalX(parent));
                // a little offset
                parent.setX(parent.getX() + 10.0); // TODO remove magic number
                // now we update the whole tree to absolute X 
                absoluteTreeCoords(parent);
            }
        } catch (Exception e) {
            // TODO properly handle exception
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
        
        myProgressMonitor.done();
        
    }
    
    
    private double findMinimalX(final ElkNode tree) {
        int numOfChildren = tree.getOutgoingEdges().size();
        if (numOfChildren == 0) {
            return tree.getX();
        } else {
            double minSubtreeX = 0.0;
            double testX = 0.0;
            for (int i = 0; i < numOfChildren; i++) {
                testX = findMinimalX((ElkNode) tree.getOutgoingEdges().get(i).getTargets().get(0));
                minSubtreeX = (testX < minSubtreeX) ? testX : minSubtreeX;
            }
            return minSubtreeX + tree.getX();
        }
    }
    
    private void absoluteTreeCoords(final ElkNode tree) {
        int numOfChildren = tree.getOutgoingEdges().size();
        // a little offset
        tree.setY(tree.getY() + 10.0); // TODO remove magic number
        if (numOfChildren > 0) {
            ElkNode child;
            for (int i = 0; i < numOfChildren; i++) {
                child = (ElkNode) tree.getOutgoingEdges().get(i).getTargets().get(0);
                child.setX(child.getX() + tree.getX());
                absoluteTreeCoords(child);
            }
        }
    }
    
    
    @Override
    public LayoutProcessorConfiguration<YconstreeLayoutPhases, ElkNode> getLayoutProcessorConfiguration(
            final ElkNode graph) {
        return null;
    }

}
