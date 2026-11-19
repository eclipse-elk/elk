/*******************************************************************************
 * Copyright (c) 2022 Kiel University and others.
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
package org.eclipse.elk.alg.topdownpacking;

import java.util.List;

import org.eclipse.elk.core.alg.ILayoutPhase;
import org.eclipse.elk.core.alg.LayoutProcessorConfiguration;
import org.eclipse.elk.core.math.ElkPadding;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.util.IElkProgressMonitor;
import org.eclipse.elk.graph.ElkNode;

/**
 * If there is additional space in the last row, the nodes are widened by an equal factor and shifted to
 * the right to make use of that space.
 *
 */
public class BottomRowEqualWhitespaceEliminator implements ILayoutPhase<TopdownPackingPhases, GridElkNode> {

    /**
     * {@inheritDoc}
     */
    @Override
    public void process(GridElkNode layoutGraph, IElkProgressMonitor progressMonitor) {
        // Start progress monitor
        progressMonitor.begin("Whitespace elimination", 1);
        progressMonitor.log("Whitespace elimination began for node " + layoutGraph.getIdentifier());
        
        if (layoutGraph.getWidth() == 0 || layoutGraph.getColumns() == 0) {
            progressMonitor.log("Parent node has no width, skipping phase");
            progressMonitor.done();
            return;
        }
        
        ElkPadding padding = layoutGraph.getProperty(CoreOptions.PADDING);
        
        // for each row check whether there is white space, if there is expand and shift nodes
        // TODO: optimization chance here, if parent width matches predicted width, we can skip directly to the last row
        for (int i = 0; i < layoutGraph.getRows(); i++) {
            List<ElkNode> row = layoutGraph.getRow(i);
            // check for whitespace next to last node
            ElkNode last = null;
            int lastIndex = row.size();
            while (last == null) {
                last = row.get(--lastIndex);
            }
            double rightBorder = last.getX() + last.getWidth();
            
            if (rightBorder + padding.right < layoutGraph.getWidth()) {
                progressMonitor.log("Eliminate white space in row " + i);
                double extraSpace = layoutGraph.getWidth() - (rightBorder + padding.right);
                double extraSpacePerNode = extraSpace / (lastIndex + 1);
                double accumulatedShift = 0;
                // go through all nodes in row, shift and enlargen them
                for (int j = 0; j <= lastIndex; j++) {
                    ElkNode node = row.get(j);
                    node.setX(node.getX() + accumulatedShift);
                    node.setWidth(node.getWidth() + extraSpacePerNode);
                    accumulatedShift += extraSpacePerNode;
                }
            }
        }

        // check whether there is vertical white space
        List<ElkNode> col = layoutGraph.getColumn(0);
        // check for whitespace below first column and expand all columns accordingly, 
        // to prevent expanding into same space as horizontal expansion
        ElkNode last = col.get(col.size() - 1);
        double bottomBorder = last.getY() + last.getHeight();
        double extraSpace = layoutGraph.getHeight() - (bottomBorder + padding.bottom);
        double extraSpacePerNode = extraSpace / (col.size() + 1);
        double accumulatedShift = 0;
        progressMonitor.log("Eliminate vertical white space");
        if (bottomBorder + padding.bottom < layoutGraph.getHeight()) {
            for (int i = 0; i < layoutGraph.getRows(); i++) {
                List<ElkNode> row = layoutGraph.getRow(i);
                // go through all nodes in row, shift and enlarge them
                for (ElkNode node : row) {
                    if (node == null) {
                        break;
                    }
                    node.setY(node.getY() + accumulatedShift);
                    node.setHeight(node.getHeight() + extraSpacePerNode);
                    accumulatedShift += extraSpacePerNode;
                }
            }
        }
        
        progressMonitor.logGraph(layoutGraph, "Graph after whitespace elimination");
        progressMonitor.done();
        
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public LayoutProcessorConfiguration<TopdownPackingPhases, GridElkNode> getLayoutProcessorConfiguration(
            GridElkNode graph) {
        return LayoutProcessorConfiguration.<TopdownPackingPhases, GridElkNode>create();
    }

}
