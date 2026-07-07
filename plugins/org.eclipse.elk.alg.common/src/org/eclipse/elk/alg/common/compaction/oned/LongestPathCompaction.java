/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.common.compaction.oned;

import java.util.Queue;

import org.eclipse.elk.core.options.Direction;

import com.google.common.collect.Lists;

/**
 * Compacts a previously calculated constraint graph (
 * {@link org.eclipse.elk.alg.layered.compaction.oned.CGraph CGraph}) by using a technique similar
 * to a longest path layering.
 * 
 * The algorithm evaluates the {@link CGroup#reposition} flag.
 */
public class LongestPathCompaction implements ICompactionAlgorithm {

    @Override
    public void compact(final OneDimensionalCompactor compactor) {

        // calculating the left-most position of any element 
        // this will be our starting point for the compaction
        double minStartPos = Double.POSITIVE_INFINITY;
        for (CNode cNode : compactor.cGraph.cNodes) {
            minStartPos = Math.min(minStartPos, 
                                   cNode.cGroup.reference.hitbox.x + cNode.cGroupOffset.x);
        }
        
        // finding the sinks of the constraint graph
        Queue<CGroup> sinks = Lists.newLinkedList();
        for (CGroup group : compactor.cGraph.cGroups) {
            group.startPos = minStartPos;
            if (group.outDegree == 0) {
                sinks.add(group);
            }
        }
        // process sinks until every node in the constraint graph was handled
        while (!sinks.isEmpty()) {
            
            CGroup group = sinks.poll();
            
            // record the movement of this group during the current compaction
            // this has to be recorded _before_ the nodes' positions are updated
            // and care has to be taken about the compaction direction. In certain 
            // scenarios nodes may move "back-and-forth". To detect this, we associate
            // a negative delta with two of the compaction directions.
            double diff = group.reference.hitbox.x;
            
            // ------------------------------------------
            // #1 final positions for this group's nodes
            // ------------------------------------------
            for (CNode node : group.cNodes) {
                // CNodes can be locked in place to avoid pulling clusters apart
                double suggestedX = group.startPos + node.cGroupOffset.x;
                if (!compactor.isLocked(group, compactor.direction)
                        // does the "fixed" position violate the constraints?
                        || (node.hitbox.x < suggestedX)) {
                    node.startPos = suggestedX;
                } else {
                    // leave the node where it was!
                    node.startPos = node.hitbox.x;
                }
            }
            
            diff -= group.reference.startPos;
            
            group.delta += diff;
            if (compactor.direction == Direction.RIGHT || compactor.direction == Direction.DOWN) {
                group.deltaNormalized += diff;
            } else {
                group.deltaNormalized -= diff;
            }
            
            // ---------------------------------------------------
            // #2 propagate start positions to constrained groups
            // ---------------------------------------------------
            for (CNode node : group.cNodes) {
                for (CNode incNode : node.constraints) {
                    // determine the required spacing
                    double spacing;
                    if (compactor.direction.isHorizontal()) {
                        spacing = compactor.spacingsHandler.getHorizontalSpacing(node, incNode);
                    } else {
                        spacing = compactor.spacingsHandler.getVerticalSpacing(node, incNode);
                    }
                    
                    incNode.cGroup.startPos = Math.max(incNode.cGroup.startPos, 
                                                   node.startPos + node.hitbox.width + spacing 
                                                   // respect the other group's node's offset
                                                   - incNode.cGroupOffset.x);
                    
                    // whether the node's current position should be preserved
                    if (compactor.isLocked(incNode, compactor.direction)) {
                        incNode.cGroup.startPos =
                                Math.max(incNode.cGroup.startPos, incNode.hitbox.x
                                        - incNode.cGroupOffset.x);
                    }
                    
                    incNode.cGroup.outDegree--;
                    if (incNode.cGroup.outDegree == 0) {
                       sinks.add(incNode.cGroup); 
                    }
                }
            }
        }
        
        // ------------------------------------------------------
        // #3 setting hitbox positions to new starting positions
        // ------------------------------------------------------
        for (CNode cNode : compactor.cGraph.cNodes) {
            cNode.hitbox.x = cNode.startPos;
        }
    }

}
