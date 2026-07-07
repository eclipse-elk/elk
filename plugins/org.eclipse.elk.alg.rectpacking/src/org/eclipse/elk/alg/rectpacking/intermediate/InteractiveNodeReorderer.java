/*******************************************************************************
 * Copyright (c) 2022 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 * 
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later 
 *******************************************************************************/
package org.eclipse.elk.alg.rectpacking.intermediate;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.eclipse.elk.alg.rectpacking.options.RectPackingOptions;
import org.eclipse.elk.core.alg.ILayoutProcessor;
import org.eclipse.elk.core.util.IElkProgressMonitor;
import org.eclipse.elk.graph.ElkNode;

/**
 * Sorts all child nodes by their desired position.
 * 
 * <dl>
 *   <dt>Precondition:</dt>
 *   <dt>Postcondition:</dt>
 *     <dd>Children are sorted as specified by their desired position in relation to their original order.</dd>
 *   <dt>Slots:</dt>
 *     <dd>Before phase 1.</dd>
 *   <dt>Same-slot dependencies:</dt>
 * </dl>
 */
public class InteractiveNodeReorderer implements ILayoutProcessor<ElkNode> {

    /* (non-Javadoc)
     * @see org.eclipse.elk.core.alg.ILayoutProcessor#process(java.lang.Object, org.eclipse.elk.core.util.IElkProgressMonitor)
     */
    @Override
    public void process(ElkNode graph, IElkProgressMonitor progressMonitor) {
        progressMonitor.begin("Interactive Node Reorderer", 1);
        List<ElkNode> rectangles = graph.getChildren();
        List<ElkNode> fixedNodes = new ArrayList<>();
        for (ElkNode elkNode : rectangles) {
            if (elkNode.hasProperty(RectPackingOptions.DESIRED_POSITION)) {
                fixedNodes.add(elkNode);
            }
        }
        for (ElkNode elkNode : fixedNodes) {
            rectangles.remove(elkNode);
        }
        Collections.sort(fixedNodes, (a, b) -> {
            int positionA = a.getProperty(RectPackingOptions.DESIRED_POSITION);
            int positionB = b.getProperty(RectPackingOptions.DESIRED_POSITION);
            if (positionA == positionB) {
                return -1;
            } else {
                return Integer.compare(positionA, positionB);
            }
        });
        for (ElkNode elkNode : fixedNodes) {
            int position = elkNode.getProperty(RectPackingOptions.DESIRED_POSITION);
            position = Math.min(position, rectangles.size());
            rectangles.add(position, elkNode);
        }

        int index = 0;
        for (ElkNode elkNode: rectangles) {
            elkNode.setProperty(RectPackingOptions.CURRENT_POSITION, index);
            index++;
        }
        progressMonitor.done();
    }

}
