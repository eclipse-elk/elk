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
package org.eclipse.elk.alg.mrtree.intermediate;

import org.eclipse.elk.alg.mrtree.graph.TGraph;
import org.eclipse.elk.alg.mrtree.graph.TNode;
import org.eclipse.elk.alg.mrtree.options.InternalProperties;
import org.eclipse.elk.alg.mrtree.options.MrTreeOptions;
import org.eclipse.elk.core.alg.ILayoutProcessor;
import org.eclipse.elk.core.options.Direction;
import org.eclipse.elk.core.util.IElkProgressMonitor;

/**
 * @author sdo
 *
 */
public class DirectionProcessor implements ILayoutProcessor<TGraph> {
    /**
     * {@inheritDoc}
     */
    public void process(final TGraph tGraph, final IElkProgressMonitor progressMonitor) {
        progressMonitor.begin("Process directions", 1);
        Direction d = tGraph.getProperty(MrTreeOptions.DIRECTION);
        
        /**
         * Swap coordinates for Direction option
         * Because this modifies the coordinate values it has to be used between Phase 3 and the NodePositionProcessor
         * to work.
         * Modifying those values is easier here because they are at the center of where the nodes will be and we want 
         * to rotate the nodes around their centers
         */
        if (d != Direction.DOWN) {
            for (TNode n : tGraph.getNodes()) {
                int x = n.getProperty(InternalProperties.XCOOR);
                int y = n.getProperty(InternalProperties.YCOOR);
                
                switch (d) {
                case UP:
                    y *= -1;
                    break;
                case RIGHT:
                    int tmp = x;
                    x = y;
                    y = tmp;
                    break;
                case LEFT:
                    int tmp2 = x;
                    x = -y;
                    y = tmp2;
                    break;
                }
                
                n.setProperty(InternalProperties.XCOOR, x);
                n.setProperty(InternalProperties.YCOOR, y);
            }
        }
        
        progressMonitor.done();
    }
}
