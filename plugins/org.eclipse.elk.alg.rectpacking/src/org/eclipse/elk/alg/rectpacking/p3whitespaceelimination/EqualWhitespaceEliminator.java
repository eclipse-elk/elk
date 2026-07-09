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
package org.eclipse.elk.alg.rectpacking.p3whitespaceelimination;

import org.eclipse.elk.alg.rectpacking.RectPackingLayoutPhases;
import org.eclipse.elk.alg.rectpacking.options.InternalProperties;
import org.eclipse.elk.alg.rectpacking.options.RectPackingOptions;
import org.eclipse.elk.core.UnsupportedConfigurationException;
import org.eclipse.elk.core.alg.ILayoutPhase;
import org.eclipse.elk.core.alg.LayoutProcessorConfiguration;
import org.eclipse.elk.core.util.IElkProgressMonitor;
import org.eclipse.elk.graph.ElkNode;

/**
 * Eliminates the whitespace in the placement of the child nodes by increasing the size of the children equally.
 * 
 * <dl>
 *   <dt>Precondition:</dt>
 *     <dd>The graph is divided into rows, stacks, blocks and subrows.</dd>
 *   <dt>Postcondition:</dt>
 *     <dd>The whitespace is eliminated and</dd>
 * </dl>
 */
public class EqualWhitespaceEliminator implements ILayoutPhase<RectPackingLayoutPhases, ElkNode> {

    /* (non-Javadoc)
     * @see org.eclipse.elk.core.alg.ILayoutProcessor#process(java.lang.Object, org.eclipse.elk.core.util.IElkProgressMonitor)
     */
    @Override
    public void process(ElkNode graph, IElkProgressMonitor progressMonitor) {
        progressMonitor.begin("Equal Whitespace Eliminator", 1);
        if (graph.hasProperty(InternalProperties.ROWS)) {
            RectangleExpansion.expand(graph.getProperty(InternalProperties.ROWS),
                    graph.getProperty(InternalProperties.DRAWING_WIDTH),
                    graph.getProperty(InternalProperties.ADDITIONAL_HEIGHT),
                    graph.getProperty(RectPackingOptions.SPACING_NODE_NODE));
        } else {
            throw new UnsupportedConfigurationException("The graph does not contain rows.");
        }
        progressMonitor.done();
        
    }

    /* (non-Javadoc)
     * @see org.eclipse.elk.core.alg.ILayoutPhase#getLayoutProcessorConfiguration(java.lang.Object)
     */
    @Override
    public LayoutProcessorConfiguration<RectPackingLayoutPhases, ElkNode> getLayoutProcessorConfiguration(
            ElkNode graph) {
        return null;
    }

}
