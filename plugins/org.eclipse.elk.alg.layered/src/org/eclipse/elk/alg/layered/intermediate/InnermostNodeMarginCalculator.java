/*******************************************************************************
 * Copyright (c) 2010, 2019 Kiel University and others.
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
package org.eclipse.elk.alg.layered.intermediate;

import org.eclipse.elk.alg.common.nodespacing.NodeDimensionCalculation;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LGraphAdapters;
import org.eclipse.elk.core.alg.ILayoutProcessor;
import org.eclipse.elk.core.util.IElkProgressMonitor;

/**
 * Computes and sets the inner node margins.
 * 
 * <p>Node margins are the space around a node that must not overlap with other diagram elements. The space can be
 * thought of as being organized like layers around the node. From inner to outer layers, they include space for the
 * following elements:</p>
 * 
 * <ul>
 *   <li>Ports and port labels.</li>
 *   <li>Self loops</li>
 *   <li>Attached comment boxes</li>
 *   <li>End labels of incident edges</li>
 * </ul>
 * 
 * <p>This processor only computes the space required for ports and port labels. The margins are extended by the
 * {@link SelfLoopRouter}, the {@link CommentNodeMarginCalculator}, and the {@link EndLabelPreprocessor}.
 * 
 * <dl>
 *   <dt>Preconditions:</dt>
 *     <dd>A layered graph.</dd>
 *     <dd>Ports have fixed port positions.</dd>
 *     <dd>Node and port labels have fixed positions.</dd>
 *   <dt>Postcondition:</dt>
 *     <dd>The node margins are properly set to form a bounding box around the node and its ports and labels. Self
 *         loops, comment boxes, and end labels of edges are not included in the margins.</dd>
 *   <dt>Slots:</dt>
 *     <dd>Before phase 4.</dd>
 *   <dt>Same-slot dependencies:</dt>
 *     <dd>{@link LabelAndNodeSizeProcessor}</dd>
 * </dl>
 */
public final class InnermostNodeMarginCalculator implements ILayoutProcessor<LGraph> {

    @Override
    public void process(final LGraph layeredGraph, final IElkProgressMonitor monitor) {
        monitor.begin("Node margin calculation", 1);

        // Calculate the margins using ELK's utility methods. What is not included in the margins yet are is space
        // required for self loops and comment boxes. We will deal with all of those later.
        NodeDimensionCalculation.getNodeMarginCalculator(LGraphAdapters.adapt(layeredGraph, false))
                .excludeEdgeHeadTailLabels().process();

        monitor.done();
    }

}
