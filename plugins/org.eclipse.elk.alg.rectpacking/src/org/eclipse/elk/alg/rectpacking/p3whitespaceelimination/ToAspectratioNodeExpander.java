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

import org.eclipse.elk.alg.rectpacking.options.InternalProperties;
import org.eclipse.elk.alg.rectpacking.options.RectPackingOptions;
import org.eclipse.elk.core.util.IElkProgressMonitor;
import org.eclipse.elk.graph.ElkNode;

/**
 * Eliminates the whitespace in the placement of the child nodes by increasing the size of the children equally.
 * 
 * <dl>
 *   <dt>Precondition:</dt>
 *     <dd>The graph is divided into rows, stacks, blocks and subrows.</dd>
 *   <dt>Postcondition:</dt>
 *     <dd>The whitespace is eliminated</dd>
 *     <dd>The drawing has the desired aspect ratio.</dd>
 * </dl>
 */
public class ToAspectratioNodeExpander extends EqualWhitespaceEliminator {

    @Override
    public void process(ElkNode graph, IElkProgressMonitor progressMonitor) {
        progressMonitor.begin("To Aspect Ratio Whitesapce Eliminator", 1);
        double width = graph.getProperty(InternalProperties.DRAWING_WIDTH);
        double height = graph.getProperty(InternalProperties.DRAWING_HEIGHT);
        double desiredAspectRatio = graph.getProperty(RectPackingOptions.ASPECT_RATIO);
        double additionalHeight = graph.getProperty(InternalProperties.ADDITIONAL_HEIGHT);
        double aspectRatio = width / height;
        if (aspectRatio < desiredAspectRatio) {
            width = height * desiredAspectRatio;
            graph.setProperty(InternalProperties.DRAWING_WIDTH, width);
        } else {
            additionalHeight += (width / desiredAspectRatio) - height;
            graph.setProperty(InternalProperties.ADDITIONAL_HEIGHT, additionalHeight);
            graph.setProperty(InternalProperties.DRAWING_HEIGHT, height + additionalHeight);
        }
        super.process(graph, progressMonitor);
        progressMonitor.done();
    }

}
