/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
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
package org.eclipse.elk.alg.spore.p3execution;

import org.eclipse.elk.alg.common.Tree;
import org.eclipse.elk.alg.common.spore.DepthFirstCompaction;
import org.eclipse.elk.alg.common.spore.InternalProperties;
import org.eclipse.elk.alg.common.spore.Node;
import org.eclipse.elk.alg.common.utils.SVGImage;
import org.eclipse.elk.alg.spore.SPOrEPhases;
import org.eclipse.elk.alg.spore.graph.Graph;
import org.eclipse.elk.core.alg.ILayoutPhase;
import org.eclipse.elk.core.alg.LayoutProcessorConfiguration;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.util.ElkUtil;
import org.eclipse.elk.core.util.IElkProgressMonitor;

/**
 * This phase compacts the {@link Graph} by shrinking the previously created {@link Tree}, i.e. shortening the
 * edges of the tree to reduce underlap of the nodes.
 */
public class ShrinkTreeCompactionPhase implements ILayoutPhase<SPOrEPhases, Graph> {
    
    @Override
    public LayoutProcessorConfiguration<SPOrEPhases, Graph> getLayoutProcessorConfiguration(final Graph graph) {
        return LayoutProcessorConfiguration.<SPOrEPhases, Graph>create();
    }

    @Override
    public void process(final Graph graph, final IElkProgressMonitor progressMonitor) {
        
        progressMonitor.begin("Shrinking tree compaction", 1);
        
        if (graph.getProperty(InternalProperties.DEBUG_SVG)) {
            debugOut(graph.tree);
            DepthFirstCompaction.compact(graph.tree, graph.orthogonalCompaction,
                    ElkUtil.debugFolderPath("spore") + "60compaction");
        } else {
            DepthFirstCompaction.compact(graph.tree, graph.orthogonalCompaction);
        }
        
        progressMonitor.done();
    }
    
    //--- debug utils --------------------------------------
    
    private void debugOut(final Tree<Node> tree) {
        SVGImage svg = new SVGImage(ElkUtil.debugFolderPath("spore") + "30Tree");
        svg.clear();
        // CHECKSTYLEOFF MagicNumber
        svg.addCircle(tree.node.vertex.x, tree.node.vertex.y, 10, "fill=\"lime\"");
        draw(tree, svg);
        svg.save();
        svg.debug = false;
    }
    
    private void draw(final Tree<Node> t, final SVGImage svg) {
        svg.addRect(t.node.rect, "fill=\"none\" stroke=\"black\"");
        t.children.forEach(c -> {
            svg.addLine(t.node.vertex.x, t.node.vertex.y, c.node.vertex.x, c.node.vertex.y, "stroke=\"blue\"");
            KVector cv = t.node.vertex.clone().sub(c.node.vertex);
            cv.scaleToLength(t.node.distance(c.node, cv));
            svg.addLine(c.node.vertex.x, c.node.vertex.y, c.node.vertex.x + cv.x, c.node.vertex.y + cv.y, 
                    "stroke=\"orange\"");
            draw(c, svg);
        });
    }
}
