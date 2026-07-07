/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.spore.p3execution;

import org.eclipse.elk.alg.common.Tree;
import org.eclipse.elk.alg.common.spore.InternalProperties;
import org.eclipse.elk.alg.common.spore.Node;
import org.eclipse.elk.alg.common.utils.SVGImage;
import org.eclipse.elk.alg.common.utils.Utils;
import org.eclipse.elk.alg.spore.SPOrEPhases;
import org.eclipse.elk.alg.spore.graph.Graph;
import org.eclipse.elk.core.alg.ILayoutPhase;
import org.eclipse.elk.core.alg.LayoutProcessorConfiguration;
import org.eclipse.elk.core.util.ElkUtil;
import org.eclipse.elk.core.util.IElkProgressMonitor;

/**
 * This class implements the GTree algorithm of Nachmanson et al..
 */
public class GrowTreePhase implements ILayoutPhase<SPOrEPhases, Graph> {
    private SVGImage svg;
    private Tree<Node> root;
    
    private boolean overlapsExisted;
    
    @Override
    public LayoutProcessorConfiguration<SPOrEPhases, Graph> getLayoutProcessorConfiguration(final Graph graph) {
        return LayoutProcessorConfiguration.<SPOrEPhases, Graph>create();
    }

    @Override
    public void process(final Graph graph, final IElkProgressMonitor progressMonitor) {
        progressMonitor.begin("Grow Tree", 1);
        
        root = graph.tree;
        
        if (graph.getProperty(InternalProperties.DEBUG_SVG)) {
            svg = new SVGImage(ElkUtil.debugFolderPath("spore") + "40or");
            svg.addGroups("n", "e", "o");
            debugOut();
        } else {
            svg = new SVGImage(null);
        }

        overlapsExisted = false;
        growAt(graph.tree);
        graph.setProperty(InternalProperties.OVERLAPS_EXISTED, overlapsExisted);
        progressMonitor.done();
    }
    
    /**
     * Traverses the Tree recursively and elongates edges where needed.
     * Also sets the overlapsExisted flag to true if any overlaps had to be removed.
     * @param r the spanning tree
     */
    private void growAt(final Tree<Node> r) {
        for (Tree<Node> c : r.children) {
            // update position of the child
            c.node.translate(r.node.vertex.clone().sub(r.node.originalVertex));
            
            // the elongation factor for an edge required to remove overlap
            double t = Utils.overlap(r.node.rect, c.node.rect);
            if (t > 1.0) {
                overlapsExisted = true;
            }
            
            // elongate the edge by factor t to remove overlap
            c.node.setCenterPosition(r.node.vertex.clone().add(
                    c.node.originalVertex.clone().sub(r.node.originalVertex)
                    .scale(t)));
            
            debugOut(r);
            
            growAt(c);
        }
    }
    
    //--- debug utils --------------------------------------
    
    private void debugOut() {
        debugOut(null);
    }
    
    private void debugOut(final Tree<Node> c) {
        svg.clearGroup("n");
        svg.clearGroup("e");
        svg.clearGroup("o");
        drawTree(root, svg);
        svg.g("n").addRect(root.node.rect, "fill=\"blue\" stroke=\"none\" opacity=\"0.2\"");
        if (c != null) {
            svg.g("o").addRect(c.node.rect, "fill=\"red\" stroke=\"none\" opacity=\"0.2\"");
        }
        svg.isave();
    }
    
    private void drawTree(final Tree<Node> t, final SVGImage img) {
        img.g("n").addRect(t.node.rect, "fill=\"none\" stroke=\"black\"");
        t.children.forEach(c -> {
            img.g("e").addLine(t.node.vertex.x, t.node.vertex.y, c.node.vertex.x, c.node.vertex.y, "stroke=\"blue\"");
            drawTree(c, img);
        });
    }
}
