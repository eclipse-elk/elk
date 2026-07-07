/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.spore.p2processingorder;

import org.eclipse.elk.alg.common.ICostFunction;
import org.eclipse.elk.alg.common.NaiveMinST;
import org.eclipse.elk.alg.common.Tree;
import org.eclipse.elk.alg.common.spore.InternalProperties;
import org.eclipse.elk.alg.spore.graph.Graph;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.util.ElkUtil;
import org.eclipse.elk.core.util.IElkProgressMonitor;

/**
 * This phase finds a maximum spanning tree for the edges of the {@link Graph} with the 
 * specified root {@link Graph#preferredRoot}.
 */
public class MaxSTPhase extends MinSTPhase {
    @Override
    public void process(final Graph graph, final IElkProgressMonitor progressMonitor) {
        progressMonitor.begin("Maximum spanning tree construction", 1);
        
        // inverted cost function
        ICostFunction invertedCF = e -> {
            return -graph.costFunction.cost(e);
        };
        
        KVector root;
        if (graph.preferredRoot != null) {
            root = graph.preferredRoot.vertex;
        } else {
            root = graph.vertices.get(0).vertex;
        }
        
        Tree<KVector> tree;
        if (graph.getProperty(InternalProperties.DEBUG_SVG)) {
            tree = NaiveMinST.createSpanningTree(graph.tEdges, root, invertedCF, 
                    ElkUtil.debugFolderPath("spore") + "20minst");
        } else {
            tree = NaiveMinST.createSpanningTree(graph.tEdges, root, invertedCF);
        }
        
        // convert result to a Tree that can be used in the execution phase
        convert(tree, graph);
        
        progressMonitor.done();
    }
}
