/*******************************************************************************
 * Copyright (c) 2013 - 2022 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.mrtree.intermediate;

import java.util.ArrayList;

import org.eclipse.elk.alg.mrtree.graph.TGraph;
import org.eclipse.elk.alg.mrtree.graph.TNode;
import org.eclipse.elk.alg.mrtree.options.InternalProperties;
import org.eclipse.elk.core.alg.ILayoutProcessor;
import org.eclipse.elk.core.util.IElkProgressMonitor;

/**
 * A processor connect all roots of a given graph to a super root which then is the new root of the
 * graph.
 * 
 * @author sor
 * @author sgu
 * @author sdo
 */
public class RootProcessor implements ILayoutProcessor<TGraph> {

    private ArrayList<TNode> roots = new ArrayList<TNode>();

    /**
     * {@inheritDoc}
     */
    public void process(final TGraph tGraph, final IElkProgressMonitor progressMonitor) {
        progressMonitor.begin("Find roots", 1);

        /** clear list of roots if processor is reused */
        roots.clear();

        /** find all roots in the graph */
        for (TNode node : tGraph.getNodes()) {
            if (node.getIncomingEdges().isEmpty()) {
                node.setProperty(InternalProperties.ROOT, true);
                roots.add(node);
            }
        }
        /**
         * if there are more than one root add a super root and set the current roots as its
         * children
         */
        switch (roots.size()) {
        case 0:
            assert tGraph.getNodes().isEmpty();
            TNode root = new TNode(0, tGraph, "DUMMY_ROOT");
            root.setProperty(InternalProperties.ROOT, true);
            root.setProperty(InternalProperties.DUMMY, true);
            tGraph.getNodes().add(root);
            break;
            
        case 1:
            // perfect, we already have only one root
            break;

        default:
            TNode superRoot = new TNode(0, tGraph, "SUPER_ROOT");

            for (TNode tRoot : roots) {
                superRoot.addChild(tRoot);
                tRoot.setProperty(InternalProperties.ROOT, false);
            }
            superRoot.setProperty(InternalProperties.ROOT, true);
            superRoot.setProperty(InternalProperties.DUMMY, true);
            tGraph.getNodes().add(superRoot);
            break;
        }
        
        progressMonitor.done();
    }
}
