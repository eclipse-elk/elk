/*******************************************************************************
 * Copyright (c) 2013 - 2022 Kiel University and others.
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

import java.util.Collections;
import java.util.Iterator;

import org.eclipse.elk.alg.mrtree.graph.TGraph;
import org.eclipse.elk.alg.mrtree.graph.TNode;
import org.eclipse.elk.alg.mrtree.options.InternalProperties;
import org.eclipse.elk.core.alg.ILayoutProcessor;
import org.eclipse.elk.core.util.IElkProgressMonitor;

import com.google.common.collect.Iterables;

/**
 * A processor which determines the neighbors and siblings for all nodes in the graph. A neighbor is
 * the current node's nearest node, at the same level. A siblings is a neighbor with the same
 * parent.
 * 
 * @author sor
 * @author sgu
 * @author sdo
 */
public class NeighborsProcessor implements ILayoutProcessor<TGraph> {

    /** the number of nodes in the given graph. */
    private int numberOfNodes;

    /**
     * {@inheritDoc}
     */
    public void process(final TGraph tGraph, final IElkProgressMonitor progressMonitor) {
        progressMonitor.begin("Processor set neighbors", 1f);

        /** save number of nodes for progress computation */
        numberOfNodes = tGraph.getNodes().isEmpty() ? 1 : tGraph.getNodes().size();

        /** find the root of the component */
        TNode root = null;
        Iterator<TNode> it = tGraph.getNodes().iterator();
        while (root == null && it.hasNext()) {
            TNode tNode = it.next();
            if (tNode.getProperty(InternalProperties.ROOT)) {
                root = tNode;
            }
        }

        /** start with the root and level down by dsf */
        if (root != null) {            
            setNeighbors(root.getChildren(), progressMonitor);
        }

        progressMonitor.done();
    }

    /**
     * Set the neighbors of each node in the current level and for their children. A neighbor is the
     * current node's nearest node, at the same level. A sibling is a neighbor with the same
     * parent.
     * 
     * @param currentLevel
     *            the list of TNode at the same level, for which the neighbors and siblings should
     *            be determined
     * @param progressMonitor
     *            the current progress monitor
     */
    private void setNeighbors(final Iterable<TNode> currentLevel,
            final IElkProgressMonitor progressMonitor) {
        /** only do something in filled levels */
        if (!Iterables.isEmpty(currentLevel)) {
            /** create subtask for recursive descent */
            IElkProgressMonitor sT = progressMonitor.subTask(Iterables.size(currentLevel)
                    / numberOfNodes);

            sT.begin("Set neighbors in level", 1f);

            /** build empty iterator */
            Iterable<TNode> nextLevel = new Iterable<TNode>() {

                public Iterator<TNode> iterator() {
                    return Collections.emptyIterator();
                }
            };

            TNode lN = null;

            /**
             * the left neighbor is the previous processed node the right neighbor of the left
             * neighbor is the current node
             */
            for (TNode cN : currentLevel) {
                /** append the children of the current node to the next level */
                nextLevel = Iterables.concat(nextLevel, cN.getChildren());
                if (lN != null) {
                    lN.setProperty(InternalProperties.RIGHTNEIGHBOR, cN);
                    cN.setProperty(InternalProperties.LEFTNEIGHBOR, lN);
                    if (cN.getParent() == lN.getParent()) {
                        lN.setProperty(InternalProperties.RIGHTSIBLING, cN);
                        cN.setProperty(InternalProperties.LEFTSIBLING, lN);
                    }
                }

                lN = cN;
            }

            /** add amount of work units to the whole task */
            sT.done();

            /** determine neighbors by bfs and for the whole graph */
            setNeighbors(nextLevel, progressMonitor);
        }

    }
}
