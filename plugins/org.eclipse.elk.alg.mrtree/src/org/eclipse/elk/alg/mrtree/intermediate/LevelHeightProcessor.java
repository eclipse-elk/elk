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

import java.util.Collections;
import java.util.Iterator;

import org.eclipse.elk.alg.mrtree.graph.TGraph;
import org.eclipse.elk.alg.mrtree.graph.TNode;
import org.eclipse.elk.alg.mrtree.options.InternalProperties;
import org.eclipse.elk.alg.mrtree.options.MrTreeOptions;
import org.eclipse.elk.core.alg.ILayoutProcessor;
import org.eclipse.elk.core.options.Direction;
import org.eclipse.elk.core.util.IElkProgressMonitor;

import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;

/**
 * A processor which determines the height for each level by setting it to the height of the tallest
 * node of the level.
 * 
 * @author sor
 * @author sgu
 * @author sdo
 */
public class LevelHeightProcessor implements ILayoutProcessor<TGraph> {

    /** number of nodes in the graph. */
    private int numberOfNodes;

    /**
     * {@inheritDoc}
     */
    public void process(final TGraph tGraph, final IElkProgressMonitor progressMonitor) {
        progressMonitor.begin("Processor determine the height for each level", 1f);

        /** save number of nodes for progress computation */
        numberOfNodes = tGraph.getNodes().isEmpty() ? 1 : tGraph.getNodes().size();

        TNode root = null;
        /** find the root of the component */
        Iterator<TNode> it = tGraph.getNodes().iterator();
        while (root == null && it.hasNext()) {
            TNode tNode = it.next();
            if (tNode.getProperty(InternalProperties.ROOT)) {
                root = tNode;
            }
        }

        /** start with the root and level down by dsf */
        if (root != null) {
            setNeighbors(Lists.newArrayList(root), progressMonitor, tGraph.getProperty(MrTreeOptions.DIRECTION));
        }

        progressMonitor.done();
    }

    /**
     * Set the height property for each node in the current level and for their children. The height
     * is the height of the tallest node in the level.
     * 
     * @param currentLevel
     *            the list of TNode at the same level, for which the neighbors and siblings should
     *            be determined
     * @param progressMonitor
     *            the current progress monitor
     */
    private void setNeighbors(final Iterable<TNode> currentLevel, 
            final IElkProgressMonitor progressMonitor, final Direction layoutDirection) {
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

            double height = 0;
            if (layoutDirection.isHorizontal()) {
                for (TNode cN : currentLevel) {
                    /** append the children of the current node to the next level */
                    nextLevel = Iterables.concat(nextLevel, cN.getChildren());
                    /** check if the node is the tallest node so far */
                    if (height < cN.getSize().x) {
                        height = cN.getSize().x;
                    }
                }
            } else {
                for (TNode cN : currentLevel) {
                    /** append the children of the current node to the next level */
                    nextLevel = Iterables.concat(nextLevel, cN.getChildren());
                    /** check if the node is the tallest node so far */
                    if (height < cN.getSize().y) {
                        height = cN.getSize().y;
                    }
                }
            }
            for (TNode cN : currentLevel) {
                /** set the level height for the node */
                cN.setProperty(InternalProperties.LEVELHEIGHT, height);
            }

            /** add amount of work units to the whole task */
            sT.done();

            /** determine neighbors by bfs and for the whole graph */
            setNeighbors(nextLevel, progressMonitor, layoutDirection);
        }

    }
}
