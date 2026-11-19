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
package org.eclipse.elk.alg.radial.intermediate.overlaps;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import org.eclipse.elk.alg.radial.InternalProperties;
import org.eclipse.elk.alg.radial.RadialUtil;
import org.eclipse.elk.alg.radial.intermediate.compaction.AbstractRadiusExtensionCompaction;
import org.eclipse.elk.alg.radial.options.RadialOptions;
import org.eclipse.elk.alg.radial.sorting.IRadialSorter;
import org.eclipse.elk.core.alg.ILayoutProcessor;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.util.IElkProgressMonitor;
import org.eclipse.elk.graph.ElkNode;

/**
 * Remove overlaps after the initial layout by extending layers with overlaps.
 */
public class RadiusExtensionOverlapRemoval extends AbstractRadiusExtensionCompaction
        implements ILayoutProcessor<ElkNode>, IOverlapRemoval {

    private IRadialSorter sorter;

    @Override
    public void removeOverlaps(final ElkNode graph, IElkProgressMonitor progressMonitor) {
        ElkNode root = graph.getProperty(InternalProperties.ROOT_NODE);
        setRoot(root);
        sorter = graph.getProperty(RadialOptions.SORTER).create();
        Double spacing = graph.getProperty(CoreOptions.SPACING_NODE_NODE);
        setSpacing(spacing);

        List<ElkNode> successors = RadialUtil.getSuccessors(root);
        extend(graph, successors, progressMonitor);
    }

    /**
     * Extend the first radius nodes until the nodes are non-overlapping.
     * 
     * @param nodes A list of the same radius.
     */
    public void extend(final ElkNode graph, final List<ElkNode> nodes, IElkProgressMonitor progressMonitor) {
        if (!nodes.isEmpty()) {
            List<KVector> oldPositions = new ArrayList<>();
            // Save old positions
            for (ElkNode node : nodes) {
                oldPositions.add(new KVector(node.getX(), node.getY()));
            }
            progressMonitor.logGraph(graph, "Before removing overlaps");
            while (overlapLayer(nodes)) {
                contractLayer(nodes, false);
            }
            progressMonitor.logGraph(graph, "After removing overlaps");

            double movedX = 0;
            double movedY = 0;
            ElkNode firstNode = null;
            if (!nodes.isEmpty()) {
                firstNode = nodes.get(0);
                movedX = firstNode.getX() - oldPositions.get(0).x;
                movedY = firstNode.getY() - oldPositions.get(0).y;
            }
            double movedDistance = Math.sqrt(movedX * movedX + movedY * movedY);
            Set<ElkNode> nextLevelNodes = RadialUtil.getNextLevelNodeSet(nodes);
            // Calculate the moved distance which is the amount all children and grandchildren have to be moved.
            int index = 1;
            if (!nextLevelNodes.isEmpty()) {
                for(ElkNode nextLevelNode : nextLevelNodes) {
                    moveNode(nextLevelNode, movedDistance);
                }
                progressMonitor.logGraph(graph, "Child movement " + index);
                index++;
            }            
            
            if (sorter != null) {
                sorter.sort(new ArrayList<>(nextLevelNodes));
            }
            extend(graph, new ArrayList<>(nextLevelNodes), progressMonitor);
        }
    }

    @Override
    public void process(final ElkNode graph, final IElkProgressMonitor progressMonitor) {
        progressMonitor.begin("Remove overlaps", 1);
        progressMonitor.logGraph(graph, "Before");
        removeOverlaps(graph, progressMonitor);
        progressMonitor.logGraph(graph, "After");
    }

}