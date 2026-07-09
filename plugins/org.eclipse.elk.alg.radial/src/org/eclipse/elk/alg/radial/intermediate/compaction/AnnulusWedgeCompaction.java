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
package org.eclipse.elk.alg.radial.intermediate.compaction;

import java.util.Arrays;
import java.util.Collection;
import java.util.List;

import org.eclipse.elk.alg.radial.InternalProperties;
import org.eclipse.elk.alg.radial.RadialUtil;
import org.eclipse.elk.alg.radial.options.RadialOptions;
import org.eclipse.elk.alg.radial.sorting.IRadialSorter;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.graph.ElkNode;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;

/**
 * The AnnulusWedgeCompaction compacts each wedge one after another. Which leads to a more compact layout, but the
 * radial property will be lost.
 */
public class AnnulusWedgeCompaction extends AbstractRadiusExtensionCompaction implements IRadialCompactor {
    /** The left contour of each wedge, saved with the key of the first node in wedge. */
    private Multimap<ElkNode, ElkNode> leftContour = HashMultimap.create();

    /** The right contour of each wedge, saved with the key of the first node in wedge. */
    private Multimap<ElkNode, ElkNode> rightContour = HashMultimap.create();

    /** The sorter which shall be used and which is needed to determine the order of the nodes per radius. */
    private IRadialSorter sorter;

    /** The root node of the graph. */
    private ElkNode root;

    @Override
    public void compact(final ElkNode graph) {
        // Get values of the graph
        root = graph.getProperty(InternalProperties.ROOT_NODE);
        setRoot(root);
        this.sorter = graph.getProperty(RadialOptions.SORTER).create();
        Integer stepSize = graph.getProperty(RadialOptions.COMPACTION_STEP_SIZE);
        if (stepSize != null) {
            setCompactionStep(stepSize);
        }
        Double spacing = graph.getProperty(CoreOptions.SPACING_NODE_NODE);
        setSpacing(spacing);

        // Calculate the first level nodes
        List<ElkNode> successors = RadialUtil.getSuccessors(root);
        if (sorter != null) {
            sorter.sort(successors);
        }
        constructContour(successors);

        // contract each wedge
        List<ElkNode> rootList = Arrays.asList(new ElkNode[] { root });
        // do it two times to assure each node is compacted as much as possible
        for (int k = 0; k < 2; k++) {
            for (int i = 0; i < successors.size(); i++) {
                List<ElkNode> nodeAsList = Arrays.asList(new ElkNode[] { successors.get(i) });
                ElkNode rightParent = i < successors.size() - 1 ? successors.get(i + 1) : successors.get(0);
                ElkNode leftParent = i == 0 ? successors.get(successors.size() - 1) : successors.get(i - 1);

                contractWedge(successors.get(i), rootList, leftParent, rightParent, nodeAsList);
            }
        }

    }

    /**
     * Contract each wedge by shortening the incoming edge as long as no overlaps occurs.
     * 
     * @param wedgeParent
     *            The node which determines the wedge.
     * @param predecessors
     *            The tree predecessor.
     * @param radialPredecessor
     *            The radial predecessor wedge node.
     * @param radialSuccessor
     *            The radial successor wedge node.
     * @param currentRadiusNodes
     *            All nodes of the current wedge from one radius.
     */
    private void contractWedge(final ElkNode wedgeParent, final List<ElkNode> predecessors,
            final ElkNode radialPredecessor, final ElkNode radialSuccessor, final List<ElkNode> currentRadiusNodes) {
        boolean isOverlapping = overlapping(predecessors, radialPredecessor, radialSuccessor, currentRadiusNodes);
        boolean wasContracted = false;

        while (!isOverlapping) {
            contractLayer(currentRadiusNodes, true);
            wasContracted = true;
            isOverlapping = overlapping(predecessors, radialPredecessor, radialSuccessor, currentRadiusNodes);
        }

        // undo last step
        if (wasContracted) {
            contractLayer(currentRadiusNodes, false);
        }
        // continue with the nodes from the next radius
        List<ElkNode> nextLevelNodes = RadialUtil.getNextLevelNodes(currentRadiusNodes);
        if (!nextLevelNodes.isEmpty()) {
            if (sorter != null) {
                sorter.sort(nextLevelNodes);
            }
            contractWedge(wedgeParent, currentRadiusNodes, radialPredecessor, radialSuccessor, nextLevelNodes);
        }
    }

    /**
     * Calculate the overlaps. Overlaps are considered between the nodes of the wedge of one radius. Furthermore the
     * overlaps of the nodes and it's predecessors and between the leftmost node and the right contour of the wedge to
     * it's left are considered. And between the overlaps between the rightmost node and left contour of the wedge to
     * it's right.
     * 
     * @param layerNodes
     * @return
     */
    private boolean overlapping(final List<ElkNode> predecessors, final ElkNode leftParent, final ElkNode rightParent,
            final List<ElkNode> layerNodes) {
        if (sorter != null) {
            sorter.sort(layerNodes);
        }
        ElkNode firstNode = layerNodes.get(0);

        // overlap with left wedge contour
        if (contourOverlap(leftParent, firstNode, false)) {
            return true;
        }

        // overlap with right wegde contour
        ElkNode lastNode = layerNodes.get(layerNodes.size() - 1);
        if (contourOverlap(rightParent, lastNode, true)) {
            return true;
        }

        // overlaps on the radius
        if (overlapLayer(layerNodes)) {
            return true;
        }

        // overlaps with the predecessors
        for (ElkNode sortedNode : layerNodes) {
            for (ElkNode predecessor : predecessors) {
                if (overlap(sortedNode, predecessor)) {
                    return true;
                }
            }
        }

        return false;
    }

    /**
     * Check if a node overlaps with a neighboring wedge contour.
     * 
     * @param neighbourWedgeParent
     * @param node
     * @param left
     * @return
     */
    private boolean contourOverlap(final ElkNode neighbourWedgeParent, final ElkNode node, final boolean left) {
        Collection<ElkNode> contour =
                left ? leftContour.get(neighbourWedgeParent) : rightContour.get(neighbourWedgeParent);
        for (ElkNode contourNode : contour) {
            if (overlap(node, contourNode)) {
                return true;
            }
        }
        return false;
    }

    /**
     * Calculate the left and right contour of each node from the first layer.
     * 
     * @param nodes
     *            The list of first layer nodes.
     */
    private void constructContour(final List<ElkNode> nodes) {
        for (ElkNode node : nodes) {
            leftContour.put(node, node);
            rightContour.put(node, node);

            List<ElkNode> successors = RadialUtil.getSuccessors(node);
            if (!successors.isEmpty()) {
                if (sorter != null) {
                    sorter.sort(successors);
                }
                leftContour.put(node, successors.get(0));
                rightContour.put(node, successors.get(successors.size() - 1));

                while (!RadialUtil.getNextLevelNodes(successors).isEmpty()) {
                    successors = RadialUtil.getNextLevelNodes(successors);
                    if (sorter != null) {
                        sorter.sort(successors);
                    }
                    leftContour.put(node, successors.get(0));
                    rightContour.put(node, successors.get(successors.size() - 1));
                }

            }
        }
    }
}
