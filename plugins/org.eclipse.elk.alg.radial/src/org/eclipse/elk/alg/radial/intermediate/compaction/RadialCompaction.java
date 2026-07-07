/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.radial.intermediate.compaction;

import java.util.List;

import org.eclipse.elk.alg.radial.InternalProperties;
import org.eclipse.elk.alg.radial.RadialUtil;
import org.eclipse.elk.alg.radial.options.RadialOptions;
import org.eclipse.elk.alg.radial.sorting.IRadialSorter;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.graph.ElkNode;

/**
 * A compaction strategy which compacts each radius one after another by edge shortening.
 */
public class RadialCompaction extends AbstractRadiusExtensionCompaction implements IRadialCompactor {

    /** The sorter which shall be used and which is needed to determine the order of the nodes per radius. */
    private IRadialSorter sorter;

    /** The value of the last radius. */
    private double lastRadius = 0;

    @Override
    public void compact(final ElkNode graph) {
        ElkNode root = graph.getProperty(InternalProperties.ROOT_NODE);
        setRoot(root);
        this.sorter = graph.getProperty(RadialOptions.SORTER).create();
        Integer stepSize = graph.getProperty(RadialOptions.COMPACTION_STEP_SIZE);
        if (stepSize != null) {
            setCompactionStep(stepSize);
        }
        Double spacing = graph.getProperty(CoreOptions.SPACING_NODE_NODE);
        setSpacing(spacing);

        List<ElkNode> firstLevelNodes = RadialUtil.getSuccessors(root);
        if (sorter != null) {
            sorter.sort(firstLevelNodes);
        }
        contract(firstLevelNodes);
    }

    /**
     * Contract each radius beginning at the inner radius. Each radius is contracted until an overlap occurs. The last
     * contraction will be undone, to remove the overlap.
     * 
     * @param nodes
     *            A list of nodes of one radius.
     */
    public void contract(final List<ElkNode> nodes) {
        if (!nodes.isEmpty()) {
            boolean isOverlapping = overlapping(nodes);
            boolean wasContracted = false;
            while (!isOverlapping) {
                contractLayer(nodes, true);
                wasContracted = true;
                isOverlapping = overlapping(nodes);
            }
            // undo last step
            if (wasContracted) {
                contractLayer(nodes, false);
            }
            List<ElkNode> nextLevelNodes = RadialUtil.getNextLevelNodes(nodes);
            if (sorter != null) {
                sorter.sort(nextLevelNodes);
            }
            lastRadius = calculateRadius(nodes.get(0));
            contract(nextLevelNodes);
        }
    }

    /**
     * Calculate the radius by help of a node.
     * 
     * @param node
     *            A node of the radius which will be determined.
     * @return the radius the given node is placed on.
     */
    private double calculateRadius(final ElkNode node) {
        double xPos = node.getX();
        double yPos = node.getY();

        ElkNode root = getRoot();
        double rootX = root.getX();
        double rootY = root.getY();

        double vectorX = xPos - rootX;
        double vectorY = yPos - rootY;
        double radius = Math.sqrt(vectorX * vectorX + vectorY * vectorY);
        return radius;
    }

    /**
     * Calculate if nodes of a radius are either overlapping their next neighbor or its predecessor on the next layer.
     * 
     * @param nodes A list of nodes of one radius.
     * @return Return if the value is overlapping or not.
     */
    private boolean overlapping(final List<ElkNode> nodes) {
        if (overlapLayer(nodes)) {
            return true;
        }
        for (ElkNode node : nodes) {
            ElkNode parent = RadialUtil.getTreeParent(node);
            if (overlap(node, parent)) {
                return true;

            }

            if (calculateRadius(node) - getSpacing() <= lastRadius) {
                return true;
            }
        }

        return false;
    }

}