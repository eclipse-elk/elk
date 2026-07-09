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
package org.eclipse.elk.alg.radial.sorting;

import java.util.Comparator;
import java.util.List;

import org.eclipse.elk.alg.radial.RadialUtil;
import org.eclipse.elk.alg.radial.options.RadialOptions;
import org.eclipse.elk.graph.ElkNode;

/**
 * The polar coordinate sorter takes a list of nodes and sort them according to their polar coordinates relative to
 * their parents. Sublists of nodes with the same parent (which are not comparable by coordinates) stay in the same
 * order as in the original list.
 */
public class PolarCoordinateSorter implements IRadialSorter {
    private static final double DEGREE_45 = 0.25 * Math.PI;
    private static final double DEGREE_90 = 0.5 * Math.PI;
    private static final double DEGREE_135 = 0.75 * Math.PI;
    private static final double DEGREE_225 = 1.25 * Math.PI;
    private static final double DEGREE_270 = 1.5 * Math.PI;
    private static final double DEGREE_315 = 1.75 * Math.PI;

    private Comparator<ElkNode> compRight = RadialUtil.createPolarComparator(0, 0);
    private Comparator<ElkNode> compLeft = RadialUtil.createPolarComparator(Math.PI, 0);
    private Comparator<ElkNode> compTop = RadialUtil.createPolarComparator(DEGREE_270, 0);
    private Comparator<ElkNode> compBottom = RadialUtil.createPolarComparator(DEGREE_90, 0);
    private IDSorter idSorter;

    /**
     * Initialize the sorter by sorting them with the polar coordinater and then assigning
     * {@link RadialOptions.ORDER_ID}.
     * 
     * @param root
     *            The root node.
     */
    @Override
    public void initialize(final ElkNode root) {
        idSorter = new IDSorter();
        List<ElkNode> successors = RadialUtil.getSuccessors(root);
        // sort the first layer, the sorting starts at degree 0 the first polar coordinate position, which is on the
        // right of the circle.
        successors.sort(compRight);
        setIDForNodes(successors, 0);
    }

    /**
     * Iterates over the whole graph and assign an order id to each node, depending on its position in the list sorted
     * by polar coordinates.
     * 
     * @param nodes
     * @param idOffset
     * @return
     */
    private int setIDForNodes(final List<ElkNode> nodes, final int idOffset) {
        int id = idOffset;
        int nextLayerId = 0;
        for (ElkNode node : nodes) {
            node.setProperty(RadialOptions.ORDER_ID, id++);
            List<ElkNode> nodeSuccessors = RadialUtil.getSuccessors(node);
            double arc = Math.atan2(node.getY() + node.getHeight() / 2, node.getX() + node.getWidth() / 2);
            arc += arc < 0 ? 2 * Math.PI : 0;

            // node is right of parent node
            if (arc < DEGREE_45 || arc > DEGREE_315) {
                nodeSuccessors.sort(compLeft);
            } else if (arc <= DEGREE_315 && arc > DEGREE_225) { // node is below parent node
                nodeSuccessors.sort(compTop);
            } else if (arc <= DEGREE_225 && arc > DEGREE_135) {// node is left
                nodeSuccessors.sort(compRight);
            } else if (arc <= DEGREE_135) {// node is top
                nodeSuccessors.sort(compBottom);
            }

            nextLayerId = setIDForNodes(nodeSuccessors, nextLayerId);
        }
        return id;
    }

    @Override
    public void sort(final List<ElkNode> nodes) {
        if (!nodes.isEmpty()) {
            if (idSorter == null) {
                ElkNode root = RadialUtil.findRootOfNode(nodes.get(0));
                initialize(root);
            }
            idSorter.sort(nodes);
        }
    }

}
