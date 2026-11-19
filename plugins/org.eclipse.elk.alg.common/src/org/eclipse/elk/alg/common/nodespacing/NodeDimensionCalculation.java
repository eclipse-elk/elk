/*******************************************************************************
 * Copyright (c) 2014, 2015 Kiel University and others.
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
package org.eclipse.elk.alg.common.nodespacing;

import org.eclipse.elk.core.util.adapters.GraphAdapters.GraphAdapter;
import org.eclipse.elk.core.util.adapters.GraphAdapters.NodeAdapter;

/**
 * Entry points to apply several methods for node dimension calculation, including positioning of
 * labels, ports, etc.
 * 
 * @author uru
 */
public final class NodeDimensionCalculation {

    /**
     * Private constructor - utility class.
     */
    private NodeDimensionCalculation() {
    }

    /**
     * Calculates label sizes and node sizes also considering ports. Make sure that the port lists
     * are sorted properly.
     * 
     * @see LabelAndNodeSizeProcessor
     * @see #sortPortLists(GraphAdapter)
     * 
     * @param adapter
     *            an instance of an adapter for the passed graph's type.
     * @param <T>
     *            the graphs type, e.g. a root KNode
     */
    public static <T> void calculateLabelAndNodeSizes(final GraphAdapter<T> adapter) {
        NodeLabelAndSizeCalculator.process(adapter);
    }

    /**
     * Calculates node margins for the nodes of the passed graph.
     * 
     * <p>
     *   If certain functionality has to be excluded, e.g., the sizes of ports should be excluded from
     *   the calculated margins, use {@link #getNodeMarginCalculator(GraphAdapter)} to retrieve an
     *   instance of a {@link NodeMarginCalculator} that can be configured further.
     * </p>
     * 
     * @param adapter
     *            an instance of an adapter for the passed graph's type.
     * @param <T>
     *            the graphs type, e.g. a root KNode
     */
    public static <T> void calculateNodeMargins(final GraphAdapter<T> adapter) {
        NodeMarginCalculator calcu = new NodeMarginCalculator(adapter);
        calcu.process();
    }

    /**
     * <p>
     *   Returns a configurable {@link NodeMarginCalculator} that can be executed using the
     *   {@link NodeMarginCalculator#process()} method.
     * </p>
     * 
     * <p>
     *   Note that {@link #calculateNodeMargins(GraphAdapter)} can be used if no detailed
     *   configuration is required.
     * </p>
     * 
     * @param adapter
     *            an instance of an adapter for the passed graph's type.
     * @param <T>
     *            the graphs type, e.g. a root KNode
     * @return an instance of a {@link NodeMarginCalculator} that can be configured to specific
     *         needs.
     */
    public static <T> NodeMarginCalculator getNodeMarginCalculator(final GraphAdapter<T> adapter) {
        return new NodeMarginCalculator(adapter);
    }

    /**
     * Sorts the port lists of all nodes of the graph clockwise. More precisely, ports are sorted by
     * side (north, east, south, west) in clockwise order, beginning at the top left corner.
     * 
     * @param adapter
     *            an instance of an adapter for the passed graph's type.
     * @param <T>
     *            the graphs type, e.g. a root KNode
     */
    public static <T> void sortPortLists(final GraphAdapter<T> adapter) {
        // Iterate through the nodes of all layers
        for (NodeAdapter<?> node : adapter.getNodes()) {
            node.sortPortList();
        }
    }

}
