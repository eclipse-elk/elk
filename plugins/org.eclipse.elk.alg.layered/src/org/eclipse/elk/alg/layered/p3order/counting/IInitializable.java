/*******************************************************************************
 * Copyright (c) 2016, 2020 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.layered.p3order.counting;

import java.util.List;

import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.p3order.LayerSweepCrossingMinimizer;

/**
 * Implement this interface to declare that data needs to be initialized before running the
 * {@link LayerSweepCrossingMinimizer}. We iterate through layers, nodes, ports and edges in the input order.
 */
public interface IInitializable {

    /**
     * Traverse graph and on each level call the respective methods for a list of initializable objects.
     * 
     * @param initializables
     *            The objects to be initialized
     */
    static void init(final List<IInitializable> initializables, final LNode[][] order) {
        for (int l = 0; l < order.length; l++) {
            for (IInitializable i : initializables) {
                i.initAtLayerLevel(l, order);
            }
            for (int n = 0; n < order[l].length; n++) {
                for (IInitializable i : initializables) {
                    i.initAtNodeLevel(l, n, order);
                }
                List<LPort> ports = order[l][n].getPorts();
                for (int p = 0; p < ports.size(); p++) {
                    for (IInitializable i : initializables) {
                        i.initAtPortLevel(l, n, p, order);
                    }
                    LPort port = ports.get(p);
                    int e = 0;
                    for (LEdge edge : port.getConnectedEdges()) {
                        for (IInitializable i : initializables) {
                            i.initAtEdgeLevel(l, n, p, e++, edge, order);
                        }
                    }
                }
            }
        }
        for (IInitializable i : initializables) {
            i.initAfterTraversal();
        }
    }

    /**
     * Initialize anything needed at layer level.
     * 
     * @param l
     *            layer index.
     * @param nodeOrder
     *            current node order for each layer.
     */
    default void initAtLayerLevel(final int l, final LNode[][] nodeOrder) {
    };

    /**
     * Initialize anything needed at node level.
     * 
     * @param l
     *            layer index.
     * @param n
     *            node index.
     * @param nodeOrder
     *            current node order for each layer.
     */
    default void initAtNodeLevel(final int l, final int n, final LNode[][] nodeOrder) {
    };

    /**
     * Initialize anything needed at port level.
     * 
     * @param l
     *            layer index.
     * @param n
     *            node index.
     * @param p
     *            port index.
     * @param nodeOrder
     *            current node order for each layer.
     */
    default void initAtPortLevel(final int l, final int n, final int p, final LNode[][] nodeOrder) {
    };

    /**
     * Initialize anything needed at edge level.
     * 
     * @param l
     *            layer index.
     * @param n
     *            node index.
     * @param p
     *            port index.
     * @param e
     *            edge index.
     * @param edge
     *            edge
     * @param nodeOrder
     *            current node order for each layer.
     */
    default void initAtEdgeLevel(final int l, final int n, final int p, final int e, final LEdge edge,
            final LNode[][] nodeOrder) {
    };

    /**
     * Finish initialization.
     */
    default void initAfterTraversal() {
    };
}
