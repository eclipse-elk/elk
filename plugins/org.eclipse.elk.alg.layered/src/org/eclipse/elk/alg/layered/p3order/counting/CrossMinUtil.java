/*******************************************************************************
 * Copyright (c) 2015 Kiel University and others.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the Eclipse Public License v1.0
 * which accompanies this distribution, and is available at
 * http://www.eclipse.org/legal/epl-v10.html
 *
 * Contributors:
 *     Kiel University - initial API and implementation
 *******************************************************************************/
package org.eclipse.elk.alg.layered.p3order.counting;

import java.util.Collections;

import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.core.options.PortSide;

import com.google.common.collect.Lists;

/**
 * Utility class for iterating over ports in any direction.
 * 
 * @author alan
 */
public final class CrossMinUtil {

    private CrossMinUtil() {
    };

    /**
     * Iterate over ports in north to south and east to west order.
     * 
     * @param node
     *            whose ports are being considered.
     * @param side
     *            of the node we are interested in
     * @return Iterable for ports on given node and side in given order.
     */
    public static Iterable<LPort> inNorthSouthEastWestOrder(final LNode node, final PortSide side) {
        switch (side) {
        case EAST:
        case NORTH:
            return node.getPorts(side);
        case SOUTH:
        case WEST:
            return Lists.reverse(node.getPorts(side));
        }
        return Collections.emptyList();
    }
    
}
