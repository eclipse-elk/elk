/*******************************************************************************
 * Copyright (c) 2012, 2015 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.layered.p3order;

import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.options.PortType;
import org.eclipse.elk.core.options.PortSide;

/**
 * Calculates port ranks in a layer-total way.
 * 
 * <p>Port ranks are calculated by giving each node {@code i} a range of values {@code [i,i+x]},
 * where {@code x} is the number of input ports or output ports. Hence nodes with many ports are
 * assigned a broader range of ranks than nodes with few ports.</p>
 *
 * @author msp
 */
public final class LayerTotalPortDistributor extends AbstractBarycenterPortDistributor {
    
    /**
     * Constructs a layer-total port distributor.
     */
    public LayerTotalPortDistributor(final int numLayers) {
        super(numLayers);
    }

    @Override
    protected float calculatePortRanks(final LNode node, final float rankSum, final PortType type) {
        float[] portRanks = getPortRanks();

        switch (type) {
        case INPUT: {
            // Count the number of input ports, and additionally the north-side input ports
            int inputCount = 0, northInputCount = 0;
            for (LPort port : node.getPorts()) {
                if (!port.getIncomingEdges().isEmpty()) {
                    inputCount++;
                    if (port.getSide() == PortSide.NORTH) {
                        northInputCount++;
                    }
                }
            }

            // Assign port ranks in the order north - west - south - east
            float northPos = rankSum + northInputCount;
            float restPos = rankSum + inputCount;
            for (LPort port : node.getPorts(PortType.INPUT)) {
                if (port.getSide() == PortSide.NORTH) {
                    portRanks[port.id] = northPos;
                    northPos--;
                } else {
                    portRanks[port.id] = restPos;
                    restPos--;
                }
            }

            // the consumed rank corresponds to the number of input ports
            return inputCount;
        }

        case OUTPUT: {
            // Iterate output ports in their natural order, that is north - east - south - west
            int pos = 0;
            for (LPort port : node.getPorts(PortType.OUTPUT)) {
                pos++;
                portRanks[port.id] = rankSum + pos;
            }
            return pos;
        }

        default:
            // this means illegal input to the method
            throw new IllegalArgumentException();
        }

    }
}
