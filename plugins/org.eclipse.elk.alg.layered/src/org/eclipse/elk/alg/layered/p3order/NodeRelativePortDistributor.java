/*******************************************************************************
 * Copyright (c) 2010, 2015 Kiel University and others.
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
 * Calculates port ranks in a node-relative way.
 * 
 * <p>Port ranks are calculated by giving each node {@code i} a range of values {@code [i,i+1]} to
 * distribute their port values in. This effectively gives a node with many ports the same weight as
 * a node with few ports as opposed to a strategy that just numbers ports from top to bottom.</p>
 * 
 * @author cds
 * @author msp
 * @author ima
 */
public final class NodeRelativePortDistributor extends AbstractBarycenterPortDistributor {

    /**
     * Constructs a node-relative port distributor.
     * 
     * @param numLayers
     *            The number of layers in the graph.
     */
    public NodeRelativePortDistributor(final int numLayers) {
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
            float incr = 1.0f / (inputCount + 1);
            float northPos = rankSum + northInputCount * incr;
            float restPos = rankSum + 1 - incr;
            for (LPort port : node.getPorts(PortType.INPUT)) {
                if (port.getSide() == PortSide.NORTH) {
                    portRanks[port.id] = northPos;
                    northPos -= incr;
                } else {
                    portRanks[port.id] = restPos;
                    restPos -= incr;
                }
            }
            break;
        }

        case OUTPUT: {
            // Count the number of output ports
            int outputCount = 0;
            for (LPort port : node.getPorts()) {
                if (!port.getOutgoingEdges().isEmpty()) {
                    outputCount++;
                }
            }

            // Iterate output ports in their natural order, that is north - east - south - west
            float incr = 1.0f / (outputCount + 1);
            float pos = rankSum + incr;
            for (LPort port : node.getPorts(PortType.OUTPUT)) {
                portRanks[port.id] = pos;
                pos += incr;
            }
            break;
        }

        default:
            // this means illegal input to the method
            throw new IllegalArgumentException("Port type is undefined");
        }
        
        // the consumed rank is always 1
        return 1;
    }

}
