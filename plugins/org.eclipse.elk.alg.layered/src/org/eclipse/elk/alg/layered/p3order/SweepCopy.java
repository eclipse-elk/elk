/*******************************************************************************
 * Copyright (c) 2016 Kiel University and others.
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
package org.eclipse.elk.alg.layered.p3order;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.HashSet;

import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.graph.LNode.NodeType;
import org.eclipse.elk.alg.layered.intermediate.PortListSorter;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.core.options.PortConstraints;
import org.eclipse.elk.core.options.PortSide;


/**
 * Stores node and port order for a sweep.
 * 
 * @author alan
 */
class SweepCopy {
    /** Saves a copy of the node order. */
    private final LNode[][] nodeOrder;
    /** Saves a copy of the orders of the ports on each node, because they are reordered in each sweep. */
    private final List<List<List<LPort>>> portOrders;

    /**
     * Copies on construction.
     * 
     * @param nodeOrderIn
     */
    SweepCopy(final LNode[][] nodeOrderIn) {
        nodeOrder = deepCopy(nodeOrderIn);
        portOrders = new ArrayList<>();
        for (LNode[] lNodes : nodeOrderIn) {
            List<List<LPort>> layer = new ArrayList<>();
            portOrders.add(layer);
            for (LNode node : lNodes) {
                layer.add(new ArrayList<>(node.getPorts()));
            }
        }
    }

    SweepCopy(final SweepCopy sc) {
        nodeOrder = deepCopy(sc.nodeOrder);
        portOrders = new ArrayList<>(sc.portOrders);
    }

    private LNode[][] deepCopy(final LNode[][] currentlyBestNodeOrder) {
        if (currentlyBestNodeOrder == null) {
            return null;
        }
        final LNode[][] result = new LNode[currentlyBestNodeOrder.length][];
        for (int i = 0; i < result.length; i++) {
            result[i] = Arrays.copyOf(currentlyBestNodeOrder[i], currentlyBestNodeOrder[i].length);
        }
        return result;
    }

    /**
     * Returns the copy of the node orders. WARNING: Do not change, or the copy will be invalid.
     * 
     * @return the nodeOrder
     */
    public LNode[][] nodes() {
        return nodeOrder;
    }

    /**
     * @param lGraph
     */
    public void transferNodeAndPortOrdersToGraph(final LGraph lGraph, final boolean setPortContstraints) {
        
        // the {@link LayeredOptions.ALLOW_NON_FLOW_PORTS_TO_SWITCH_SIDES} option allows the crossing minimizer 
        // to decide the side a corresponding dummy node is placed on in order to reduce the number of crossings
        // as a consequence the configured port side may not be valid anymore and has to be corrected
        List<LNode> northSouthPortDummies = new ArrayList<>();
        Set<LNode> updatePortOrder = new HashSet<>();
        
        // iterate the layers
        List<Layer> layers = lGraph.getLayers();
        for (int i = 0; i < layers.size(); i++) {
            List<LNode> nodes = layers.get(i).getNodes();
            northSouthPortDummies.clear();
            
            // iterate and order the nodes within the layer
            for (int j = 0; j < nodes.size(); j++) {
                LNode node = nodeOrder[i][j];
                // use the id field to remember the order within the layer
                node.id = j;
                if (node.getType() == NodeType.NORTH_SOUTH_PORT) {
                    northSouthPortDummies.add(node);
                }
                
                lGraph.getLayers().get(i).getNodes().set(j, node);
                // order ports as computed
                node.getPorts().clear();
                node.getPorts().addAll(portOrders.get(i).get(j));
                if (setPortContstraints) {
                    if (!node.getProperty(LayeredOptions.PORT_CONSTRAINTS).isOrderFixed()) {
                        node.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_ORDER);
                    }
                }
            }
            
            // assert that the port side is set properly
            for (LNode dummy : northSouthPortDummies) {
                LNode origin = assertCorrectPortSides(dummy);
                updatePortOrder.add(origin);
                updatePortOrder.add(dummy);
            }
        }

        // since the side of certain ports may have changed at this point, 
        // the list of ports must be re-sorted (see PortListSorter) 
        // and the port list views must be re-cached.
        for (LNode node : updatePortOrder) {
            Collections.sort(node.getPorts(), PortListSorter.CMP_COMBINED);
            node.cachePortSides();
        }
    }

    /**
     * Corrects the {@link PortSide} of dummy's origin.  
     * @return The {@link LNode} ('origin') whose port {@code dummy} represents. 
     */
    private LNode assertCorrectPortSides(final LNode dummy) {
        assert dummy.getType() == NodeType.NORTH_SOUTH_PORT;

        LNode origin = dummy.getProperty(InternalProperties.IN_LAYER_LAYOUT_UNIT);

        // a north south port dummy has exactly one port
        List<LPort> dummyPorts = dummy.getPorts();
        LPort dummyPort = dummyPorts.get(0);

        // find the corresponding port on the regular node
        for (LPort port : origin.getPorts()) {
            if (port.equals(dummyPort.getProperty(InternalProperties.ORIGIN))) {
                // switch the port's side if necessary
                if ((port.getSide() == PortSide.NORTH) && (dummy.id > origin.id)) {
                    port.setSide(PortSide.SOUTH);
                    if (port.isExplicitlySuppliedPortAnchor()) {
                        // Set new coordinates for port anchor since it was switched from NORTH to SOUTH.
                        // The y coordinate is updated by mirroring the y coordinate
                        double portHeight =  port.getSize().y;
                        double anchorY = port.getAnchor().y;
                        port.getAnchor().y = portHeight - anchorY;
                    }
                } else if ((port.getSide() == PortSide.SOUTH) && (origin.id > dummy.id)) {
                    port.setSide(PortSide.NORTH);
                    if (port.isExplicitlySuppliedPortAnchor()) {
                        // Set new coordinates for port anchor since it was switched from NORTH to SOUTH.
                        // The y coordinate is updated by mirroring the y coordinate
                        double portHeight =  port.getSize().y;
                        double anchorY = port.getAnchor().y;
                        port.getAnchor().y =  - (portHeight - anchorY);
                    }
                }
                break;
            }
        }
        return origin;
    }
    
}
