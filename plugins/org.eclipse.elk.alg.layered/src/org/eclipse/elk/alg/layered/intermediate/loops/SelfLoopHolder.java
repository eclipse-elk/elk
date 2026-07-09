/*******************************************************************************
 * Copyright (c) 2019 Kiel University and others.
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
package org.eclipse.elk.alg.layered.intermediate.loops;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Queue;
import java.util.stream.StreamSupport;

import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LNode.NodeType;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.core.options.PortSide;

/**
 * Holds all the information required to route self loops of a particular node. Each node with self loops has one of
 * these attached to it through the {@link InternalProperties#SELF_LOOP_HOLDER} property.
 */
public final class SelfLoopHolder {

    /** The node this instance belongs to. */
    private final LNode lNode;
    /** List of the node's {@link SelfHyperLoop}s. */
    private final List<SelfHyperLoop> slHyperLoops = new ArrayList<>();
    /**
     * Map of original ports to the {@link SelfLoopPort}s created to represent them. This is a {@link LinkedHashMap}
     * because we need a stable iteration order over the values later on so that we create {@link SelfHyperLoop}s in a
     * predictable order, thus preventing self loops from reordering between layout runs.
     */
    private final Map<LPort, SelfLoopPort> slPorts = new LinkedHashMap<>();
    
    /** Whether at least one self loop port is currently hidden from its node. */
    private boolean arePortsHidden = false;
    /** The number of routing slots on each side of the node, indexed by {@link PortSide} ordinal. */
    private int[] routingSlotCount = new int[PortSide.values().length];

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Creation

    private SelfLoopHolder(final LNode node) {
        this.lNode = node;
    }

    /**
     * Creates a new {@link SelfLoopHolder} for the given node and sets its {@link InternalProperties#SELF_LOOP_HOLDER}
     * property accordingly. This method should only be called if {@link #needsSelfLoopProcessing(LNode)} returns
     * {@code true}.
     * 
     * @param lNode
     *            the node to possibly install a {@link SelfLoopHolder} on.
     * @return the created {@link SelfLoopHolder}.
     */
    public static SelfLoopHolder install(final LNode lNode) {
        assert needsSelfLoopProcessing(lNode);
        
        SelfLoopHolder holder = new SelfLoopHolder(lNode);
        lNode.setProperty(InternalProperties.SELF_LOOP_HOLDER, holder);
        
        holder.initialize();

        return holder;
    }

    /**
     * Checks if the given node is a regular node and has at least one self loop.
     * 
     * @param lNode
     *            the node to check.
     * @return {@code true} if a {@link SelfLoopHolder} should be created for this node.
     */
    public static boolean needsSelfLoopProcessing(final LNode lNode) {
        if (lNode.getType() != NodeType.NORMAL) {
            return false;
        }

        return StreamSupport.stream(lNode.getOutgoingEdges().spliterator(), false).anyMatch(edge -> edge.isSelfLoop());
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Initialization
    
    /** ID of a port that hasn't been visited yet. */
    private static final int UNVISITED = 0;
    /** ID of a port that has already been visited. */
    private static final int VISITED = 1;

    /**
     * Initializes everything by populating our data model.
     */
    private void initialize() {
        List<SelfLoopEdge> slEdges = new ArrayList<>();

        // Create self loop edges and ports for every self loop
        for (LEdge lEdge : lNode.getOutgoingEdges()) {
            if (lEdge.isSelfLoop()) {
                slEdges.add(new SelfLoopEdge(lEdge, selfLoopPortFor(lEdge.getSource()),
                        selfLoopPortFor(lEdge.getTarget())));
            }
        }
        
        // Reset port IDs for the BFS we're about to run
        for (SelfLoopPort slPort : slPorts.values()) {
            slPort.getLPort().id = UNVISITED;
        }
        
        // Run BFS at every port to gather the edges into hyperloops
        for (SelfLoopPort slPort : slPorts.values()) {
            if (slPort.getLPort().id == UNVISITED) {
                slHyperLoops.add(initializeHyperLoop(slPort));
            }
        }
    }

    /**
     * Returns the {@link SelfLoopPort} representation of the given port. Creates one if none exists yet.
     */
    private SelfLoopPort selfLoopPortFor(final LPort lport) {
        if (!slPorts.containsKey(lport)) {
            slPorts.put(lport, new SelfLoopPort(lport));
        }

        return slPorts.get(lport);
    }
    
    /**
     * Collects all self loops reachable from the given port and merges them into a hyper loop.
     */
    private SelfHyperLoop initializeHyperLoop(final SelfLoopPort slPort) {
        SelfHyperLoop slLoop = new SelfHyperLoop(this);
        
        // Run a BFS starting at the port
        Queue<SelfLoopPort> bfsQueue = new LinkedList<>();
        bfsQueue.add(slPort);
        
        while (!bfsQueue.isEmpty()) {
            SelfLoopPort currentSLPort = bfsQueue.poll();
            currentSLPort.getLPort().id = VISITED;
            
            // Add each outgoing edge to our hyper loop
            for (SelfLoopEdge slEdge : currentSLPort.getOutgoingSLEdges()) {
                slLoop.addSelfLoopEdge(slEdge);
                
                SelfLoopPort slTargetPort = slEdge.getSLTarget();
                if (slTargetPort.getLPort().id == UNVISITED) {
                    bfsQueue.add(slTargetPort);
                }
            }
            
            // Add each incoming edge to our hyper loop
            for (SelfLoopEdge slEdge : currentSLPort.getIncomingSLEdges()) {
                slLoop.addSelfLoopEdge(slEdge);
                
                SelfLoopPort slSourcePort = slEdge.getSLSource();
                if (slSourcePort.getLPort().id == UNVISITED) {
                    bfsQueue.add(slSourcePort);
                }
            }
        }
        
        return slLoop;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Accessors
    
    /**
     * Returns the {@link LNode} represented by this instance.
     */
    public LNode getLNode() {
        return lNode;
    }
    
    /**
     * Returns the {@link SelfHyperLoop}s.
     */
    public List<SelfHyperLoop> getSLHyperLoops() {
        return slHyperLoops;
    }
    
    /**
     * Returns the map of {@link LPort}s to the {@link SelfLoopPort}s created for them.
     */
    public Map<LPort, SelfLoopPort> getSLPortMap() {
        return slPorts;
    }
    
    /**
     * Returns whether at least one self loop port is currently hidden from its node. Which one it is can be found by
     * iterating over our {@link SelfLoopPort}s and calling {@link SelfLoopPort#isHidden()}.
     */
    public boolean arePortsHidden() {
        return arePortsHidden;
    }
    
    /**
     * Sets whether at least one self loop port is currently hidden.
     */
    public void setPortsHidden(final boolean hidden) {
        arePortsHidden = hidden;
    }
    
    /**
     * Number of routing slots per side, indexed by {@link PortSide} ordinal. To be updated by clients.
     */
    public int[] getRoutingSlotCount() {
        return routingSlotCount;
    }

}
