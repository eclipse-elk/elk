/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.layered.options;

import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.p4nodes.NetworkSimplexPlacer;
import org.eclipse.elk.graph.properties.ExperimentalPropertyValue;

/**
 * Represents the flexibility of a node during node placement with the {@link NetworkSimplexPlacer}. Here, flexibility
 * means that either the ports are allowed to move on a node's border or that the node's size may change. Both points
 * have the goal to further reduce the overall edge length, possibly straightening a larger number of edges. 
 */
public enum NodeFlexibility {

    /** Node size must not be altered, ports remain where they were placed prior to node placement. */
    NONE,
    
    /** Ports are allowed to move on the node's border. */
    PORT_POSITION,
    
    /** Similar to {@link #NODE_SIZE} but tries to avoid enlarging nodes excessively. */ 
    @ExperimentalPropertyValue
    NODE_SIZE_WHERE_SPACE_PERMITS,
    
    /** The size of the node may be changed. This implies that also the ports may be re-positioned. */
    NODE_SIZE;
    
    /**
     * @return true if the flexibility is {@link #NODE_SIZE}.
     */
    public boolean isFlexibleSize() {
        return this == NODE_SIZE;
    }
    
    /**
     * @return true if the flexibility is <em>at least<em> {@link #NODE_SIZE_WHERE_SPACE_PERMITS}.
     */
    public boolean isFlexibleSizeWhereSpacePermits() {
        return this == NODE_SIZE_WHERE_SPACE_PERMITS || this == NODE_SIZE;
    }
    
    /**
     * @return @return true if the flexibility is <em>at least<em> {@link #PORT_POSITION}.
     */
    public boolean isFlexiblePorts() {
        return this == PORT_POSITION || this == NODE_SIZE_WHERE_SPACE_PERMITS || this == NODE_SIZE;
    }

    /**
     * @param nf
     *            the desired {@link NodeFlexibility}
     * @return true if the flexibility is <em>at least<em> {@code nf}.
     */
    public boolean isAtLeast(final NodeFlexibility nf) {
        switch (this) {
        case NODE_SIZE: return nf.isFlexibleSize();
        case NODE_SIZE_WHERE_SPACE_PERMITS: nf.isFlexibleSizeWhereSpacePermits();
        case PORT_POSITION: return nf.isFlexiblePorts();
        default: return true; // NONE
        }
    }
    
    /**
     * @param lNode
     *            the node for which the {@link NodeFlexibility} is requested.
     * @return the {@link LayeredOptions#NODE_PLACEMENT_NETWORK_SIMPLEX_NODE_FLEXIBILITY} that is specified for
     *         {@value lNode}, or if it is unspecified for {@value lNode}, the
     *         {@link LayeredOptions#NODE_PLACEMENT_NETWORK_SIMPLEX_NODE_FLEXIBILITY_DEFAULT} value of {@value lNode}'s
     *         parent.
     */
    public static NodeFlexibility getNodeFlexibility(final LNode lNode) {
        NodeFlexibility nf;
        if (lNode.getAllProperties().containsKey(LayeredOptions.NODE_PLACEMENT_NETWORK_SIMPLEX_NODE_FLEXIBILITY)) {
            nf = lNode.getProperty(LayeredOptions.NODE_PLACEMENT_NETWORK_SIMPLEX_NODE_FLEXIBILITY);
        } else {
            nf = lNode.getGraph().getProperty(LayeredOptions.NODE_PLACEMENT_NETWORK_SIMPLEX_NODE_FLEXIBILITY_DEFAULT);
        }
        return nf;
    }
}
