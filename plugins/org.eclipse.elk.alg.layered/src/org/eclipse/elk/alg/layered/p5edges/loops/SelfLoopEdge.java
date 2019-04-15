/*******************************************************************************
 * Copyright (c) 2018 Kiel University and others.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the Eclipse Public License v1.0
 * which accompanies this distribution, and is available at
 * http://www.eclipse.org/legal/epl-v10.html
 *******************************************************************************/
package org.eclipse.elk.alg.layered.p5edges.loops;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.core.options.PortSide;

/**
 * An edge which is part of a {@link SelfLoopComponent}. It holds the source and the target port, the original
 * {@link LEdge} and the component it belongs to.
 */
public class SelfLoopEdge {
    
    /** The self loop's source port. */
    private final SelfLoopPort source;
    /** The self loop's target port. */
    private final SelfLoopPort target;
    /** The self loop's original edge. */
    private final LEdge edge;
    /** The component the self loop belongs to. */
    private final SelfLoopComponent component;
    /** Other self loops this self loop depends on. */
    private final Map<PortSide, List<SelfLoopEdge>> dependencyEdges = new HashMap<>();
    /** TODO Document. */
    private final Map<PortSide, Integer> edgeOrders = new HashMap<>();

    /**
     * Create a new self loop.
     */
    public SelfLoopEdge(final SelfLoopComponent component, final SelfLoopPort source, final SelfLoopPort target,
            final LEdge edge) {
        this.component = component;
        this.source = source;
        this.target = target;
        this.edge = edge;
    }

    /**
     * Returns the self loop's source port.
     */
    public SelfLoopPort getSource() {
        return source;
    }

    /**
     * Returns the self loop's target port.
     */
    public SelfLoopPort getTarget() {
        return target;
    }

    /**
     * Returns the edge represented by this self loop.
     */
    public LEdge getEdge() {
        return edge;
    }

    /**
     * Returns the component this self loop belongs to.
     */
    public SelfLoopComponent getComponent() {
        return component;
    }

    /**
     * Returns the self loops this self loop depends on for each port side.
     */
    public Map<PortSide, List<SelfLoopEdge>> getDependencyEdges() {
        return dependencyEdges;
    }

    /**
     * TODO Document.
     */
    public Map<PortSide, Integer> getEdgeOrders() {
        return edgeOrders;
    }

    @Override
    public String toString() {
        return source + " -> " + target;
    }
}
