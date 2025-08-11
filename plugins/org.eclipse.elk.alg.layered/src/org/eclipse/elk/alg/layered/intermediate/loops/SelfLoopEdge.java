/*******************************************************************************
 * Copyright (c) 2019 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.alg.layered.intermediate.loops;

import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LLabel;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.core.options.PortSide;

/**
 * Represents a single self loop edge.
 */
public class SelfLoopEdge {

    /** The edge represented by this instance. */
    private final LEdge lEdge;
    /** The self hyper loop this edge belongs to. */
    private SelfHyperLoop slHyperLoop;
    /** The edge's source port. */
    private final SelfLoopPort slSource;
    /** The edge's target port. */
    private final SelfLoopPort slTarget;

    /**
     * Creates a new instance that represents the given edge. Also adds the instance to the source and target port's
     * list of incoming and outgoing edges.
     */
    SelfLoopEdge(final LEdge lEdge, final SelfLoopPort slSource, final SelfLoopPort slTarget) {
        assert lEdge.isSelfLoop();
        assert slSource.getLPort() == lEdge.getSource();
        assert slTarget.getLPort() == lEdge.getTarget();

        this.lEdge = lEdge;
        this.slSource = slSource;
        this.slTarget = slTarget;
        
        slSource.getOutgoingSLEdges().add(this);
        slTarget.getIncomingSLEdges().add(this);
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Accessors
    
    /**
     * Returns the edge represented by this instance.
     */
    public LEdge getLEdge() {
        return lEdge;
    }
    
    /**
     * Returns the {@link SelfHyperLoop} this edge belongs to. This will only be {@code null} shortly after creation,
     * but gets set to a non-{@code null} value during self loop initialization.
     */
    public SelfHyperLoop getSLHyperLoop() {
        return slHyperLoop;
    }
    
    /**
     * Sets the {@link SelfHyperLoop} this edge should belong to.
     */
    void setSLHyperLoop(final SelfHyperLoop slLoop) {
        assert this.slHyperLoop == null;
        
        this.slHyperLoop = slLoop;
    }
    
    /**
     * Returns the edge's source port.
     */
    public SelfLoopPort getSLSource() {
        return slSource;
    }
    
    /**
     * Returns the edge's target port.
     */
    public SelfLoopPort getSLTarget() {
        return slTarget;
    }
    
    public boolean isInline() {
        for (LLabel label : lEdge.getLabels()) {
            if (label.getProperty(LayeredOptions.EDGE_LABELS_INLINE)) {
                return true;
            }
        }
        return false;
    }
    
    public PortSide getLabelSide() {
        return this.slHyperLoop.getSLLabels().getSide();
    }

}
