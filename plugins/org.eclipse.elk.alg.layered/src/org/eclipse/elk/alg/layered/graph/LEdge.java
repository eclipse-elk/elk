/*******************************************************************************
 * Copyright (c) 2010, 2019 Kiel University and others.
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
package org.eclipse.elk.alg.layered.graph;

import java.util.List;
import java.util.ArrayList;

import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.PortType;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.core.math.KVectorChain;
import org.eclipse.elk.core.options.EdgeLabelPlacement;
import org.eclipse.elk.core.options.PortSide;


/**
 * An edge in a layered graph. Edges may only be connected to ports of a node, which represent the point where the edge
 * touches the node.
 */
public final class LEdge extends LGraphElement {

    /** the serial version UID. */
    private static final long serialVersionUID = 1429497419118554817L;

    /** the bend points. */
    private KVectorChain bendPoints = new KVectorChain();
    /** the source port. */
    private LPort source;
    /** the target port. */
    private LPort target;
    /** labels assigned to this edge. */
    private final List<LLabel> labels = new ArrayList<>(3);

    /**
     * Reverses the edge, including its bend points. Negates the {@code REVERSED} property. (an edge that was marked as
     * being reversed is then unmarked, and the other way around) This does not change any properties on the connected
     * ports. End labels are reversed as well ( {@code HEAD} labels become {@code TAIL} labels and vice versa).
     * 
     * @param layeredGraph
     *            the layered graph
     * @param adaptPorts
     *            If true and a connected port is a collector port (a port used to merge edges), the corresponding
     *            opposite port is used instead of the original one.
     */
    public void reverse(final LGraph layeredGraph, final boolean adaptPorts) {
        LPort oldSource = getSource();
        LPort oldTarget = getTarget();

        setSource(null);
        setTarget(null);
        if (adaptPorts && oldTarget.getProperty(InternalProperties.INPUT_COLLECT)) {
            setSource(
                    LGraphUtil.provideCollectorPort(layeredGraph, oldTarget.getNode(), PortType.OUTPUT, PortSide.EAST));
        } else {
            setSource(oldTarget);
        }
        if (adaptPorts && oldSource.getProperty(InternalProperties.OUTPUT_COLLECT)) {
            setTarget(
                    LGraphUtil.provideCollectorPort(layeredGraph, oldSource.getNode(), PortType.INPUT, PortSide.WEST));
        } else {
            setTarget(oldSource);
        }

        // Switch end labels
        for (LLabel label : labels) {
            EdgeLabelPlacement labelPlacement = label.getProperty(LayeredOptions.EDGE_LABELS_PLACEMENT);

            if (labelPlacement == EdgeLabelPlacement.TAIL) {
                label.setProperty(LayeredOptions.EDGE_LABELS_PLACEMENT, EdgeLabelPlacement.HEAD);
            } else if (labelPlacement == EdgeLabelPlacement.HEAD) {
                label.setProperty(LayeredOptions.EDGE_LABELS_PLACEMENT, EdgeLabelPlacement.TAIL);
            }
        }

        boolean reversed = getProperty(InternalProperties.REVERSED);
        setProperty(InternalProperties.REVERSED, !reversed);

        bendPoints = KVectorChain.reverse(bendPoints);
    }

    /**
     * Returns the source port.
     * 
     * @return the source port
     */
    public LPort getSource() {
        return source;
    }

    /**
     * Sets the source port of this edge and adds itself to the port's list of edges. If the edge previously had another
     * source, it is removed from the original port's list of edges. Be careful not to use this method while iterating
     * through the edges list of the old port nor of the new port, since that could lead to
     * {@link java.util.ConcurrentModificationException}s.
     * 
     * @param source
     *            the source port to set
     */
    public void setSource(final LPort source) {
        if (this.source != null) {
            this.source.getOutgoingEdges().remove(this);
        }

        this.source = source;

        if (this.source != null) {
            this.source.getOutgoingEdges().add(this);
        }
    }

    /**
     * Returns the target port.
     * 
     * @return the target port
     */
    public LPort getTarget() {
        return target;
    }

    /**
     * Sets the target port of this edge and adds itself to the port's list of edges. If the edge previously had another
     * target, it is removed from the original port's list of edges. Be careful not to use this method while iterating
     * through the edges list of the old port nor of the new port, since that could lead to
     * {@link java.util.ConcurrentModificationException}s.
     * 
     * @param target
     *            the target port to set
     */
    public void setTarget(final LPort target) {
        if (this.target != null) {
            this.target.getIncomingEdges().remove(this);
        }

        this.target = target;

        if (this.target != null) {
            this.target.getIncomingEdges().add(this);
        }
    }

    /**
     * The same as {@link #setTarget(LPort)} with the exception that the index the edge is inserted at in the target
     * port's list of incoming edges can be specified.
     * 
     * <p>
     * <em>If you need to think about whether you want to use this or {@link #setTarget(LPort)}, chances are you want
     * the latter.</em>
     * </p>
     *
     * @param targetPort
     *            the target port to set
     * @param index
     *            the index to insert the edge at in the port's list of edges.
     * @throws IndexOutOfBoundsException
     *             if {@code target != null} and the index is invalid.
     */
    public void setTargetAndInsertAtIndex(final LPort targetPort, final int index) {
        if (this.target != null) {
            this.target.getIncomingEdges().remove(this);
        }

        this.target = targetPort;

        if (this.target != null) {
            // The insertion index below is the only difference to setTarget(LPort)
            this.target.getIncomingEdges().add(index, this);
        }
    }

    /**
     * Determines if this edge is a self-loop or not. An edge is considered a self-loop if both, source and target port
     * are defined and belong to the same non-null node.
     * 
     * @return {@code true} if this edge is a self-loop.
     */
    public boolean isSelfLoop() {
        if (this.source == null || this.target == null) {
            return false;
        }

        return this.source.getNode() != null && this.source.getNode() == this.target.getNode();
    }

    /**
     * @return {@code true} if this edge is not a self loop, the source node and target node of this edge reside in the
     *         same layer.
     * @see #isSelfLoop()
     */
    public boolean isInLayerEdge() {
        return !isSelfLoop() && (source.getNode().getLayer() == target.getNode().getLayer());
    }

    /**
     * Returns the list of bend points, with coordinates relative to the {@code LayeredGraph}'s origin. The list is
     * initially empty.
     * 
     * @return the bend points
     */
    public KVectorChain getBendPoints() {
        return bendPoints;
    }

    /**
     * Returns the list of edge labels.
     * 
     * @return all labels
     */
    public List<LLabel> getLabels() {
        return this.labels;
    }

    /**
     * @param port
     *            one of the ports of this edge.
     * @return the other port of this edge. That is, if {@code port} is the source port of this edge, the target port is
     *         returned and vice versa.
     * @throws IllegalArgumentException
     *             if {@code port} is neither target nor source of this edge.
     */
    public LPort getOther(final LPort port) {
        if (port == source) {
            return target;
        } else if (port == target) {
            return source;
        } else {
            throw new IllegalArgumentException("'port' must be either the source port or target port of the edge.");
        }
    }

    /**
     * @param node
     *            one of the nodes of this edge.
     * @return the other node of this edge. That is, if {@code node} is the source node of this edge, the target node is
     *         returned and vice versa.
     * @throws IllegalArgumentException
     *             if {@code node} is neither target nor source of this edge.
     */
    public LNode getOther(final LNode node) {
        if (node == source.getNode()) {
            return target.getNode();
        } else if (node == target.getNode()) {
            return source.getNode();
        } else {
            throw new IllegalArgumentException("'node' must either be the source node or target node of the edge.");
        }
    }

    @Override
    public String getDesignation() {
        if (!labels.isEmpty() && !(labels.get(0).getText() == null || labels.get(0).getText().isEmpty())) {
            return labels.get(0).getText();
        }
        return super.getDesignation();
    }

    @Override
    public String toString() {
        StringBuilder result = new StringBuilder();
        result.append("e_");
        String designation = getDesignation();
        if (designation != null) {
            result.append(designation);
        }
        if (source != null && target != null) {
            result.append(" ").append(source.getDesignation());
            result.append("[").append(source.getNode()).append("]");
            result.append(" -> ").append(target.getDesignation());
            result.append("[").append(target.getNode()).append("]");
        }
        return result.toString();
    }

}
