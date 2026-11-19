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
package org.eclipse.elk.alg.common.networksimplex;

/**
 * An edge in the graph processed by the {@link NetworkSimplex} algorithm. It has a source and
 * target {@link NNode}, a weight and a minimum length (delta).
 */
public class NEdge {

    // SUPPRESS CHECKSTYLE NEXT 22 VisibilityModifier
    /** A public id, unused internally, use it for whatever you want. */
    public int id;
    /** Internally set and used id to index arrays. */
    protected int internalId;
    /** An object from which this edge is derived. */
    public Object origin;
    /** The source node of this edge. */
    public NNode source;
    /** The target node of this edge. */
    public NNode target;
    /** The weight of this edge. */
    public double weight;
    /** The minimum length of this edge. */
    public int delta = 1;

    /**
     * A flag indicating whether a specified edge is part of the spanning tree determined by
     * {@code tightTree()}.
     * 
     * @see NetworkSimplex#tightTreeDFS(NNode)
     */
    protected boolean treeEdge = false;

    /**
     * @return an {@link NEdgeBuilder} to create a new edge.
     */
    public static NEdgeBuilder of() {
        return new NEdgeBuilder();
    }

    /**
     * @param origin
     *            an object from which this edge originates.
     * @return an {@link NEdgeBuilder} to create a new edge.
     */
    public static NEdgeBuilder of(final Object origin) {
        return new NEdgeBuilder();
    }
    
    /**
     * @return the source
     */
    public NNode getSource() {
        return source;
    }
    
    /**
     * @return the target
     */
    public NNode getTarget() {
        return target;
    }
    
    /**
     * @param some
     *            One of the source and target nodes of this edge.
     * @return the opposite node of {@code some}. That is, if {@code some} is the source node of
     *         this edge, the target node is returned and vice versa.
     * @throws IllegalArgumentException
     *             if some is neither target nor source of this edge.
     */
    public NNode getOther(final NNode some) {
        if (some == source) {
            return target;
        } else if (some == target) {
            return source;
        } else {
            throw new IllegalArgumentException("Node " + some + " not part of edge " + this);
        }
    }
    
    /**
     * Reversed this edge, i.e. the target becomes the source and the source becomes the target.
     * Also, the lists of incoming and outgoing edges of the source and target node are updated.
     * 
     * @return this.
     */
    public NEdge reverse() {
        NNode tmp = source;
        source = target;
        target = tmp;
        
        target.getOutgoingEdges().remove(this);
        target.getIncomingEdges().add(this);
        
        source.getIncomingEdges().remove(this);
        source.getOutgoingEdges().add(this);
        return this;
    }
    
    @Override
    public String toString() {
        return "NEdge[id=" + id + " w=" + weight + " d=" + delta + "]";
    }
    
    /**
     * Builder class for an {@link NEdge}. Allows to conveniently construct new edges.
     */
    public static final class NEdgeBuilder {
        
        private NEdge edge;
        
        private NEdgeBuilder() {
            edge = new NEdge();
        }

        /**
         * Sets the id field.
         * 
         * @param id
         *            an id
         * @return this builder.
         */
        public NEdgeBuilder id(final int id) {
            edge.id = id;
            return this;
        }

        /**
         * Sets the origin field.
         * 
         * @param origin
         *            any object.
         * @return this builder.
         */
        public NEdgeBuilder origin(final Object origin) {
            edge.origin = origin;
            return this;
        }

        /**
         * Sets the weight of this edge.
         * 
         * @param weight
         *            some positive weight.
         * @return this builder.
         */
        public NEdgeBuilder weight(final double weight) {
            edge.weight = weight;
            return this;
        }

        /**
         * Sets the minimal length of this edge.
         * 
         * @param delta
         *            a non-negative integer.
         * @return this builder.
         */
        public NEdgeBuilder delta(final int delta) {
            edge.delta = delta;
            return this;
        }
        
        /**
         * Sets the source node.
         * 
         * @param source
         *            a {@link NNode}.
         * @return this builder.
         */
        public NEdgeBuilder source(final NNode source) {
            edge.source = source;
            return this;
        }

        /**
         * Sets the target node.
         * 
         * @param target
         *            a {@link NNode}.
         * @return this builder.
         */
        public NEdgeBuilder target(final NNode target) {
            edge.target = target;
            return this;
        }

        /**
         * Finally returns the {@link NEdge} instance. As a side effect the edge is added to the
         * outgoing edges of the source {@link NNode} and to the incoming edges of the target
         * {@link NNode}.
         * 
         * @return the newly created {@link NEdge}.
         * 
         * @throws IllegalStateException
         *             if either no source or target was specified, or if source equals target.
         */
        public NEdge create() {
            
            if (edge.source == null || edge.target == null) {
                throw new IllegalStateException(NEdge.class.getSimpleName()
                        + " must have a source and target " + NNode.class.getSimpleName()
                        + " specified.");
            }
            
            if (edge.source == edge.target) {
                throw new IllegalStateException("Network simplex does not support self-loops: "
                        + edge + " " + edge.source + " " + edge.target);
            }
            
            edge.source.getOutgoingEdges().add(edge);
            edge.target.getIncomingEdges().add(edge);
            
            return edge;
        }
    }
}
