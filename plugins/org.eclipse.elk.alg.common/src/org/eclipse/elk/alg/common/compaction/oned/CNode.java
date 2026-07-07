/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.common.compaction.oned;

import java.util.List;
import java.util.function.Function;

import org.eclipse.elk.core.math.ElkRectangle;
import org.eclipse.elk.core.math.KVector;

import com.google.common.collect.Lists;

/**
 * Representation of a node/box in the constraint graph.
 * 
 * Use {@link CNode#of()} to construct a new instance. The method returns a builder object that is finalized using
 * {@link CNodeBuilder#create(CGraph)}.
 */
public final class CNode {
    
    // Variables are public for convenience reasons.
    // SUPPRESS CHECKSTYLE NEXT 42 VisibilityModifier
    
    /* -------------- Publicly Usable Fields -------------- */
    /** An optional id for public use (not used internally). There is no warranty, use at your own risk. */
    public int id;
    /** Element this node originates from; for public use (not used internally). */
    public Object origin;
    /** A more precise specification of what this node represents (not used internally). */
    public String type;
    /** A function called by {@link #toString()}. May be used for debugging. */
    public Function<CNode, String> toStringDelegate;
    
    /* -------------- Group Information -------------- */
    /** containing {@link CGroup}. */
    public CGroup cGroup;
    /** offset to the root position of the containing {@link CGroup} . */
    public KVector cGroupOffset = new KVector();

    /* -------------- Specification -------------- */
    /** the area occupied by this node prior to compaction. */
    public ElkRectangle hitboxPreCompaction;
    /** the area occupied by this node. */
    public ElkRectangle hitbox;

    /* -------------- Transient Data -------------- */
    /** Constraints of this node as determined by a {@link IConstraintCalculationAlgorithm}. */
    public List<CNode> constraints = Lists.newArrayList();
    /** Leftmost possible position for this {@link CNode} to be drawn. 
     *  This position may be intermediate as it is modified during the compaction process. */
    public double startPos = Double.NEGATIVE_INFINITY;
    
    
    private CNode() { }
    
    /**
     * @return a {@link CNodeBuilder} to construct a new node.
     */
    public static CNodeBuilder of() {
        return new CNodeBuilder();
    }
    
    /**
     * @return an svg representation of this {@link CNode} to the output for debugging.
     */
    public String getDebugSVG() {
        StringBuffer sb = new StringBuffer();
        sb.append("<rect x=\"" + this.hitbox.x + "\" y=\"" + this.hitbox.y + "\" width=\""
                + Math.max(1, this.hitbox.width) + "\" height=\"" + Math.max(1, this.hitbox.height)
                + "\" fill=\"" + "green" // (this.reposition ? "green" : "orange")
                + "\" stroke=\"black\" opacity=\"0.5\"/>");
        sb.append("<text x=\"" + (this.hitbox.x + 2) + "\" y=\""
                + (this.hitbox.y + 2 * 2 * 2 + 2 + 1) + "\">" + "(" + Math.round(this.hitbox.x)
                + ", " + Math.round(this.hitbox.y) + ")\n" + this + "</text>");
        return sb.toString();
    }
    
    @Override
    public String toString() {
        if (toStringDelegate != null) {
            return toStringDelegate.apply(this);
        }
        return super.toString();
    }
    
    /**
     * Builder class for {@link CNode}s. 
     */
    public static final class CNodeBuilder {
        
        /** The node currently being constructed. */  
        private CNode node;
        /** A node into which's {@link CGroup} this node should be placed upon creation. */
        private CNode groupParent;
        
        private CNodeBuilder() {
            node = new CNode();
        }

        /**
         * Sets the id field.
         * 
         * @param id
         *            a non-negative integer.
         * @return this builder.
         */
        public CNodeBuilder id(final int id) {
            node.id = id;
            return this;
        }

        /**
         * Sets the origin field.
         * 
         * @param origin
         *            any object.
         * @return this builder.
         */
        public CNodeBuilder origin(final Object origin) {
            node.origin = origin;
            return this;
        }

        /**
         * Sets the type field.
         * 
         * @param type
         *            any string, has no semantic meaning, can be used for debugging.
         * @return this builder.
         */
        public CNodeBuilder type(final String type) {
            node.type = type;
            return this;
        }
        
        /**
         * Sets the hitbox.
         * 
         * @param hitbox
         *            a rectangle defining the occupied space of this node.
         * @return this builder.
         */
        public CNodeBuilder hitbox(final ElkRectangle hitbox) {
            node.hitbox = hitbox;
            return this;
        }

        /**
         * Allows the customization of the {@link CNode}'s {@link #toString()} method.
         * 
         * @param delegate
         * 
         * @return this builder.
         */
        public CNodeBuilder toStringDelegate(final Function<CNode, String> delegate) {
            node.toStringDelegate = delegate;
            return this;
        }

        /**
         * Add this node to {@code parent}'s {@link CGroup}.
         * @param parent the node which should be the parent. May be {@code null}, in which case no grouping is applied.
         * @return this builder.
         */
        public CNodeBuilder groupWith(final CNode parent) {
            groupParent = parent;
            return this;
        }
        
        /**
         * Finally creates this node. That is, the node is added to the passed {@link CGraph}.
         * 
         * @param graph
         *            the {@link CGraph} this node belongs to.
         * @return the created {@link CNode}.
         */
        public CNode create(final CGraph graph) {
            graph.cNodes.add(node);
            if (groupParent != null) {
                if (groupParent.cGroup == null) {
                    CGroup.of()
                          .nodes(groupParent, node)
                          .create(graph);
                } else {
                    groupParent.cGroup.addCNode(node);
                }

            }
            return node;
        }

    }
}
