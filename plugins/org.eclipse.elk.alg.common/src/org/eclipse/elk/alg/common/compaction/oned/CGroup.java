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

import java.util.Set;

import com.google.common.collect.Sets;

/**
 * Represents a group of {@link CNode}s whose relative distances to each other are preserved.
 * For instance, when compacting a layered graph, {@link CGroup}s are used to ensure that vertical edge segments
 * that are connected to north/south ports, are kept at the position of the port (relative to the node).
 */
public final class CGroup {
    
    // Variables are public for convenience reasons
    // SUPPRESS CHECKSTYLE NEXT 34 VisibilityModifier

    /* -------------- Publicly Usable Fields -------------- */
    /** An id for public use. There is no warranty, use at your own risk. */
    public int id;
    /** Some node that was elected to be the master. The exact meaning of this is to be determined by the user. */
    public CNode master;

    /* -------------- Contained CNodes -------------- */
    /** grouped {@link CNode}s. */
    public Set<CNode> cNodes;
    
    /* -------------- Transient Data -------------- */
    /** root position of the {@link CGroup}. */
    public double startPos = Double.NEGATIVE_INFINITY;
    /** constraints pointing from within the {@link CGroup} to CNodes outside the {@link CGroup}. */
    public Set<CNode> incomingConstraints;
    /** number of constraints originating from within the {@link CGroup}. May change during compaction. */
    public int outDegree;
    /** number of constraints originating from within the {@link CGroup}. As opposed to {@link #outDegree} this value 
     * must not be altered by a compaction algorithm. */
    public int outDegreeReal;
    /** the reference node of this group, i.e. the reference for the group offset of other nodes.
     * The reference node is not fixed, it may be changed by the {@link OneDimensionalCompactor}. */
    protected CNode reference; 
    
    /** The field can be used to determine whether a group has moved during compaction. It has to be reset externally 
     * and is updated during compaction, for instance by the {@link LongestPathCompaction}. Bear in mind that not every
     * compaction algorithm updates this field. */
    public double delta = 0;
    /** Similar to {@link #delta} with the difference that it does not represent the direction-less
     * sum of a this group's movements but instead considers the compaction direction. Thus, if a
     * node moves back and forth the same distance this field's value will be zero. */
    public double deltaNormalized = 0;

    private CGroup() {
        cNodes = Sets.newLinkedHashSet();
        incomingConstraints = Sets.newHashSet();
        outDegree = 0;
        outDegreeReal = 0;
    }
    
    /**
     * @return a new builder for a {@link CGroup}.
     */
    public static CGroupBuilder of() {
        return new CGroupBuilder();
    }

    /**
     * Adds a {@link CNode} to the {@link CGroup} and updates the incoming constraints.
     * 
     * @param cNode
     *            the {@link CNode} to add
     */
    public void addCNode(final CNode cNode) {
        if (cNode.cGroup != null) {
            throw new RuntimeException("CNode belongs to another CGroup.");
        }
        cNodes.add(cNode);
        cNode.cGroup = this;
        if (reference == null) {
            reference = cNode;
        }
    }
    
    /**
     * Removes the passed {@link CNode} from this group and sets the {@link CNode#cGroup} field to
     * null (if it belongs to this group).
     * 
     * @param cNode
     *            the {@link CNode} to remove.
     * @return true if this {@link CGroup} actually contained the passed {@link CNode}, false
     *         otherwise.
     */
    public boolean removeCNode(final CNode cNode) {
        boolean removed = cNodes.remove(cNode);
        if (removed) {
            cNode.cGroup = null;
        }
        return removed;
    }
    
    /**
     * Builder class for {@link CGroup}s. 
     */
    public static final class CGroupBuilder {
        
        /** The group currently being constructed. */  
        private CGroup group;
        
        private CGroupBuilder() {
            group = new CGroup();
        }

        /**
         * Adds the passed nodes to this group.
         * 
         * @param nodes
         *            a non-negative integer.
         * @return this builder.
         */
        public CGroupBuilder nodes(final CNode... nodes) {
            for (CNode n : nodes) {
                group.addCNode(n);
            }
            return this;
        }
        
        /**
         * Sets the reference node of this group.
         * 
         * @param referenceNode
         * @return this builder.
         */
        public CGroupBuilder reference(final CNode referenceNode) {
            group.reference = referenceNode;
            return this;
        }

        /**
         * Sets the passed node as 'master'.
         * 
         * @param master
         * @return this builder.
         */
        public CGroupBuilder master(final CNode master) {
            group.master = master;
            return this;
        }
        

        /**
         * Finally creates this group. That is, the node is added to the passed {@link CGraph}.
         * 
         * @param graph
         *            the {@link CGraph} this node belongs to.
         * @return the created {@link CGroup}.
         */
        public CGroup create(final CGraph graph) {
            graph.cGroups.add(group);
            return group;
        }

    }
}
