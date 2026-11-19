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
package org.eclipse.elk.alg.layered.compaction.oned;

import java.util.Set;

import com.google.common.collect.Sets;

/**
 * Represents a group of {@link CNode}s whose relative distances to each other are preserved.
 * For instance, when compacting a layered graph, CGroups are used to ensure that vertical edge segments,
 * that are connected to north/south ports, are kept at the position of the port.
 * 
 * @see org.eclipse.elk.alg.layered.intermediate.compaction.HorizontalGraphCompactor
 */
public final class CGroup {
    
    // Variables are public for convenience reasons since this class is used internally only.
    // SUPPRESS CHECKSTYLE NEXT 30 VisibilityModifier
    /** root position of the {@link CGroup}. */
    public double startPos = Double.NEGATIVE_INFINITY;
    /**
     * The field can be used to determine whether a group has moved during compaction. It has to be
     * reset externally and is updated during compaction, for instance by the
     * {@link org.eclipse.elk.alg.layered.compaction.oned.algs.LongestPathCompaction
     * LongestPathCompaction}. Bear in mind that not every compaction algorithm updates this field.
     */
    public double delta = 0;
    /**
     * Similar to {@link #delta} with the difference that it does not represent the direction-less
     * sum of a this group's movements but instead considers the compaction direction. Thus, if a
     * node moves back and forth the same distance this field's value will be zero.
     */
    public double deltaNormalized = 0;
    /** grouped {@link CNode}s. */
    public Set<CNode> cNodes;
    /** constraints pointing from within the {@link CGroup} to CNodes outside the {@link CGroup}. */
    public Set<CNode> incomingConstraints;
    /** number of constraints originating from within the {@link CGroup}. */
    public int outDegree;
    /** the reference node of this group, i.e. the reference for the group offset of other nodes. */
    public CNode reference; 
    /** An id for public use. There is no warranty, use at your own risk. */
    public int id;
    /** flags this group to be repositioned in the case of left/right balanced compaction. */
    public boolean reposition = true;

    /**
     * The constructor for a {@link CGroup} receives {@link CNode}s to group.
     * 
     * @param inputCNodes
     *            the {@link CNode}s to add
     */
    public CGroup(final CNode... inputCNodes) {
        cNodes = Sets.newLinkedHashSet();
        incomingConstraints = Sets.newHashSet();
        outDegree = 0;
        for (CNode cNode : inputCNodes) {
            if (reference == null) {
                reference = cNode;
            }
            addCNode(cNode);
        }
    }

    /**
     * Adds a {@link CNode} to the {@link CGroup} and updates the incoming constraints.
     * 
     * @param cNode
     *            the {@link CNode} to add
     */
    public void addCNode(final CNode cNode) {
        cNodes.add(cNode);
        if (cNode.cGroup != null) {
            throw new RuntimeException("CNode belongs to another CGroup.");
        }
        cNode.cGroup = this;
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
}
