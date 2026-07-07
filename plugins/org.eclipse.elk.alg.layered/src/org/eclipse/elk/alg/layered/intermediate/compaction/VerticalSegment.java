/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.layered.intermediate.compaction;

import java.util.List;

import org.eclipse.elk.alg.common.compaction.oned.CNode;
import org.eclipse.elk.alg.common.compaction.oned.CompareFuzzy;
import org.eclipse.elk.alg.common.compaction.oned.Quadruplet;
import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.core.math.ElkRectangle;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.math.KVectorChain;

import com.google.common.base.Joiner;
import com.google.common.collect.Lists;
import com.google.common.math.DoubleMath;

/**
 * Represents a vertical segment on a single {@link LEdge} that is merged with intersecting
 * {@link VerticalSegment}s.
 * <p>
 * There are two types of vertical segments that are handled differently and are distinguished by
 * the presence of a parent node. The usual segment stems from a step in an edge and has no parent
 * node. The other type is the vertical part of an edge that is connected to the north or south side
 * of a node. In this case the node is the parent of the segment and neither impose any constraints on
 * each other because they are supposed to connect and therefore collide.
 * </p>
 */
public final class VerticalSegment implements Comparable<VerticalSegment> {
    // SUPPRESS CHECKSTYLE NEXT 34 VisibilityModifier
    
    /** Nodes that may become the parent of the {@link CNode} that will represent this segment. */
    public List<CNode> potentialGroupParents = Lists.newArrayList();
    
    /** The edges that contribute at least partly to this vertical segment. */
    public List<LEdge> representedLEdges = Lists.newArrayList();
    /** The bendpoints that are located within the hitbox of this vertical segment.
     *  After compaction they should be adjusted. */
    public List<KVector> affectedBends = Lists.newArrayList();
    /** Bounding boxes, e.g. of splines, that are to be adjusted after compaction. */
    public List<ElkRectangle> affectedBoundingBoxes = Lists.newArrayList();

    /** The area occupied by this vertical segment. */
    public ElkRectangle hitbox = new ElkRectangle();
    
    /** The junction points of the original {@link LEdge}, that are between the bend points. */
    public KVectorChain junctionPoints = new KVectorChain();
    
    /** Whether spacing on a particular side should be ignored. */
    public Quadruplet ignoreSpacing = new Quadruplet();
    
    /** Pre-computed constraints that are to be added to the {@link CNode} which represents this segment. */
    public List<VerticalSegment> constraints = Lists.newArrayList();

    /** Orthogonal edges may have vertical segments that leave/enter a north/south port. Such a port is stored here. */
    public LPort aPort;
    
    /** Segments that have been joined with this one. */
    public List<VerticalSegment> joined = Lists.newArrayList();
    
    /**
     * Creates a new instance setting all fields.
     * 
     * @param bend1
     *            first bend point
     * @param bend2
     *            second bend point
     * @param cNode
     *            parent node
     * @param lEdge
     *            parent edge
     */
    public VerticalSegment(final KVector bend1, final KVector bend2, final CNode cNode,
            final LEdge lEdge) {

        affectedBends.add(bend1);
        affectedBends.add(bend2);

        hitbox.x = Math.min(bend1.x, bend2.x);
        hitbox.y = Math.min(bend1.y, bend2.y);
        hitbox.width = Math.abs(bend1.x - bend2.x);
        hitbox.height = Math.abs(bend1.y - bend2.y);
        
        // adding junction points
        KVectorChain inJPs = lEdge.getProperty(LayeredOptions.JUNCTION_POINTS);
        // this property could be null
        if (inJPs != null) {
            for (KVector jp : inJPs) {
                if (CompareFuzzy.eq(jp.x, bend1.x)) {
                    junctionPoints.add(jp);
                }
            }
        }

        if (cNode != null) {
            potentialGroupParents.add(cNode);
        }
//        parentNode = cNode;

        representedLEdges.add(lEdge);
    }
    
    /**
     * Joins this segment with the {@code other} segment. As such, this segment is modified to capture the union of the
     * two segments. The other segment remains unaltered.
     * 
     * @param other
     *            another vertical segment.
     * @return this instance.
     */
    public VerticalSegment joinWith(final VerticalSegment other) {
       
        representedLEdges.addAll(other.representedLEdges);
        affectedBends.addAll(other.affectedBends);
        affectedBoundingBoxes.addAll(other.affectedBoundingBoxes);
        junctionPoints.addAll(other.junctionPoints);
        constraints.addAll(other.constraints);
        potentialGroupParents.addAll(other.potentialGroupParents);
        
        double newX = Math.min(hitbox.x, other.hitbox.x);
        double newY = Math.min(hitbox.y, other.hitbox.y);
        double maxX = Math.max(hitbox.x + hitbox.width, other.hitbox.x + other.hitbox.width);
        double newW = maxX - newX;
        double maxY = Math.max(hitbox.y + hitbox.height, other.hitbox.y + other.hitbox.height);
        double newH = maxY - newY;
        hitbox.setRect(newX, newY, newW, newH);

        ignoreSpacing.applyOr(other.ignoreSpacing);
        
        if (aPort == null) {
            aPort = other.aPort;
        }
        
        joined.addAll(other.joined);
        joined.add(other);
        
        return this;
    }

    /**
     * Checks whether this segment intersects another.
     * 
     * @param o
     *            the other segment
     * @return {@code true} if the segments intersect
     */
    public boolean intersects(final VerticalSegment o) {
        return CompareFuzzy.eq(this.hitbox.x, o.hitbox.x)
                && !(CompareFuzzy.lt(this.hitbox.getBottomLeft().y, o.hitbox.y) 
                        || CompareFuzzy.lt(o.hitbox.getBottomLeft().y, this.hitbox.y));
    }

    @Override
    public int compareTo(final VerticalSegment o) {
        int d = DoubleMath.fuzzyCompare(hitbox.x, o.hitbox.x, CompareFuzzy.TOLERANCE);
        if (d == 0) {
            return Double.compare(this.hitbox.y, o.hitbox.y);
        }
        return d;
    }
    
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append("VerticalSegment ");
        sb.append(hitbox);
        sb.append(" ");
        sb.append(Joiner.on(", ").join(representedLEdges));
        return sb.toString();
    }
}
