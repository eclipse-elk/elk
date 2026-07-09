/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
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
package org.eclipse.elk.alg.layered.intermediate.compaction;

import java.util.List;

import org.eclipse.elk.alg.common.compaction.oned.CGroup;
import org.eclipse.elk.alg.common.compaction.oned.CNode;
import org.eclipse.elk.alg.common.compaction.oned.OneDimensionalCompactor;
import org.eclipse.elk.alg.common.compaction.oned.ScanlineConstraintCalculator;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LNode.NodeType;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.layered.options.Spacings;
import org.eclipse.elk.core.UnsupportedConfigurationException;
import org.eclipse.elk.core.options.EdgeRouting;

import com.google.common.collect.Lists;

/**
 * An extended scanline procedure that is able to properly handle different spacing values between nodes and edges.
 */
public class EdgeAwareScanlineConstraintCalculation extends ScanlineConstraintCalculator {

    private static final double EPSILON = .5;
    private static final double SMALL_EPSILON = 0.01;
    
    private double verticalEdgeEdgeSpacing;
    
    private EdgeRouting edgeRouting;

    /**
     * 
     */
    public EdgeAwareScanlineConstraintCalculation(final LGraph graph) {
        this.verticalEdgeEdgeSpacing = graph.getProperty(LayeredOptions.SPACING_EDGE_EDGE);
        this.edgeRouting = graph.getProperty(LayeredOptions.EDGE_ROUTING);
    }
    
    @Override
    public void calculateConstraints(final OneDimensionalCompactor theCompactor) {
        this.compactor = theCompactor;

        switch (edgeRouting) {
        case ORTHOGONAL:
            calculateForOrthogonal();
            break;
        case SPLINES:
            calculateForSpline();
            break;
        default:
            throw new UnsupportedConfigurationException();
        }

    }
    
    // -------------------------------------- SPLINES --------------------------------------
    
    private void calculateForSpline() {
        final List<Runnable> schedule = Lists.newArrayList();
        
        /* -------------------- Vertical Segments -------------------- */
        // Note that some constraints between subsequent vertical segments of the same spline 
        //  have been precalculated during import 
        // Apart from that the boxes are not enlarged since this risks to introduce overlaps 
        
        sweep(n -> n.origin instanceof VerticalSegment);
        
        /* -------------------- Everything -------------------- */
        
        // find the minimum spacing
        double minSpacing = compactor.cGraph.cNodes.stream()
            .mapToDouble(n -> {
                // ignore external ports since their spacing is internally set to 0
                if (n.origin instanceof LNode && ((LNode) n.origin).getType() == NodeType.EXTERNAL_PORT) {
                    return Double.POSITIVE_INFINITY;
                }
                
                VerticalSegment vs = HorizontalGraphCompactor.getVerticalSegmentOrNull(n);
                if (vs != null) {
                    return Math.max(0, verticalEdgeEdgeSpacing / 2 - EPSILON);
                }
                LNode lNode = HorizontalGraphCompactor.getLNodeOrNull(n);
                if (lNode != null) {
                    double spacing = Spacings.getIndividualOrDefault(lNode, LayeredOptions.SPACING_NODE_NODE);
                    return Math.max(0, spacing / 2 - EPSILON);
                }
                
                return Double.POSITIVE_INFINITY;
            }).min().orElseGet(() -> 0.0);

        compactor.cGraph.cNodes.stream()
            .filter(n -> n.origin instanceof LNode)
            .forEach(n -> {
                alterHitbox(n, minSpacing, 1);
                schedule.add(() -> alterHitbox(n, minSpacing, -1));
            });

        sweep(n -> true);
        
        schedule.forEach(r -> r.run());
        schedule.clear();
        
    }
    
    // -------------------------------------- ORTHOGONAL --------------------------------------
    private void calculateForOrthogonal() {
        
        final List<Runnable> schedule = Lists.newArrayList();
        
        /* -------------------- Vertical Segments -------------------- */
        compactor.cGraph.cNodes.stream()
            .filter(n -> n.origin instanceof VerticalSegment)
            .forEach(n -> {
                double spacing = Math.max(0, verticalEdgeEdgeSpacing / 2 - EPSILON);
                alterHitbox(n, spacing, 1);
                
                schedule.add(() -> alterHitbox(n, spacing, -1));
            });
        
        sweep(n -> n.origin instanceof VerticalSegment);
        
        schedule.forEach(r -> r.run());
        schedule.clear();

        /* -------------------- Nodes -------------------- */
        compactor.cGraph.cNodes.stream()
            .filter(n -> n.origin instanceof LNode)
            .forEach(n -> {
                LNode lNode = HorizontalGraphCompactor.getLNodeOrNull(n);
                // One would want to use node-to-node spacing here but there's a problem:
                // if edge-node spacing is smaller than node-to-node/2, two nodes can be placed closer than they 
                // actually may by the previous phases ...
                //double spacing = Spacings.getIndividualOrDefault(lNode, LayeredOptions.SPACING_NODE_NODE);
                double spacing = Spacings.getIndividualOrDefault(lNode, LayeredOptions.SPACING_EDGE_EDGE);
                double finalSpacing = Math.max(0, spacing / 2 - EPSILON);
                alterHitbox(n, finalSpacing, 1);
             
                schedule.add(() -> alterHitbox(n, finalSpacing, -1));
            });
        
        sweep(n -> n.origin instanceof LNode);
        
        schedule.forEach(r -> r.run());
        schedule.clear();
        
        /* -------------------- Everything -------------------- */
        
        // find the minimum spacing
        double minSpacing = compactor.cGraph.cNodes.stream()
            .mapToDouble(n -> {
                // ignore external ports since their spacing is internally set to 0
                if (n.origin instanceof LNode && ((LNode) n.origin).getType() == NodeType.EXTERNAL_PORT) {
                    return Double.POSITIVE_INFINITY;
                }
                
                VerticalSegment vs = HorizontalGraphCompactor.getVerticalSegmentOrNull(n);
                if (vs != null) {
                    return Math.max(0, verticalEdgeEdgeSpacing / 2 - EPSILON);
                }
                LNode lNode = HorizontalGraphCompactor.getLNodeOrNull(n);
                if (lNode != null) {
                    double spacing = Spacings.getIndividualOrDefault(lNode, LayeredOptions.SPACING_NODE_NODE);
                    return Math.max(0, spacing / 2 - EPSILON);
                }
                
                return Double.POSITIVE_INFINITY;
            }).min().orElseGet(() -> 0.0);

        
        compactor.cGraph.cGroups.stream()
            .forEach(g -> {
                alterGroupedHitboxOrthogonal(g, minSpacing, 1);
                schedule.add(() -> alterGroupedHitboxOrthogonal(g, minSpacing, -1));
            });
        
        sweep(n -> true);
        
        schedule.forEach(r -> r.run());
        schedule.clear();
        
    }

    /**
     * Alters the hitboxes of the whole {@link CGroup}. That is, for orthogonal edge routes the {@link CGroup} is
     * constructed such that the group's {@link CGroup#master master} is a node and any further {@link CNode node} of
     * the group represents a {@link VerticalSegment}. When enlarging the vertical segments, care has to be taken to 
     * regard the {@link VerticalSegment#ignoreSpacing spacing ignores}.
     */
    private void alterGroupedHitboxOrthogonal(final CGroup g, final double spacing, final double fac) {
        CNode master = g.master;
        if (master == null) {
            master = g.cNodes.iterator().next();
        }

        // modify the master
        alterHitbox(master, spacing, fac);
        if (g.cNodes.size() == 1) {
            return;
        }

        // adjust all vertical segments
        double delta = spacing * fac;
        for (CNode n : g.cNodes) {
            if (n != master) {
                // n represents a vertical segment (by construction of the groups)
                VerticalSegment vs = HorizontalGraphCompactor.getVerticalSegmentOrNull(n);

                if (vs.ignoreSpacing.up) {
                    // it's a southward segment 
                    n.hitbox.y += delta + SMALL_EPSILON;
                    n.hitbox.height -= delta + SMALL_EPSILON;
                } else if (vs.ignoreSpacing.down) {
                    // northward
                    n.hitbox.height -= delta + SMALL_EPSILON;
                }
            }
        }
    }
    
    // -------------------------------------- COMMON --------------------------------------

    /**
     * Alter the hitbox of a {@link CNode} that either represents an {@link LNode} or a {@link VerticalSegment}.
     * For {@link VerticalSegment}s the {@link VerticalSegment#ignoreSpacing ignoreSpacing} is regarded.  
     */
    private void alterHitbox(final CNode n, final double spacing, final double fac) {
        double delta = spacing * fac;
        
        if (n.origin instanceof VerticalSegment) {
            VerticalSegment vs = HorizontalGraphCompactor.getVerticalSegmentOrNull(n);
            if (!vs.ignoreSpacing.up) {
                n.hitbox.y -= delta + SMALL_EPSILON;
                n.hitbox.height += delta + SMALL_EPSILON;
            } else if (!vs.ignoreSpacing.down) {
                n.hitbox.height += delta + SMALL_EPSILON;
            }
        } else if (n.origin instanceof LNode) {
            n.hitbox.y -= delta;
            n.hitbox.height += 2 * delta;
        }
    }
    
}
