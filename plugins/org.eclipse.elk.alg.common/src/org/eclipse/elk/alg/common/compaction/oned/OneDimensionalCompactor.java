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
package org.eclipse.elk.alg.common.compaction.oned;

import java.util.List;
import java.util.Map;

import org.eclipse.elk.core.math.ElkRectangle;
import org.eclipse.elk.core.options.Direction;
import org.eclipse.elk.core.util.Pair;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

/**
 * Implements the compaction of a {@link CGraph}. This includes the creation of a constraint graph
 * and its subsequent compaction.
 */
public final class OneDimensionalCompactor {
    
    /** Longest path-based compaction strategy. */
    public static final ICompactionAlgorithm LONGEST_PATH_COMPACTION = new LongestPathCompaction();
    /** Currently used compaction algorithm. */
    private ICompactionAlgorithm compactionAlgorithm = LONGEST_PATH_COMPACTION;

    /** Constraint calculation using a scanline technique. */
    public static final IConstraintCalculationAlgorithm SCANLINE_CONSTRAINTS =
            new ScanlineConstraintCalculator();
    /** Constraint calculation by pair-wise comparison of CNodes. */
    public static final IConstraintCalculationAlgorithm QUADRATIC_CONSTRAINTS =
            new QuadraticConstraintCalculation();
    /** Currently used instance of the constraint calculation algorithm. */
    private IConstraintCalculationAlgorithm constraintAlgorithm = SCANLINE_CONSTRAINTS;

    // SUPPRESS CHECKSTYLE NEXT 20 VisibilityModifier
    /** the {@link CGraph}. */
    public CGraph cGraph;
    /** A function evaluating for the passed direction if the node may be moved in that direction. */
    public ILockFunction lockFun;
    /** How spacings are determined. */
    public ISpacingsHandler spacingsHandler = ISpacingsHandler.DEFAULT_SPACING_HANDLER;
   
    /** compacting in this direction. */
    public Direction direction = Direction.UNDEFINED;
    /** flag indicating whether the {@link #finish()} method has been called. */
    private boolean finished = false;

    /**
     * Initializes the fields of the {@link OneDimensionalCompactor}.
     * 
     * @param cGraph
     *          the graph to compact
     */
    public OneDimensionalCompactor(final CGraph cGraph) {
        this.cGraph = cGraph;
        
        // for any pre-specified groups, deduce the offset of the elements
        calculateGroupOffsets();
        
        // the compaction operates solely on CGroups, 
        // thus, we wrap any plain CNodes into a CGroup
        for (CNode n : cGraph.cNodes) {
            if (n.cGroup == null) {
                CGroup.of()
                      .nodes(n)
                      .create(cGraph);
            }
            
            // remember the situation before compaction
            n.hitboxPreCompaction = new ElkRectangle(n.hitbox);
        }
    }
    
    /**
     * Sets the spacing handler to the passed handler. Not using this method results in the
     * {@link ISpacingsHandler#DEFAULT_SPACING_HANDLER} being used.
     * 
     * @param handler
     *            the {@link ISpacingsHandler} to use.
     * @return this instance of {@link OneDimensionalCompactor}
     */
    public OneDimensionalCompactor setSpacingsHandler(final ISpacingsHandler handler) {
        this.spacingsHandler = handler;
        return this;
    }

    /**
     * @param compactor
     *            the compaction algorithm to be used.
     * @return this instance of {@link OneDimensionalCompactor}
     */
    public OneDimensionalCompactor setCompactionAlgorithm(final ICompactionAlgorithm compactor) {
        this.compactionAlgorithm = compactor;
        return this;
    }
    
    /**
     * @param theConstraintAlgorithm
     *            the constraintAlgorithm to set.
     * @return this instance of {@link OneDimensionalCompactor}.
     */
    public OneDimensionalCompactor setConstraintAlgorithm(
            final IConstraintCalculationAlgorithm theConstraintAlgorithm) {
        this.constraintAlgorithm = theConstraintAlgorithm;
        return this;
    }
    
    /**
     * Compacts the graph in the specified direction.
     * 
     * @return
     *          this instance of {@link OneDimensionalCompactor}
     */
    public OneDimensionalCompactor compact() {
        
        if (finished) {
            throw new IllegalStateException("The " + getClass().getSimpleName()
                    + " instance has been finished already.");
        }
        
        // if no direction was specified the direction defaults to LEFT
        if (direction == Direction.UNDEFINED) {
            changeDirection(Direction.LEFT);
        }
        
        // prepare compaction by resetting initial outDegree value for groups
        for (CGroup g : cGraph.cGroups) {
            g.outDegree = g.outDegreeReal;
        }
        // ... and resetting nodes' positions
        for (CNode n : cGraph.cNodes) {
            n.startPos = Double.NEGATIVE_INFINITY;
        }
        
        // perform the actual compaction
        compactionAlgorithm.compact(this);

        return this;
    }
    
    /**
     * Call this method to indicate to the compacter that the compaction is finished now. No further
     * compaction steps are allowed. As a result, the direction is changed to the compactors natural
     * LEFT direction. In other words, if coordinates have been mirrored or transposed, they are
     * back to the original orientation now.
     * 
     * @return this instance of {@link OneDimensionalCompactor}. Note that the only use case after
     *         this call is to print debug output.
     */
    public OneDimensionalCompactor finish() {
        changeDirection(Direction.LEFT);
        finished = true;
        
        return this;
    }
    
    /**
     * Changes the direction for compaction by transforming the hitboxes.
     * 
     * @param dir
     *          the new direction
     * @return
     *          this instance of {@link OneDimensionalCompactor}
     */
    public OneDimensionalCompactor changeDirection(final Direction dir) {
        
        if (finished) {
            throw new IllegalStateException("The " + getClass().getSimpleName()
                    + " instance has been finished already.");
        }
        
        // prohibits vertical compaction if the graph has edges, because vertical and horizontal
        // segments would interchange
        if (!cGraph.supports(dir)) {
            throw new RuntimeException(
                    "The direction " + dir + " is not supported by the CGraph instance.");
        }

        // short-circuit if nothing todo
        if (dir == direction) {
            return this;
        }
        
        Direction oldDirection = direction;
        direction = dir;
        
        // executes required transformations and calculates constraints
        // if the graph is compacted in opposing directions the constraints are reversed instead of
        // being recalculated and CNodes are locked
        switch (oldDirection) {
        case UNDEFINED:
            switch (dir) {
            case LEFT:
                calculateConstraints();
                break;
                
            case RIGHT:
                mirrorHitboxes(); calculateConstraints();
                break;
                
            case UP:
                transposeHitboxes(); calculateConstraints();
                break;
                
            case DOWN:
                transposeHitboxes(); mirrorHitboxes(); calculateConstraints();
                break;

            default:
                break;
            }
            break;
            
        case LEFT:
            switch (dir) {
            case RIGHT:
                mirrorHitboxes(); reverseConstraints();
                break;
                
            case UP:
                transposeHitboxes(); calculateConstraints();
                break;
                
            case DOWN:
                transposeHitboxes(); mirrorHitboxes(); calculateConstraints();
                break;

            default:
                break;
            }
            break;
            
        case RIGHT:
            switch (dir) {
            case LEFT:
                mirrorHitboxes(); reverseConstraints();
                break;
                
            case UP:
                mirrorHitboxes(); transposeHitboxes(); calculateConstraints();
                break;
                
            case DOWN:
                mirrorHitboxes(); transposeHitboxes(); mirrorHitboxes(); calculateConstraints();
                break;

            default:
                break;
            }
            break;
            
        case UP:
            switch (dir) {
            case LEFT:
                transposeHitboxes(); calculateConstraints();
                break;
                
            case RIGHT:
                transposeHitboxes(); mirrorHitboxes(); calculateConstraints();
                break;
                
            case DOWN:
                mirrorHitboxes(); reverseConstraints();
                break;

            default:
                break;
            }
            break;
    
        case DOWN:
            switch (dir) {
            case LEFT:
                mirrorHitboxes(); transposeHitboxes(); calculateConstraints();
                break;
                
            case RIGHT:
                mirrorHitboxes(); transposeHitboxes(); mirrorHitboxes(); calculateConstraints();
                break;
                
            case UP:
                mirrorHitboxes(); reverseConstraints();
                break;

            default:
                break;
            }
            break;

        default:
            break;
        }
        
        return this;
    }
    
    /**
     * Iterates the {@link CNode}s of this group and returns {@code true} as soon as a node is locked in the passed
     * {@code direction}. Otherwise returns {@code false}.
     * 
     * @param group
     * @param dir
     * @return {@code true} if any {@link CNode} of this group is locked in the passed direction.
     */
    public boolean isLocked(final CGroup group, final Direction dir) {
        for (CNode n : group.cNodes) {
            if (isLocked(n, dir)) {
                return true;
            }
        }
        return false;
    }
    
    /**
     * @param node
     * @param dir
     * @return {@code true} if the node may <b>not</b> move in the passed direction. {@code false} otherwise.
     */
    public boolean isLocked(final CNode node, final Direction dir) {
        if (lockFun != null) {
            return lockFun.isLocked(node, dir);
        }
        return false;
    }
    
    /**
     * Sets a function specifying locks. May be {@code null} in which case all nodes 
     * are allowed to move arbitrarily.
     * 
     * @param fun
     * @return this instance of {@link OneDimensionalCompactor}.
     */
    public OneDimensionalCompactor setLockFunction(final ILockFunction fun) {
        lockFun = fun;
        return this;
    }
    
    /**
     * Initiates a recalculation of all constraints. Usually this method does not have to be called
     * manually since {@link #changeDirection(Direction)} either recalculates or reversed
     * constraints as necessary.
     * 
     * @return this instance of a {@link OneDimensionalCompactor}.
     */
    public OneDimensionalCompactor forceConstraintsRecalculation() {
        calculateConstraints();
        return this;
    }

    /**
     * For each {@link CGroup} g of the {@link CGroup} to be compacted, calculates the offsets of
     * g's nodes. It shouldn't be necessary to call this method in the common case since it's either
     * called by the constructor or for certain direction changes in
     * {@link #changeDirection(Direction)}.
     * 
     * The method assures that a group's {@link CGroup#reference} is always the left-most
     * element (i.e. the one with the lowest x-coordinate). As a consequence, all offsets
     * (in x) are positive. 
     * 
     * @return this instance of a {@link OneDimensionalCompactor}.
     */
    public OneDimensionalCompactor calculateGroupOffsets() {
        // for any pre-specified groups, deduce the offset of the elements
        for (CGroup group : cGraph.cGroups) {
            group.reference = null;

            // find left-most element
            for (CNode n : group.cNodes) {
                n.cGroupOffset.reset();
                if (group.reference == null || n.hitbox.x < group.reference.hitbox.x) {
                    group.reference = n;
                }
            }
            
            // calculate offsets
            for (CNode n : group.cNodes) {
                n.cGroupOffset.x = n.hitbox.x - group.reference.hitbox.x;
                n.cGroupOffset.y = n.hitbox.y - group.reference.hitbox.y;
            }
        }
        
        return this;
    }

    //////////////////////////////////////////// PRIVATE API ////////////////////////////////////////////
    
    /**
     * Mirrors hitboxes horizontally.
     */
    private void mirrorHitboxes() {
        for (CNode cNode : cGraph.cNodes) {
            cNode.hitbox.x = -cNode.hitbox.x - cNode.hitbox.width;
        }
        
        calculateGroupOffsets();
    }
    
    /**
     * Transposes hitboxes.
     */
    private void transposeHitboxes() {
        double tmp;
        for (CNode cNode : cGraph.cNodes) {
            tmp = cNode.hitbox.x;
            cNode.hitbox.x = cNode.hitbox.y;
            cNode.hitbox.y = tmp;
            
            tmp = cNode.hitbox.width;
            cNode.hitbox.width = cNode.hitbox.height;
            cNode.hitbox.height = tmp;

            tmp = cNode.cGroupOffset.x;
            cNode.cGroupOffset.x = cNode.cGroupOffset.y;
            cNode.cGroupOffset.y = tmp;
        }
        
        calculateGroupOffsets();
    }
    
    private void calculateConstraints() {
        
        // resetting constraints
        for (CNode cNode : cGraph.cNodes) {
            cNode.constraints.clear();
        }
        
        // apply any precalculated constraints
        List<Pair<CNode, CNode>> cstrs;
        if (direction.isHorizontal()) {
            cstrs = cGraph.predefinedHorizontalConstraints;
        } else {
            cstrs = cGraph.predefinedVerticalConstraints;
        }
        cstrs.forEach(p -> {
            if (direction == Direction.LEFT || direction == Direction.UP) {
                p.getFirst().constraints.add(p.getSecond());   
            } else {
                p.getSecond().constraints.add(p.getFirst());
            }
        });
        
        // run the specified constraint calculation algorithm
        constraintAlgorithm.calculateConstraints(this);
        
        // update the "external" constraints of the groups
        calculateConstraintsForCGroups();
    }
    
    /**
     * Sets the constraints of {@link CGroup}s according to the constraints of the included
     * {@link CNode}s.
     */
    private void calculateConstraintsForCGroups() {
        
        for (CGroup group : cGraph.cGroups) {
            group.outDegree = 0;
            group.outDegreeReal = 0;
            group.incomingConstraints.clear();
        }
        
        for (CGroup group : cGraph.cGroups) {
            for (CNode cNode : group.cNodes) {
                for (CNode inc : cNode.constraints) {

                    if (inc.cGroup != group) {
                        // adding constraints to each CGroup, that refer to CNodes outside of the
                        // CGroup
                        group.incomingConstraints.add(inc);
                        // setting the outDegree of CGroups
                        // (number of CNodes outside the group that the group is constrained by)
                        inc.cGroup.outDegree++;
                        inc.cGroup.outDegreeReal++;
                    }
                }
            }
        }
        
    }
    
    /**
     * If the graph is compacted twice in opposing directions, the constraints can be reversed to
     * avoid recalculating them. Also the {@link CNode}s' starting position is reset to be ready for
     * another compaction.
     */
    private void reverseConstraints() {

        // maps CNodes to temporary lists of incoming constraints
        Map<CNode, List<CNode>> incMap = Maps.newHashMap();
        for (CNode cNode : cGraph.cNodes) {
            incMap.put(cNode, Lists.newArrayList());
        }
        
        // resetting fields of CNodes and reversing constraints
        for (CNode cNode : cGraph.cNodes) {
            cNode.startPos = Double.NEGATIVE_INFINITY;
            for (CNode inc : cNode.constraints) {
                incMap.get(inc).add(cNode);
            }
        }
        
        // write back
        for (CNode cNode : cGraph.cNodes) {
            cNode.constraints.clear();
            cNode.constraints = incMap.get(cNode);
        }
        
        // resetting constraints for CGroups
        calculateConstraintsForCGroups();
    }
   
}
