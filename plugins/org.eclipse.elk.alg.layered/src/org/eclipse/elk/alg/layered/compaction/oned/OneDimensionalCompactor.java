/*******************************************************************************
 * Copyright (c) 2016 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.alg.layered.compaction.oned;

import java.util.List;
import java.util.Map;

import org.eclipse.elk.alg.layered.compaction.oned.algs.ICompactionAlgorithm;
import org.eclipse.elk.alg.layered.compaction.oned.algs.IConstraintCalculationAlgorithm;
import org.eclipse.elk.alg.layered.compaction.oned.algs.LongestPathCompaction;
import org.eclipse.elk.alg.layered.compaction.oned.algs.QuadraticConstraintCalculation;
import org.eclipse.elk.alg.layered.compaction.oned.algs.ScanlineConstraintCalculator;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.options.Direction;
import org.eclipse.elk.core.util.IElkProgressMonitor;
import org.eclipse.elk.core.util.Pair;
import org.eclipse.elk.core.util.LoggedGraph.Type;

import com.google.common.base.Function;
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
    /** compacting in this direction. */
    public Direction direction = Direction.UNDEFINED;
    /** a function that sets the {@link CNode#reposition} flag according to the direction. */
    private Function<Pair<CNode, Direction>, Boolean> lockingStrategy;
    /** flag indicating whether the {@link #finish()} method has been called. */
    private boolean finished = false;
    /** How spacings are determined. */
    public ISpacingsHandler<? super CNode> spacingsHandler = ISpacingsHandler.DEFAULT_SPACING_HANDLER;
    
    /**
     * Initializes the fields of the {@link OneDimensionalCompactor}.
     * 
     * @param cGraph
     *          the graph to compact
     */
    public OneDimensionalCompactor(final CGraph cGraph) {
        
        this.cGraph = cGraph;
        // the default locking strategy locks CNodes if they are not constrained
        setLockingStrategy((pair) -> !(pair.getFirst().cGroup.outDegree == 0));
        
        // for any pre-specified groups, deduce the offset of the elements
        calculateGroupOffsets();
        
        // the compaction operates solely on CGroups, 
        // thus, we wrap any plain CNodes into a CGroup
        for (CNode n : cGraph.cNodes) {
            if (n.cGroup == null) {
                CGroup group = new CGroup(n);
                cGraph.cGroups.add(group);
            }
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
    public OneDimensionalCompactor setSpacingsHandler(final ISpacingsHandler<? super CNode> handler) {
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
        
        // prepare compaction by setting initial outDegree value for groups
        // iterates once over all constraints
        for (CGroup g : cGraph.cGroups) {
            g.outDegree = 0;
        }

        for (CNode n : cGraph.cNodes) {
            // reset any previously calculated start positions
            n.startPos = Double.NEGATIVE_INFINITY;
            
            for (CNode incN : n.constraints) {
                incN.cGroup.outDegree++;
            }
        }
        
        // perform the actual compaction
        compactionAlgorithm.compact(this);
        
        // remove all locks
        for (CNode node : cGraph.cNodes) {
            node.reposition = true;
        }
        
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
     * Sets the function to lock cGraph.cNodes. Does not apply the strategy yet, to do so, call
     * {@link #applyLockingStrategy(Direction)}.
     * 
     * @param strategy
     *            a function that returns the required state for the reposition flag of a
     *            {@link CNode} according to the direction
     * @return this instance of {@link OneDimensionalCompactor}
     */
    public OneDimensionalCompactor setLockingStrategy(
            final Function<Pair<CNode, Direction>, Boolean> strategy) {
        lockingStrategy = strategy;
        return this;
    }
    
    /**
     * Locks nodes according to the currently set locking strategy and the current direction. After
     * the next call to {@link #compact()}, all locks are removed again.
     * 
     * @return this instance of {@link OneDimensionalCompactor}
     * 
     * @see #setLockingStrategy(Function)
     */
    public OneDimensionalCompactor applyLockingStrategy() {
        applyLockingStrategy(direction);
        return this;
    }

    /**
     * Locks nodes according to the currently set locking strategy and the current direction. After
     * the next call to {@link #compact()}, all locks are removed again.
     * 
     * @param dir
     *            the {@link Direction} for the locks to be set.
     * @return this instance of {@link OneDimensionalCompactor}
     * 
     * @see #setLockingStrategy(Function)
     */
    public OneDimensionalCompactor applyLockingStrategy(final Direction dir) {

        for (CGroup cGroup : cGraph.cGroups) {
            cGroup.reposition = true;
        }
        
        for (CNode cNode : cGraph.cNodes) {
            cNode.reposition = lockingStrategy.apply(Pair.of(cNode, dir));
            
            cNode.cGroup.reposition &= lockingStrategy.apply(Pair.of(cNode, dir));
        }
        
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
     * The method assures that a groups {@link CGroup#reference} is always the left-most
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
            // mirroring the offsets inside CGroups
            // FIXME do we still require this?!
            if (cNode.parentNode != null) {
                cNode.cGroupOffset.x = -cNode.cGroupOffset.x + cNode.parentNode.hitbox.width;
            }
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
    
    /* -----------------------------------------------------------------------------------------------
     * ////////////////////////////////////////// DEBUGGING //////////////////////////////////////////
     * ----------------------------------------------------------------------------------------------- */
    
    // elkjs-exclude-start
    /**
     * For debugging. Writes hitboxes to svg.
     * 
     * @param name
     *          filename
     * @return
     *          this instance of {@link OneDimensionalCompactor}
     */
    public OneDimensionalCompactor drawHitboxes(final IElkProgressMonitor monitor, final String name) {
        //determine viewBox
        KVector topLeft = new KVector(Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY);
        KVector bottomRight = new KVector(Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY);
        for (CNode cNode : cGraph.cNodes) {
            topLeft.x = Math.min(topLeft.x, cNode.hitbox.x);
            topLeft.y = Math.min(topLeft.y, cNode.hitbox.y);
            bottomRight.x = Math.max(bottomRight.x, cNode.hitbox.x + cNode.hitbox.width);
            bottomRight.y = Math.max(bottomRight.y, cNode.hitbox.y + cNode.hitbox.height);
        }
        
        KVector size = bottomRight.clone().sub(topLeft);

        // drawing hitboxes to svg
        StringBuilder out = new StringBuilder();
        
        out.append("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
        out.append("<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"100%\" height=\"100%\""
                + "  viewBox=\"" + (topLeft.x) + " " + (topLeft.y) 
                + " " + size.x + " " + size.y + "\">\n");

        out.append("<defs><marker id=\"markerArrow\" markerWidth=\"13\" "
                + "markerHeight=\"13\" refX=\"2\" refY=\"6\" "
                + "orient=\"auto\">"
                + "  <path d=\"M2,2 L2,11 L10,6 L2,2\" style=\"fill: #000000;\" />"
                + "</marker></defs>\n");

        for (CNode cNode : cGraph.cNodes) {
            // the node's representation
            out.append(cNode.getDebugSVG()).append("\n");

            // the constraints
            for (CNode inc : cNode.constraints) {
                out.append("<line x1=\"" + (cNode.hitbox.x + cNode.hitbox.width / 2)
                        + "\" y1=\"" + (cNode.hitbox.y + cNode.hitbox.height / 2) + "\" x2=\""
                        + (inc.hitbox.x + inc.hitbox.width / 2) + "\" y2=\""
                        + (inc.hitbox.y + inc.hitbox.height / 2)
                        + "\" stroke=\"grey\" opacity=\"0.2\""
                        + " style=\"marker-start: url(#markerArrow);\" />\n");
            }
        }

        out.append("</svg>");

        monitor.logGraph(out.toString(), name, Type.SVG);
        
        return this;
    }
    // elkjs-exclude-end
    
}
