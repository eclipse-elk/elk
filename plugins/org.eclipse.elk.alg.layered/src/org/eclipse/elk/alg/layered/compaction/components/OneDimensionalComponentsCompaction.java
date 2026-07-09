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
package org.eclipse.elk.alg.layered.compaction.components;

import java.util.List;
import java.util.Map.Entry;
import java.util.Set;

import org.eclipse.elk.alg.layered.compaction.oned.CGraph;
import org.eclipse.elk.alg.layered.compaction.oned.CGroup;
import org.eclipse.elk.alg.layered.compaction.oned.CNode;
import org.eclipse.elk.alg.layered.compaction.oned.OneDimensionalCompactor;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.options.Direction;
import org.eclipse.elk.core.util.IElkProgressMonitor;
import org.eclipse.elk.core.util.Pair;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;

/**
 * Compacts connected components using a {@link OneDimensionalCompactor}. Consecutively applies
 * horizontal and vertical compaction until no improvement can be made anymore.
 * 
 * @author uru
 * @param <N>
 *            the type of the contained nodes
 * @param <E>
 *            the type of the contained edges
 */
public final class OneDimensionalComponentsCompaction<N, E> {

    private CGraph compactionGraph;

    private ComponentsToCGraphTransformer<N, E> transformer;
    
    private List<Pair<CGroup, CNode>> verticalExternalExtensions;
    private List<Pair<CGroup, CNode>> horizontalExternalExtensions;
    
    private OneDimensionalCompactor compactor;
    
    private KVector topLeft;
    private KVector bottomRight;
    
    private static final int MAX_ITERATION = 10;
    private static final double EPSILON = 0.0001;
    
    /** To differentiate horizontal and vertical directions. */
    private static enum Dir {
        HORZ,
        VERT
    }
    private static final Set<Direction> LEFT_RIGHT = Sets.newHashSet(Direction.LEFT, Direction.RIGHT);
    private static final Set<Direction> UP_DOWN = Sets.newHashSet(Direction.UP, Direction.DOWN);
    
    /** Use {@link #init(IConnectedComponents, double)}. */
    private OneDimensionalComponentsCompaction() { }

    /**
     * Creates a new compaction object for the passed elements.
     * 
     * @param ccs
     *            a list with components to be compacted.
     * @param spacing
     *            the spacing to be preserved between any pair of components.
     * @return a new compactor object.
     * @param <N>
     *            the type of the contained nodes
     * @param <E>
     *            the type of the contained edges
     */
    public static <N, E> OneDimensionalComponentsCompaction<N, E> init(
            final IConnectedComponents<N, E> ccs, final double spacing) {

        OneDimensionalComponentsCompaction<N, E> compaction =
                new OneDimensionalComponentsCompaction<N, E>();
        compaction.transformer = new ComponentsToCGraphTransformer<N, E>(spacing);
        compaction.compactionGraph = compaction.transformer.transform(ccs);
        

        return compaction;
    }
    
    /**
     * Executes the compaction.
     * 
     * @param monitor
     *            a progress monitor.
     */
    public void compact(final IElkProgressMonitor monitor) {
        
        // separate vertical and horizontal extensions
        List<CNode> allNodes = Lists.newArrayList();
        verticalExternalExtensions = Lists.newArrayList();
        horizontalExternalExtensions = Lists.newArrayList();

        for (Entry<IExternalExtension<E>, Pair<CGroup, CNode>> entry : transformer
                .getExternalExtensions().entrySet()) {
            allNodes.add(entry.getValue().getSecond());
            if (entry.getKey().getDirection().isHorizontal()) {
                horizontalExternalExtensions.add(entry.getValue());
            } else {
                verticalExternalExtensions.add(entry.getValue());
            }
        }
        
        // add the external edge representations prior to the creation 
        // of the one dimensional compactor.
        // it is necessary for all group offsets to be calculated properly
        addExternalEdgeRepresentations(horizontalExternalExtensions);
        addExternalEdgeRepresentations(verticalExternalExtensions);

        // create our favorite compactor instance
        compactor = new OneDimensionalCompactor(compactionGraph);
        compactor.setSpacingsHandler(ComponentsToCGraphTransformer.SPACING_HANDLER);

        removeExternalEdgeRepresentations(horizontalExternalExtensions);
        removeExternalEdgeRepresentations(verticalExternalExtensions);
        
        allNodes.addAll(compactor.cGraph.cNodes);
        
        // calculating new graph size and offset
        topLeft = new KVector(Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY);
        bottomRight = new KVector(Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY);
        
        // determine the current drawing's dimensions
        for (CNode cNode : allNodes) {
            topLeft.x = Math.min(topLeft.x, cNode.hitbox.x);
            topLeft.y = Math.min(topLeft.y, cNode.hitbox.y);
            bottomRight.x = Math.max(bottomRight.x, cNode.hitbox.x + cNode.hitbox.width);
            bottomRight.y = Math.max(bottomRight.y, cNode.hitbox.y + cNode.hitbox.height);
        }

        // we "lock" nodes based on outgoing constraints, i.e. if there's a 
        // node in the "opposite" shadow we do not move it away
        // remember do execute an initial unlocked compaction 
        // in #compact() first, though
        compactor.setLockingStrategy((pair) -> pair.getFirst().cGroup.outDegree != 0);
        
        // now execute compaction until no improvement is made or we hit the maximum number of iterations
        int run = 0;
        double delta;
        do {
            delta = compact(run);
            run++;
        } while ((run < 2 || delta > EPSILON) && run < MAX_ITERATION);

        // "align" the compaction top left
        compactor.setLockingStrategy((pair) -> 
                    pair.getFirst().lock.get(pair.getSecond())
                || (pair.getFirst().cGroup.outDegree != 0 
                    && pair.getFirst().lock.get(pair.getSecond())));
        compact(0);
        
        // finish it!
        // make sure that the external edge representations have been removed from the 
        // compaction graph at this point for correct graph size calculation
        compactor.finish();
        
        // finally apply the determined positions 
        transformer.applyLayout();
    }

    /**
     * Internal compaction method, performs 3 compactions into each direction (horizontal first,
     * then vertical). The first compaction is unconstrained, the latter two consider the node locks.
     */
    private double compact(final int run) {
        
        double delta = 0;
        
        // reset all groups' deltas
        for (CGroup g : compactionGraph.cGroups) {
            g.delta = 0;
            g.deltaNormalized = 0;
        }
        
        // ----------------------------------------------------
        // #1 We want to perform horizontal compaction first.
        //    For this we have to add the vertical external 
        //    edges to the hulls to prevent edge node overlaps
        // ----------------------------------------------------
        
        addPlaceholders(Dir.HORZ);
        addExternalEdgeRepresentations(verticalExternalExtensions);
        compactor.calculateGroupOffsets();
        // there are new nodes in the compaction graph, we have to update the constraints
        compactor.forceConstraintsRecalculation();
        
        Direction direction = Direction.LEFT;
        // compact horizontally
        compactor
            .changeDirection(direction)
            .compact()
            .changeDirection(direction.opposite())
            .applyLockingStrategy()
            .compact()
            .changeDirection(direction)
            .applyLockingStrategy()
            .compact();
        
        // very important to transform back to LEFT 
        // because the hitboxes representing external edges 
        // that will be added in a second have not been transformed yet
        compactor.changeDirection(Direction.LEFT);

        // remove the vertical external edges again
        removeExternalEdgeRepresentations(verticalExternalExtensions);
        removePlaceholders(Dir.HORZ);

        // .. and offset the horizontal external edges (that will be added next)
        // according to the just finished horizontal compaction
        updateExternalExtensionDimensions(Dir.HORZ);
        updatePlaceholders(Dir.VERT);
        
        // ----------------------------------------------------
        // #2 We want to perform vertical compaction.
        // ----------------------------------------------------
        
        // now add them
        addPlaceholders(Dir.VERT);
        addExternalEdgeRepresentations(horizontalExternalExtensions);
        compactor.calculateGroupOffsets();
        // remember delta from horizontal compaction
        for (CGroup g : compactionGraph.cGroups) {
            delta += Math.abs(g.deltaNormalized);
        }
        // reset all groups' deltas
        for (CGroup g : compactionGraph.cGroups) {
            g.delta = 0;
            g.deltaNormalized = 0;
        }

        direction = Direction.UP;
        // compact vertically
        compactor
            .changeDirection(direction)
            .forceConstraintsRecalculation()
            .compact()
            .changeDirection(direction.opposite())
            .applyLockingStrategy()
            .compact()
            .changeDirection(direction)
            .applyLockingStrategy()
            .compact();
        
        // transform back to LEFT
        compactor.changeDirection(Direction.LEFT);
        
        // remove the horizontal external edges
        removeExternalEdgeRepresentations(horizontalExternalExtensions);
        removePlaceholders(Dir.VERT);
        
        // ... offset the vertical external edges 
        // (which have been excluded during the last compaction)
        updateExternalExtensionDimensions(Dir.VERT);
        updatePlaceholders(Dir.HORZ);
        
        compactor.forceConstraintsRecalculation();
        
        // ... and the delta from vertical compaction
        for (CGroup g : compactionGraph.cGroups) {
            delta += Math.abs(g.deltaNormalized);
        }
        
        return delta;
    }
    
    /**
     * @return After {@link #compact(IKielerProgressMonitor)} has been executed, returns the new
     *         size of the graph. That is, the size of a rectangle that represents the bounding box
     *         of all contained {@link IComponent}s.
     * 
     *         Note that this bounding box considers only real nodes and no contribution of external
     *         edges.
     */
    public KVector getGraphSize() {
        return transformer.getGraphSize();
    }

    /**
     * @param c
     *            the {@link IComponent} for which the offset is requested.
     * @return After {@link #compact(IKielerProgressMonitor)} has been executed returns the offset
     *         of this {@link IComponent} compared to the component's position prior to compaction.
     */
    public KVector getOffset(final IComponent<N, E> c) {
        KVector individual = transformer.getOffset(c);
        return individual.negate().add(transformer.getGlobalOffset());
    }

    /**
     * Adds the passed external edges to the {@link #compactionGraph}.
     */
    private void addExternalEdgeRepresentations(final List<Pair<CGroup, CNode>> ees) {
        for (Pair<CGroup, CNode> p : ees) {
            compactionGraph.cNodes.add(p.getSecond());
            p.getFirst().addCNode(p.getSecond());
        }
    }
    
    /**
     * Remove the passed external edges from the {@link #compactionGraph}.
     */
    private void removeExternalEdgeRepresentations(final List<Pair<CGroup, CNode>> ees) {
        for (Pair<CGroup, CNode> p : ees) {
            compactionGraph.cNodes.remove(p.getSecond());
            p.getFirst().removeCNode(p.getSecond());
        }
    }
    
    private void addPlaceholders(final Dir dir) {
        Set<Direction> dirs = (dir == Dir.VERT) ? UP_DOWN : LEFT_RIGHT;
        for (Direction d : dirs) {
            for (Pair<CGroup, CNode> pair : transformer.getExternalPlaceholder().get(d)) {
                compactionGraph.cNodes.add(pair.getSecond());
                compactionGraph.cGroups.add(pair.getSecond().cGroup);
            }
        }
    }
    
    private void removePlaceholders(final Dir dir) {
        Set<Direction> dirs = (dir == Dir.VERT) ? UP_DOWN : LEFT_RIGHT;
        for (Direction d : dirs) {
            for (Pair<CGroup, CNode> pair : transformer.getExternalPlaceholder().get(d)) {
                compactionGraph.cNodes.remove(pair.getSecond());
                compactionGraph.cGroups.remove(pair.getSecond().cGroup);
            }
        }
    }
    
    private void updatePlaceholders(final Dir dir) {
        Set<Direction> dirs = (dir == Dir.VERT) ? UP_DOWN : LEFT_RIGHT;
        for (Direction d : dirs) {
            for (Pair<CGroup, CNode> pair : transformer.getExternalPlaceholder().get(d)) {
                CNode cNode = pair.getSecond();
                CGroup parentComponentGroup = pair.getFirst();
                
                // deltaNormalized is negative if a group was moved to the left,
                //  positive if the group was moved to the right
                // the size of each component stays the same!
                double adelta = parentComponentGroup.deltaNormalized;
                switch (d) {
                    case LEFT:
                    case RIGHT:
                        cNode.hitbox.y += adelta;
                        break;
        
                    case UP:
                    case DOWN:
                        // y sticks to the top of the diagram
                        cNode.hitbox.x += adelta;
                        break;
                }
            }
        }
    }
    
    /**
     * During compaction of the x dimension the horizontal external extensions are not part of the
     * graph that is being compacted. As a consequence, their positions are not updated. We make up
     * for this in this method. Additionally, the dimensions of the diagram are respected and
     * external extensions are guaranteed to always extend to the very border of the diagram
     * (preventing components from overlapping with edges).
     * 
     * @param dir
     *            whether to update the vertical or the horizontal extensions.
     */
    private void updateExternalExtensionDimensions(final Dir dir) {

        for (Entry<IExternalExtension<E>, Pair<CGroup, CNode>> entry : transformer
                .getExternalExtensions().entrySet()) {

            IExternalExtension<E> ee = entry.getKey();

            if (dir == Dir.VERT) {
                if (ee.getDirection() != Direction.UP && ee.getDirection() != Direction.DOWN) {
                    continue;
                }
            } else {
                if (ee.getDirection() != Direction.LEFT && ee.getDirection() != Direction.RIGHT) {
                    continue;
                }
            }

            CNode cNode = entry.getValue().getSecond();
            // the node may not be associated with the group atm, thus cNode.group may be null
            CGroup group = entry.getValue().getFirst();

            // deltaNormalized is negative if a group was moved to the left,
            // positive if the group was moved to the right
            double adelta = group.deltaNormalized;
            switch (ee.getDirection()) {
                case LEFT:
                    // x sticks to the left of the diagram
                    cNode.hitbox.x = topLeft.x;
                    cNode.hitbox.width = Math.max(1, cNode.hitbox.width + adelta);
                    break;
                case RIGHT:
                    cNode.hitbox.x = cNode.hitbox.x + adelta;
                    cNode.hitbox.width = Math.max(1, cNode.hitbox.width - adelta);
                    break;
    
                case UP:
                    // y sticks to the top of the diagram
                    cNode.hitbox.y = topLeft.y;
                    cNode.hitbox.height = Math.max(1, cNode.hitbox.height + adelta);
                    break;
                case DOWN:
                    cNode.hitbox.y = cNode.hitbox.y + adelta;
                    cNode.hitbox.height = Math.max(1,  cNode.hitbox.height - adelta);
                    break;
            }
        }
        
    }

}
