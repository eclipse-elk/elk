/*******************************************************************************
 * Copyright (c) 2022 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.alg.mrtree.intermediate;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Optional;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.stream.Collectors;

import org.eclipse.elk.alg.mrtree.TreeUtil;
import org.eclipse.elk.alg.mrtree.graph.TGraph;
import org.eclipse.elk.alg.mrtree.graph.TNode;
import org.eclipse.elk.alg.mrtree.options.EdgeRoutingMode;
import org.eclipse.elk.alg.mrtree.options.InternalProperties;
import org.eclipse.elk.alg.mrtree.options.MrTreeOptions;
import org.eclipse.elk.core.alg.ILayoutProcessor;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.options.Direction;
import org.eclipse.elk.core.util.IElkProgressMonitor;
import org.eclipse.elk.core.util.Pair;
import org.eclipse.elk.core.util.Triple;

/**
 * Applies one dimensional compaction to the graph.
 * 
 * @author jnc, sdo
 */
public class CompactionProcessor implements ILayoutProcessor<TGraph> {
    private List<Pair<Double, Double>> levels;
    
    /**
     * {@inheritDoc}
     */
    public void process(final TGraph tGraph, final IElkProgressMonitor progressMonitor) {
        progressMonitor.begin("Process compaction", 1);
        
        if (!tGraph.getProperty(MrTreeOptions.COMPACTION)) {
            return; // leave if option is not set
        }
        
        Direction dir = tGraph.getProperty(MrTreeOptions.DIRECTION);
        double nodeNodeSpacing = tGraph.getProperty(MrTreeOptions.SPACING_NODE_NODE);
        
        setUpLevels(tGraph, dir);
        
        computeNodeConstraints(tGraph, nodeNodeSpacing / 2 / 2);
        
        // Simple one dimensional compaction \w level preservation
        List<TNode> nodes = tGraph.getNodes();
        nodes.sort((x, y) -> Double.compare(
                TreeUtil.getDirectionVector(dir).dotProduct(new KVector(x.getPosition().x, x.getPosition().y)), 
                TreeUtil.getDirectionVector(dir).dotProduct(new KVector(y.getPosition().x, y.getPosition().y))));
        for (TNode n : nodes) {
            // Root nodes aren't compactable
            if (!n.getProperty(InternalProperties.ROOT)) {
                TNode d = getLowestDependentNode(n, dir);
                TNode p = TreeUtil.getLowestParent(n, tGraph);
                double newPos = 0, newPosSize = 0;
                if (d != null) {
                    // n has a dependent node
                    KVector pos = d.getPosition();
                    switch (dir) {
                    case LEFT:
                        newPos = pos.x - nodeNodeSpacing - n.getSize().x;
                        if (p.getPosition().x - nodeNodeSpacing - n.getSize().x < newPos) {
                            newPos = p.getPosition().x - nodeNodeSpacing - n.getSize().x;
                        }
                        newPosSize = newPos + n.getSize().x;
                        break;
                    case RIGHT:
                        newPos = pos.x + d.getSize().x + nodeNodeSpacing;
                        if (p.getPosition().x + nodeNodeSpacing > newPos) {
                            newPos = p.getPosition().x + p.getSize().x + nodeNodeSpacing;
                        }
                        newPosSize = newPos + n.getSize().x;
                        break;
                    case UP:
                        newPos = pos.y - nodeNodeSpacing - n.getSize().y;
                        if (p.getPosition().y - nodeNodeSpacing - n.getSize().y < newPos) {
                            newPos = p.getPosition().y - nodeNodeSpacing - n.getSize().y;
                        }
                        newPosSize = newPos + n.getSize().y;
                        break;
                    case DOWN:
                        newPos = pos.y + d.getSize().y + nodeNodeSpacing;
                        if (p.getPosition().y + nodeNodeSpacing > newPos) {
                            newPos = p.getPosition().y + p.getSize().y + nodeNodeSpacing;
                        }
                        newPosSize = newPos + n.getSize().y;
                        break;
                    }
                } else if (p != null) { 
                    // n does not have a dependent node but a parent
                    switch (dir) {
                    case LEFT:
                        newPos = p.getPosition().x - nodeNodeSpacing - n.getSize().x;
                        newPosSize = newPos + n.getSize().x;
                        break;
                    case RIGHT:
                        newPos = p.getPosition().x + p.getSize().x + nodeNodeSpacing;
                        newPosSize = newPos + n.getSize().x;
                        break;
                    case UP:
                        newPos = p.getPosition().y - nodeNodeSpacing - n.getSize().y;
                        newPosSize = newPos + n.getSize().y;
                        break;
                    case DOWN:
                        newPos = p.getPosition().y + p.getSize().y + nodeNodeSpacing;
                        newPosSize = newPos + n.getSize().y;
                        break;
                    }
                }
                
                if (tGraph.getProperty(MrTreeOptions.EDGE_ROUTING_MODE) == EdgeRoutingMode.AVOID_OVERLAP) {
                    final double finalNewPos = newPos, finalNewPosSize = newPosSize;
                    Optional<Pair<Double, Double>> level = levels.stream().
                            filter(x -> x.getFirst() <= finalNewPos && x.getSecond() >= finalNewPosSize).
                            findFirst();
                    if (level.isPresent()) { 
                        // The Node ended up within a level
                        if (dir.isHorizontal()) {
                            n.getPosition().x = newPos;
                        } else {
                            n.getPosition().y = newPos;
                        }
                    } else {
                        if (dir == Direction.LEFT || dir == Direction.UP) {
                            // Skip the first level as it only contains the SUPER_ROOT here
                            level = levels.stream().skip(1).
                                    filter(x -> x.getFirst() <= finalNewPos).findFirst();
                        } else {
                            // Skip the first level as it only contains the SUPER_ROOT here
                            level = levels.stream().skip(1).
                                    filter(x -> x.getFirst() >= finalNewPos).findFirst();
                        }
                        
                        // Force n into the found level
                        if (level.isPresent()) {
                            if (dir.isHorizontal()) {
                                n.getPosition().x = level.get().getFirst();
                            } else {
                                n.getPosition().y = level.get().getFirst();
                            }
                        }
                    }
                    
                    // Update Tree Level of node
                    if (level.isPresent()) {
                        int newIndex = levels.indexOf(level.get());
                        if (newIndex > 0 && newIndex != n.getProperty(MrTreeOptions.TREE_LEVEL)) {
                            n.setProperty(InternalProperties.COMPACT_LEVEL_ASCENSION, true);
                            n.setProperty(MrTreeOptions.TREE_LEVEL, newIndex);
                        }
                    }
                } else {
                    // In case of Aggressive Compaction just set the parameters
                    if (dir.isHorizontal()) {
                        n.getPosition().x = newPos;
                    } else {
                        n.getPosition().y = newPos;
                    }
                }
            }
        }
        progressMonitor.done();
    }
    
    /**
     * Sets up the bounds for each level.
     * 
     * @param tGraph the Graph the levels are in
     * @param dir the graphs layout direction
     */
    void setUpLevels(final TGraph tGraph, final Direction dir) {
        levels = new ArrayList<Pair<Double, Double>>();
        
        for (TNode n : tGraph.getNodes()) {
            // Adapt levels size to the current level
            while (n.getProperty(MrTreeOptions.TREE_LEVEL) > levels.size() - 1) {
                levels.add(new Pair<>(Double.MAX_VALUE, -Double.MAX_VALUE));
            }
            
            // Update level bounds
            int curLevel = n.getProperty(MrTreeOptions.TREE_LEVEL);
            if (dir.isHorizontal()) {
                if (n.getPosition().x < levels.get(curLevel).getFirst()) {
                    levels.get(curLevel).setFirst(n.getPosition().x);
                }
                if (n.getPosition().x + n.getSize().x > levels.get(curLevel).getSecond()) {
                    levels.get(curLevel).setSecond(n.getPosition().x + n.getSize().x);
                }
            } else {
                if (n.getPosition().y < levels.get(curLevel).getFirst()) {
                    levels.get(curLevel).setFirst(n.getPosition().y);
                }
                if (n.getPosition().y + n.getSize().y > levels.get(curLevel).getSecond()) {
                    levels.get(curLevel).setSecond(n.getPosition().y + n.getSize().y);
                }
            }
        }
    }
    
    /**
     * Compute constrains using scanline algorithm.
     * 
     * @param tGraph the Graph to work on
     */
    void computeNodeConstraints(final TGraph tGraph, final double nodeNodeSpacing) {
        Direction d = tGraph.getProperty(MrTreeOptions.DIRECTION);
        Direction right = d.isHorizontal() ? Direction.DOWN : Direction.RIGHT;
        
        // Get a filtered list of all relevant nodes
        List<TNode> actualNodes = tGraph.getNodes().stream().filter(x -> !x.getLabel().contains("SUPER_ROOT")).
                collect(Collectors.toList());
        // Build the point list by adding the upper left node edges, then the lower left ones and sorting all of them
        List<Triple<TNode, KVector, Boolean>> points = actualNodes.stream().
                map(x -> new Triple<TNode, KVector, Boolean>(x, x.getPosition().clone().
                        sub(nodeNodeSpacing, nodeNodeSpacing), true)).
                collect(Collectors.toList());
        points.addAll(actualNodes.stream().
                map(x -> new Triple<TNode, KVector, Boolean>(
                        x, x.getPosition().clone().add(x.getSize().x + nodeNodeSpacing, 
                                x.getSize().y + nodeNodeSpacing), false)).
                collect(Collectors.toList()));
        points.sort((x, y) -> Double.compare(
                TreeUtil.getDirectionVector(right).dotProduct(x.getSecond().clone()), 
                TreeUtil.getDirectionVector(right).dotProduct(y.getSecond().clone())));
        
        // Set up set sorted by the node positions projections onto a 1D space given by the direction vector
        SortedSet<TNode> s = new TreeSet<TNode>((x, y) -> Double.compare(
                TreeUtil.getDirectionVector(d).dotProduct(x.getPosition().clone()), 
                TreeUtil.getDirectionVector(d).dotProduct(y.getPosition().clone())));
        
        HashMap<TNode, TNode> cand = new HashMap<TNode, TNode>();
        
        // Scanline
        for (Triple<TNode, KVector, Boolean> p : points) {
            TNode r = p.getFirst();
            if (p.getThird()) { 
                s.add(r);
                if (s.headSet(r).size() > 0) {
                    cand.put(r, s.headSet(r).last());
                }
                if (s.tailSet(r).size() > 1) {
                    cand.put(getRightElement(s, r), r);
                }
            } else {
                // This deviates a bit from the proposed pseudo code but we need to check if right and left even exist
                if (s.headSet(r).size() > 0) {
                    TNode leftNode = s.headSet(r).last();
                    if (leftNode == cand.get(r)) {
                        r.getProperty(InternalProperties.COMPACT_CONSTRAINTS).add(leftNode); 
                    }
                }
                if (s.tailSet(r).size() > 1) {
                    TNode rightNode = getRightElement(s, r);
                    if (cand.get(rightNode) == r) {
                        rightNode.getProperty(InternalProperties.COMPACT_CONSTRAINTS).add(r);
                    }
                }
                s.remove(r);
            }
        }
    }
    
    /**
     * Gets the right element of n in the sorted set s.
     * 
     * @param s the sorted set
     * @param n the current node
     * @return the right element of n
     */
    TNode getRightElement(final SortedSet<TNode> s, final TNode n) {
        SortedSet<TNode> tailSet = s.tailSet(n);
        if (tailSet.size() <= 1) {
            throw new NullPointerException();
        }
        // Get the second element of the set using the iterator 
        Iterator<TNode> i = tailSet.iterator();
        i.next(); 
        return i.next();
    }
    
    /**
     * Returns the lowest Node in the graph in the given layout direction that is dependent on n.
     * 
     * @param n the current node
     * @return the lowest Node in the graph in the given layout direction that is dependent on n
     */
    TNode getLowestDependentNode(final TNode n, final Direction d) {
        List<TNode> cons = n.getProperty(InternalProperties.COMPACT_CONSTRAINTS);
        if (cons == null || cons.size() < 1) {
            return null;
        } else if (cons.size() == 1) {
            return cons.get(0);
        }
        
        TNode lowestCons = null;
        switch (d) {
        case LEFT:
            lowestCons = cons.stream().min((x, y) -> Double.compare(
                    x.getPosition().x, 
                    y.getPosition().x)).get();
            break;
        case RIGHT:
            lowestCons = cons.stream().max((x, y) -> Double.compare(
                    x.getPosition().x + x.getSize().x, 
                    y.getPosition().x + y.getSize().x)).get();
            break;
        case UP:
            lowestCons = cons.stream().min((x, y) -> Double.compare(
                    x.getPosition().y, 
                    y.getPosition().y)).get();
            break;
        case DOWN:
            lowestCons = cons.stream().max((x, y) -> Double.compare(
                    x.getPosition().y + x.getSize().y, 
                    y.getPosition().y + y.getSize().y)).get();
            break;
        }
        
        return lowestCons;
    }
}
