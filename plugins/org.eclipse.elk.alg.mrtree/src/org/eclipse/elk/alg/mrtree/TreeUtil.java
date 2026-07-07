/*******************************************************************************
 * Copyright (c) 2013 - 2022 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.mrtree;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.stream.Collectors;

import org.eclipse.elk.alg.mrtree.graph.TEdge;
import org.eclipse.elk.alg.mrtree.graph.TGraph;
import org.eclipse.elk.alg.mrtree.graph.TNode;
import org.eclipse.elk.alg.mrtree.options.InternalProperties;
import org.eclipse.elk.alg.mrtree.options.MrTreeOptions;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.options.Direction;

import com.google.common.collect.Iterables;

/**
 * Utility class for Mr Tree.
 * 
 * @author sgu
 * @author sdo
 */
public final class TreeUtil {

    /**
     * Hidden constructor to avoid instantiation.
     */
    private TreeUtil() {
    }
    
    /**
     * Gets the root of the graph.
     * @param tGraph the graph
     * @return the root
     */
    public static TNode getRoot(final TGraph tGraph) {
        return tGraph.getNodes().stream().
                filter(x -> x.getProperty(InternalProperties.ROOT)).
                findFirst().get();
    }
    
    /**
     * Gets the depth of the tree.
     * @param root the root
     * @return the depth of the elk tree
     */
    public static int depth(final TNode root) {
        List<TNode> childs = getChildren(root);
        if (childs.size() == 0) {
            return 1;
        } else {
            return childs.stream().map(x -> depth(x)).max(Integer::compare).get() + 1;
        }
    }
    
    /**
     * Gets the child nodes of n within a tree.
     * @param n the n
     * @return the child nodes of n within a tree
     */
    public static List<TNode> getChildren(final TNode n) {
        List<TNode> re = new ArrayList<TNode>();
        for (TEdge out : n.getOutgoingEdges()) {
            re.add(out.getTarget());
        }
        re = re.stream().distinct().collect(Collectors.toList());
        return re;
    }
    
    /**
     * Gets the distance from n to the root in a tree.
     * @param n the n
     * @param root the root
     * @return the distance from n to the root in an elk tree
     */
    public static int rootDistance(final TNode n, final TNode root) {
        if (n == root) {
            return 0;
        }
        // We should be in a tree so we assume that there is only one parent
        TNode p = (TNode) n.getIncomingEdges().get(0).getSource(); 
        return rootDistance(p, root) + 1;
    }
    
    /**
     * Gets all edges from a graph with n as their target.
     * The difference between this method and n.getIncomingEdges() is that the latter may in some cases be outdated
     * and not contain certain edges.
     * @param n the node we want the incoming edges from
     * @param graph the graph the edges are in
     * @return all the edges that go into n
     */
    public static List<TEdge> getAllIncomingEdges(final TNode n, final TGraph graph) {
        List<TEdge> re = new ArrayList<TEdge>(); 
        
        for (TEdge e : graph.getEdges()) {
            if (e.getTarget().id == n.id 
                    && e.getSource().getProperty(MrTreeOptions.TREE_LEVEL) 
                    != e.getTarget().getProperty(MrTreeOptions.TREE_LEVEL) 
                    && !re.stream().anyMatch(x -> x.toString().equals(e.toString()))) {
                re.add(e);
            }
        }
        
        re.sort((x, y) -> Double.compare(x.getSource().getPosition().x, y.getSource().getPosition().x));
        return re;
    }
    
    /**
     * Returns all nodes including n of the subtree which has n as its root.
     * 
     * @param n The root node of the subtree
     * @return All nodes of the subtree
     */
    public static List<TNode> getSubtree(final TNode n) {
        if (TreeUtil.getChildren(n).size() == 0) {
            List<TNode> re = new ArrayList<TNode>();
            re.add(n);
            return re;
        }
        
        List<TNode> re = TreeUtil.getChildren(n).stream().map(x -> getSubtree(x)).reduce((x, y) -> {
            List<TNode> ree = new ArrayList<TNode>();
            ree.addAll(x);
            ree.addAll(y);
            return ree;
        }).get();
        re.add(n);
        return re;
    }
    
    /**
     * Gets all edges from a graph with n as their source.
     * The difference between this method and n.getOutgoingEdges() is that the latter may in some cases be outdated
     * and not contain certain edges.
     * 
     * @param n the node we want the outgoing edges of
     * @param graph the graph the edges are in
     * @return all the edges that come out of n
     */
    public static List<TEdge> getAllOutgoingEdges(final TNode n, final TGraph graph) {
        List<TEdge> re = new ArrayList<TEdge>();
        
        for (TEdge e : graph.getEdges()) {
            if (e.getSource().id == n.id && !e.getSource().getLabel().equals("SUPER_ROOT") 
                    && e.getSource().getProperty(MrTreeOptions.TREE_LEVEL) 
                        != e.getTarget().getProperty(MrTreeOptions.TREE_LEVEL)
                    && !re.stream().anyMatch(x -> x.toString().equals(e.toString()))) {
                re.add(e);
            }
        }
        
        re.sort((x, y) -> Double.compare(x.getTarget().getPosition().x, y.getTarget().getPosition().x));
        return re;
    }
    
    /**
     * Gets the first point the edge will currently be routed towards.
     * 
     * @param e the edge
     * @return the first point the edge will currently be routed towards
     */
    public static KVector getFirstPoint(final TEdge e) {
        if (e.getBendPoints().size() == 0) {
            return new KVector(e.getTarget().getPosition().x, 
                    e.getTarget().getPosition().y);
        } else {
            return e.getBendPoints().getFirst();
        }
    }
    
    /**
     * Gets the last point the edge will currently be routed towards.
     * 
     * @param e the edge
     * @return the last point the edge will currently be routed towards
     */
    public static KVector getLastPoint(final TEdge e) {
        if (e.getBendPoints().size() == 0) {
            return new KVector(e.getSource().getPosition().x, 
                    e.getSource().getPosition().y);
        } else {
            return e.getBendPoints().getLast();
        }
    }
    
    /**
     * Gets the layout direction of the graph.
     * 
     * @param g the graph
     * @return the first point the edge will currently be routed towards
     */
    public static Direction getDirection(final TGraph g) {
        return g.getProperty(MrTreeOptions.DIRECTION);
    }
    
    /**
     * Gets the unit vector facing the specified direction.
     * @param d the direction
     * @return the unit vector facing the specified direction
     */
    public static KVector getDirectionVector(final Direction d) {
        switch (d) {
        case UP:
            return new KVector(0, -1);
        case RIGHT:
            return new KVector(1, 0);
        case LEFT:
            return new KVector(-1, 0);
        default:
            return new KVector(0, 1);
        }
    }
    
    /**
     * Returns the node size from the middle of the node in the specified Direction.
     * 
     * @param n the node
     * @param d the direction
     * @return the node size in the specified Direction
     */
    public static double getNodeSizeInDirection(final TNode n, final Direction d) {
        if (d == Direction.LEFT) {
            return -n.getSize().x / 2;
        } else if (d == Direction.UP) {
            return -n.getSize().y / 2;
        } else if (d == Direction.RIGHT) {
            return n.getSize().x / 2;
        } else {
            return n.getSize().y / 2;
        }
    }
    
    /**
     * Returns the node size vector from the middle of the node in the specified Direction.
     * 
     * @param n the node
     * @param d the direction
     * @return the node size vector in the specified Direction
     */
    public static KVector getNodeSizeVectorInDirection(final TNode n, final Direction d) {
        if (d == Direction.LEFT) {
            return new KVector(-n.getSize().x / 2, 0);
        } else if (d == Direction.UP) {
            return new KVector(0, -n.getSize().y / 2);
        } else if (d == Direction.RIGHT) {
            return new KVector(n.getSize().x / 2, 0);
        } else {
            return new KVector(0, n.getSize().y / 2);
        }
    }
    
    /**
     * Returns the Direction which is 90 degrees right of d.
     * 
     * @param d the input direction
     * @return the Direction which is 90 degrees right of d
     */
    public static Direction turnRight(final Direction d) {
        if (d == Direction.LEFT) {
            return Direction.UP;
        } else if (d == Direction.UP) {
            return Direction.RIGHT;
        } else if (d == Direction.RIGHT) {
            return Direction.DOWN;
        } else {
            return Direction.LEFT;
        }
    }
    
    /**
     * Returns the Direction which is 90 degrees left of d.
     * 
     * @param d the input direction
     * @return the Direction which is 90 degrees left of d
     */
    public static Direction turnLeft(final Direction d) {
        if (d == Direction.RIGHT) {
            return Direction.UP;
        } else if (d == Direction.UP) {
            return Direction.LEFT;
        } else if (d == Direction.LEFT) {
            return Direction.DOWN;
        } else {
            return Direction.RIGHT;
        }
    } 
    
    /**
     * Sets center outside of the Node.
     * 
     * @param center The Vector to mutate
     * @param next The next bendpoint after center
     * @param size The size of the node
     */
    public static void toNodeBorder(final KVector center, final KVector next, final KVector size) {
        double wh = size.x / 2, hh = size.y / 2;
        double absx = Math.abs(next.x - center.x), absy = Math.abs(next.y - center.y);
        double xscale = 1, yscale = 1;
        if (absx > wh) {
            xscale = wh / absx;
        }
        if (absy > hh) {
            yscale = hh / absy;
        }
        double scale = Math.min(xscale, yscale);
        center.x += scale * (next.x - center.x);
        center.y += scale * (next.y - center.y);
    }
    
    /**
     * Checks whether the edge induces a cycle into the graph.
     * 
     * @param g the graph
     * @param e the edge
     * @return true if the edge goes against the layout direction
     */
    public static boolean isCycleInducing(final TEdge e, final TGraph g) {
        return getDirectionVector(getDirection(g)).
                dotProduct(new KVector(e.getTarget().getPosition().x - e.getSource().getPosition().x, 
                        e.getTarget().getPosition().y - e.getSource().getPosition().y)) <= 0;
    }
    
    /**
     * Returns a unique long for two integers by placing their bits after each other into the long.
     * 
     * @param a the first 32 bit number
     * @param b the second 32 bit number
     * @return the unique key
     */
    public static long getUniqueLong(final int a, final int b) {
        return (long) a << 32 | b & 0xFFFFFFFFL; // SUPPRESS CHECKSTYLE MagicNumber
    }
    
    /**
     * Gets the lowest parent of the graph. 
     * Low is in this context defined as a coordinate with a high dot product with the direction vector of the graph.
     * 
     * @param n the node to get the parent of
     * @param g the graph the parent is in
     * @return the parent
     */
    public static TNode getLowestParent(final TNode n, final TGraph g) {
        KVector dirVec = getDirectionVector(getDirection(g));
        if (n.getIncomingEdges().size() == 0) {
            return null;
        }
        List<TNode> sources = n.getIncomingEdges().stream().map(x -> x.getSource()).collect(Collectors.toList());
        List<TNode> parents = g.getNodes().stream().filter(x -> sources.contains(x)).collect(Collectors.toList());
        Double lowestParentPos = parents.stream().
            map(x -> new KVector(x.getPosition().x + x.getSize().x / 2, 
                    x.getPosition().y + x.getSize().y / 2).dotProduct(dirVec)).
            max(Comparator.naturalOrder()).get();
        TNode lowestParent = parents.stream().
            filter(x -> new KVector(x.getPosition().x + x.getSize().x / 2, x.getPosition().y + x.getSize().y / 2).
                    dotProduct(dirVec) == lowestParentPos).
            findFirst().get();
            
        return lowestParent;
    }
    
    /**
     * This method returns the leftmost node at the deepest level. This is implemented using a
     * postorder walk of the subtree under given level.
     * 
     * @param currentlevel
     *            a list of nodes at level one
     * @return the leftmost node at the deepest level.
     */
    public static TNode getLeftMost(final Iterable<TNode> currentlevel) {
        return getLeftMost(currentlevel, -1);
    }

    /**
     * This method returns the leftmost node at the given level. This is implemented using a
     * postorder walk of the subtree under given level, depth levels down. Depth here refers to the
     * level below where the leftmost descendant is being found.
     * 
     * If given level is negative it returns the leftmost node at the deepest level.
     * 
     * @param currentlevel
     *            a list of nodes at level one
     * @param depth
     *            the depth to search for
     * @return the leftmost descendant at depth levels down
     */
    public static TNode getLeftMost(final Iterable<TNode> currentlevel, final int depth) {
        if (0 < Iterables.size(currentlevel)) {
            int d = depth;

            // the leftmost descendant at depth levels down
            if (1 < d) {
                d--;
                // build empty iterator
                Iterable<TNode> nextLevel = new Iterable<TNode>() {

                    public Iterator<TNode> iterator() {
                        return Collections.emptyIterator();
                    }
                };

                for (TNode cN : currentlevel) {
                    // append the children of the current node to the next level
                    nextLevel = Iterables.concat(nextLevel, cN.getChildren());
                }
                return getLeftMost(nextLevel, d);
            }

            // the leftmost node at the deepest level
            if (d < 0) {
                // build empty iterator
                Iterable<TNode> nextLevel = new Iterable<TNode>() {

                    public Iterator<TNode> iterator() {
                        return Collections.emptyIterator();
                    }
                };

                for (TNode cN : currentlevel) {
                    // append the children of the current node to the next level
                    nextLevel = Iterables.concat(nextLevel, cN.getChildren());
                }

                //
                if (0 < Iterables.size(nextLevel)) {
                    return getLeftMost(nextLevel, d);
                }
            }
        }
        // return the leftmost node at the current level
        return Iterables.getFirst(currentlevel, null);
    }

}
