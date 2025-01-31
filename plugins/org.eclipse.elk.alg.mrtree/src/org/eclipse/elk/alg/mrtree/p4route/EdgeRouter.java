/*******************************************************************************
 * Copyright (c) 2013 - 2022 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.alg.mrtree.p4route;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.stream.Collectors;

import org.eclipse.elk.alg.mrtree.TreeLayoutPhases;
import org.eclipse.elk.alg.mrtree.TreeUtil;
import org.eclipse.elk.alg.mrtree.graph.TEdge;
import org.eclipse.elk.alg.mrtree.graph.TGraph;
import org.eclipse.elk.alg.mrtree.graph.TNode;
import org.eclipse.elk.alg.mrtree.intermediate.IntermediateProcessorStrategy;
import org.eclipse.elk.alg.mrtree.options.EdgeRoutingMode;
import org.eclipse.elk.alg.mrtree.options.InternalProperties;
import org.eclipse.elk.alg.mrtree.options.MrTreeOptions;
import org.eclipse.elk.core.alg.ILayoutPhase;
import org.eclipse.elk.core.alg.LayoutProcessorConfiguration;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.math.KVectorChain;
import org.eclipse.elk.core.options.Direction;
import org.eclipse.elk.core.util.IElkProgressMonitor;
import org.eclipse.elk.core.util.Pair;
import org.eclipse.elk.core.util.Triple;

/**
 * This class implements a dull edge routing by setting just source and target of a edge in mode MiddleToMiddle 
 * as well as a edge routing algorithm that is designed to prevent edge node and edge edge overlaps in more complex 
 * graphs.
 * 
 * @author sor
 * @author sgu
 * @author jnc
 * @author sdo
 * 
 */
public class EdgeRouter implements ILayoutPhase<TreeLayoutPhases, TGraph> {
    
    private final double oneHalf = 0.5;
    private final double steepEndEdgeTheresholdDistance = 50;
    private final double steepEndEdgeRatio = 5.3;
    private final double steepEndEdgeSampleHeight = 40;

    /** intermediate processing configuration. */
    private static final LayoutProcessorConfiguration<TreeLayoutPhases, TGraph> INTERMEDIATE_PROCESSING_CONFIG =
            LayoutProcessorConfiguration.<TreeLayoutPhases, TGraph>create()
                    .before(TreeLayoutPhases.P4_EDGE_ROUTING)
                        .add(IntermediateProcessorStrategy.LEVEL_COORDS)
                        .add(IntermediateProcessorStrategy.COMPACTION_PROC)
                        .add(IntermediateProcessorStrategy.GRAPH_BOUNDS_PROC);

    /**
     * {@inheritDoc}
     */
    @Override
    public LayoutProcessorConfiguration<TreeLayoutPhases, TGraph> getLayoutProcessorConfiguration(final TGraph graph) {
        return INTERMEDIATE_PROCESSING_CONFIG;
    }
    
    /**
     * {@inheritDoc}
     */
    public void process(final TGraph tGraph, final IElkProgressMonitor progressMonitor) {
        progressMonitor.begin("Edge routing", 1);

        EdgeRoutingMode mode = tGraph.getProperty(MrTreeOptions.EDGE_ROUTING_MODE);
        if (mode == EdgeRoutingMode.MIDDLE_TO_MIDDLE) {
            middleToMiddle(tGraph);
        } else if (mode == EdgeRoutingMode.AVOID_OVERLAP) {
            avoidOverlap(tGraph);
            
            // fallback
            for (TEdge e : tGraph.getEdges()) {
                if (e.getBendPoints().size() < 2) {
                    middleToMiddleEdgeRoute(e);
                }
            }
        }
        
        progressMonitor.done();
    }
    
    /**
     * Handles MiddleToMiddle edge routing.
     * 
     * @param tGraph the graph to work on
     */
    private void middleToMiddle(final TGraph tGraph) {
        for (TEdge tEdge : tGraph.getEdges()) {
            middleToMiddleEdgeRoute(tEdge);
        }
    }
    
    /**
     * Routes a single edge from the middle of one node to the other.
     * 
     * @param tEdge the edge to route
     */
    private void middleToMiddleEdgeRoute(final TEdge tEdge) {
        KVectorChain bendPoints = tEdge.getBendPoints();

        // add the source port and target points to the vector chain
        TNode source = tEdge.getSource();
        TNode target = tEdge.getTarget();
        
        KVector sourcePoint = new KVector(source.getPosition().x + source.getSize().x / 2, 
                source.getPosition().y + source.getSize().y / 2);
        KVector targetPoint = new KVector(target.getPosition().x + target.getSize().x / 2, 
                target.getPosition().y + target.getSize().y / 2);
        
        bendPoints.addFirst(sourcePoint);
        bendPoints.addLast(targetPoint);
        
        // correct the source and target points
        TreeUtil.toNodeBorder(sourcePoint, bendPoints.get(1),
                tEdge.getSource().getSize());
        TreeUtil.toNodeBorder(targetPoint, bendPoints.get(bendPoints.size() - 2),
                tEdge.getTarget().getSize());
    }
    
    /**
     * Handles AvoidOverlap Edge Routing.
     * 
     * @param tGraph
     */
    private void avoidOverlap(final TGraph tGraph) {
        TNode root = TreeUtil.getRoot(tGraph);
        double nodeBendpointPadding = tGraph.getProperty(MrTreeOptions.SPACING_EDGE_NODE);
        double edgeEndTexturePadding = tGraph.getProperty(MrTreeOptions.EDGE_END_TEXTURE_LENGTH);
        Direction d = tGraph.getProperty(MrTreeOptions.DIRECTION);
        
        avoidOverlapSetStartPoints(tGraph, d, nodeBendpointPadding);
        
        avoidOverlapSpecialEdges(tGraph, root, d, nodeBendpointPadding, edgeEndTexturePadding);
        
        avoidOverlapSetEndPoints(tGraph, d, nodeBendpointPadding, edgeEndTexturePadding);
    }
    
    /**
     * Handles special edges such as multi level edges and cycle inducing edges for the AvoidOverlap router.
     * Multi level edges are edges that don't just go from one level to the next but through multiple levels.
     * Cycle inducing edges start lower in the tree than where they end.
     * 
     * @param tGraph the Graph to work on
     * @param root the Graph's root node
     * @param d the current layout direction
     * @param nodeBendpointPadding the Node BendPoint Padding
     * @param edgeEndTexturePadding the length of the texture at the end of edges set in the MrTreeOptions
     */
    private void avoidOverlapSpecialEdges(final TGraph tGraph, final TNode root, final Direction d, 
            final double nodeBendpointPadding, final double edgeEndTexturePadding) {
        int sideOneEdges = 0, sideTwoEdges = 0; // Counts how many edges are routed along the sides of the tree
        HashMap<Long, MultiLevelEdgeNodeNodeGap> nodeGaps = new HashMap<Long, MultiLevelEdgeNodeNodeGap>();
        
        int maxLevel = tGraph.getNodes().stream().
                map(x -> x.getProperty(MrTreeOptions.TREE_LEVEL)).
                max(Integer::compare).get() + 1;
        int[] outsPerLevel = new int[maxLevel];
        int[] insPerLevel = new int[maxLevel];
        for (int i = 0; i < maxLevel; i++) {
            outsPerLevel[i] = 0;
            insPerLevel[i] = 0;
        }
        
        List<TEdge> distictEdges = tGraph.getEdges().stream().distinct().collect(Collectors.toList());
        for (TEdge e : distictEdges) {
            int sourceLevel = e.getSource().getProperty(MrTreeOptions.TREE_LEVEL);
            int targetLevel = e.getTarget().getProperty(MrTreeOptions.TREE_LEVEL);
            int levelDiff = targetLevel - sourceLevel;
            if (levelDiff > 1) {
                // Multi level edges
                for (int curLevel = sourceLevel + 1; curLevel < targetLevel; curLevel++) {
                    final int finalCurlevel = curLevel;
                    List<TNode> nextLevelNodes = tGraph.getNodes().stream().filter(x -> 
                                    x.getProperty(MrTreeOptions.TREE_LEVEL) == finalCurlevel)
                                    .collect(Collectors.toList());
                    // Find the node gap in the next level through which we can route our multi level edge
                    int i = 0;
                    if (d.isHorizontal()) {
                        nextLevelNodes.sort((x, y) -> Double.compare(x.getPosition().y, y.getPosition().y));
                        for (i = 0; i < nextLevelNodes.size(); i++) {
                            double interpolation = (curLevel - sourceLevel) / (double) (targetLevel - sourceLevel);
                            if (nextLevelNodes.get(i).getPosition().y > (e.getSource().getPosition().y 
                                    * (1 - interpolation) + e.getTarget().getPosition().y * interpolation)) {
                                break;
                            }
                        }
                        
                        // Skip unnecessary level side bendPoints
                        if (nextLevelNodes.size() > 0) {
                            KVector start = e.getBendPoints().size() == 0 
                                    ? e.getSource().getPosition().clone() : e.getBendPoints().getLast();
                            KVector last = nextLevelNodes.get(nextLevelNodes.size() - 1).getPosition().clone().
                                    add(nextLevelNodes.get(nextLevelNodes.size() - 1).getSize());
                            KVector first = nextLevelNodes.get(0).getPosition().clone().
                                    add(nextLevelNodes.get(0).getSize());
                            if (i >= nextLevelNodes.size() - 1 && start.y > last.y 
                                    && e.getTarget().getPosition().y > last.y) {
                                continue;
                            }
                            if (i <= 0 && start.y < first.x 
                                    && e.getTarget().getPosition().y < first.y) {
                                continue;
                            }
                        }
                    } else {
                        nextLevelNodes.sort((x, y) -> Double.compare(x.getPosition().x, y.getPosition().x));
                        for (i = 0; i < nextLevelNodes.size(); i++) {
                            double interpolation = (curLevel - sourceLevel) / (double) (targetLevel - sourceLevel);
                            if (nextLevelNodes.get(i).getPosition().x > (e.getSource().getPosition().x 
                                    * (1 - interpolation) + e.getTarget().getPosition().x * interpolation)) {
                                break;
                            }
                        }
                        
                        // Skip unnecessary level side bendPoints
                        if (nextLevelNodes.size() > 0) {
                            KVector start = e.getBendPoints().size() == 0 
                                    ? e.getSource().getPosition().clone() : e.getBendPoints().getLast();
                            KVector last = nextLevelNodes.get(nextLevelNodes.size() - 1).getPosition().clone().
                                    add(nextLevelNodes.get(nextLevelNodes.size() - 1).getSize());
                            KVector first = nextLevelNodes.get(0).getPosition().clone().
                                    add(nextLevelNodes.get(0).getSize());
                            if (i >= nextLevelNodes.size() - 1 && start.x > last.x 
                                    && e.getTarget().getPosition().x > last.x) {
                                continue;
                            }
                            if (i <= 0 && start.x < first.x 
                                    && e.getTarget().getPosition().x < first.x) {
                                continue;
                            }
                        }
                    }
                    
                    // Add multi level edge gap to nodeGaps
                    KVector bend1 = new KVector(), bend2 = new KVector();
                    e.getBendPoints().add(bend1); e.getBendPoints().add(bend2);
                    MultiLevelEdgeNodeNodeGap gap;
                    Triple<KVector, KVector, TEdge> bendTriple = new Triple<KVector, KVector, TEdge>(bend1, bend2, e);
                    long key = TreeUtil.getUniqueLong(curLevel, i);
                    if (!nodeGaps.containsKey(key)) {
                        gap = new MultiLevelEdgeNodeNodeGap(i == 0 ? null : nextLevelNodes.get(i - 1), 
                                i == nextLevelNodes.size() ? null : nextLevelNodes.get(i), bendTriple, tGraph);
                        nodeGaps.put(key, gap);
                    } else {
                        gap = nodeGaps.get(key);
                        gap.addBendPoints(bendTriple);
                    }
                    
                    // Increment end of level edge counters
                    if (!d.isHorizontal()) {
                        if (gap.isOnFirstNodeSide() && gap.getNeighborTwo().getPosition().x 
                                <= tGraph.getProperty(InternalProperties.GRAPH_XMIN)) {
                            sideOneEdges++;
                        }
                        if (gap.isOnLastNodeSide() 
                                && gap.getNeighborOne().getPosition().x + gap.getNeighborOne().getSize().x
                                >= tGraph.getProperty(InternalProperties.GRAPH_XMAX)) {
                            sideTwoEdges++;
                        }
                    } else {
                        if (gap.isOnFirstNodeSide() && gap.getNeighborTwo().getPosition().y 
                                <= tGraph.getProperty(InternalProperties.GRAPH_YMIN)) {
                            sideOneEdges++;
                        }
                        if (gap.isOnLastNodeSide() 
                                && gap.getNeighborOne().getPosition().y + gap.getNeighborOne().getSize().y
                                >= tGraph.getProperty(InternalProperties.GRAPH_YMAX)) {
                            sideTwoEdges++;
                        }
                    }
                }
            } else if (levelDiff == 0) {
                // Fallback for edges that come from and go to the same level
                middleToMiddleEdgeRoute(e); 
            } else if (levelDiff < 0) {
                outsPerLevel[sourceLevel]++;
                insPerLevel[targetLevel]++;
                Pair<Integer, Integer> sides = avoidOverlapHandleCycleInducingEdges(e, d, tGraph, 
                        new Pair<Integer, Integer>(sideOneEdges, sideTwoEdges), nodeBendpointPadding, 
                        edgeEndTexturePadding, new Pair<Integer, Integer>(
                                insPerLevel[targetLevel], outsPerLevel[sourceLevel]));
                sideOneEdges = sides.getFirst();
                sideTwoEdges = sides.getSecond();
            }
        }
    }
    
    private Pair<Integer, Integer> avoidOverlapHandleCycleInducingEdges(final TEdge e, final Direction d, 
            final TGraph tGraph, final Pair<Integer, Integer> sideEdges, final double nodeBendpointPadding, 
            final double edgeEndTexturePadding, final Pair<Integer, Integer> inOuts) {
        int sideOneEdges = sideEdges.getFirst();
        int sideTwoEdges = sideEdges.getSecond();
        
        TNode s = e.getSource(), t = e.getTarget();
        // Compute bendTmp which is the coordinate at the end of the graph the edge will be routed along
        double bendTmp = 0, middleTree = 0;
        if (d.isHorizontal()) {
            middleTree = tGraph.getNodes().stream().
                map(x -> x.getPosition().y + x.getSize().y / 2).
                mapToDouble(Double::doubleValue).average().getAsDouble();
            if (s.getPosition().y + s.getSize().y / 2 > middleTree) {
                final int finalSideTwoEdges = ++sideTwoEdges;
                bendTmp = tGraph.getNodes().stream().
                        map(x -> x.getPosition().y + x.getSize().y + nodeBendpointPadding * finalSideTwoEdges).
                        max(Double::compare).get();
            } else {
                final int finalSideOneEdges = ++sideOneEdges;
                bendTmp = tGraph.getNodes().stream().
                        map(x -> x.getPosition().y - nodeBendpointPadding * finalSideOneEdges).
                        min(Double::compare).get();
            }
        } else {
            middleTree = tGraph.getNodes().stream().
                map(x -> x.getPosition().x + x.getSize().x / 2).
                mapToDouble(Double::doubleValue).average().getAsDouble();
            if (s.getPosition().x + s.getSize().x / 2 > middleTree) {
                final int finalSideTwoEdges = ++sideTwoEdges;
                bendTmp = tGraph.getNodes().stream().
                        map(x -> x.getPosition().x + x.getSize().x + nodeBendpointPadding * finalSideTwoEdges).
                        max(Double::compare).get();
            } else {
                final int finalSideOneEdges = ++sideOneEdges;
                bendTmp = tGraph.getNodes().stream().
                        map(x -> x.getPosition().x - nodeBendpointPadding * finalSideOneEdges).
                        min(Double::compare).get();
            }
        }
        // Set bendPoints
        if (d == Direction.LEFT) {
            e.getBendPoints().add(s.getProperty(InternalProperties.LEVELMIN) - nodeBendpointPadding, bendTmp);
            e.getBendPoints().add(t.getPosition().x + t.getSize().x 
                    + nodeBendpointPadding + edgeEndTexturePadding, bendTmp);
            e.getBendPoints().add(t.getPosition().x + t.getSize().x 
                    + nodeBendpointPadding + edgeEndTexturePadding, t.getPosition().y + t.getSize().y / 2);
            e.getBendPoints().add(t.getPosition().x + t.getSize().x, t.getPosition().y + t.getSize().y / 2);
        } else if (d == Direction.RIGHT) {
            e.getBendPoints().add(s.getProperty(InternalProperties.LEVELMAX) + nodeBendpointPadding, 
                    s.getPosition().y + s.getSize().y / 2);
            e.getBendPoints().add(s.getPosition().x + s.getSize().x + nodeBendpointPadding, bendTmp);
            e.getBendPoints().add(t.getPosition().x - nodeBendpointPadding - edgeEndTexturePadding, bendTmp);
            e.getBendPoints().add(t.getPosition().x - nodeBendpointPadding - edgeEndTexturePadding, 
                    t.getPosition().y + t.getSize().y / 2);
            e.getBendPoints().add(t.getPosition().x, t.getPosition().y + t.getSize().y / 2);
        } else if (d == Direction.UP) {
            e.getBendPoints().add(bendTmp, s.getProperty(InternalProperties.LEVELMIN) - nodeBendpointPadding);
            e.getBendPoints().add(bendTmp,
                    t.getPosition().y + t.getSize().y + nodeBendpointPadding + edgeEndTexturePadding);
            e.getBendPoints().add(t.getPosition().x + t.getSize().x / 2, 
                    t.getPosition().y + t.getSize().y + nodeBendpointPadding + edgeEndTexturePadding);
            e.getBendPoints().add(t.getPosition().x + t.getSize().x / 2, 
                    t.getPosition().y + t.getSize().y + nodeBendpointPadding);
        } else {
            if (!e.getBendPoints().isEmpty()) {
                e.getBendPoints().getLast().y = s.getProperty(InternalProperties.LEVELMAX) + nodeBendpointPadding 
                        * inOuts.getSecond();
            }
            e.getBendPoints().add(bendTmp, 
                    s.getProperty(InternalProperties.LEVELMAX) + nodeBendpointPadding * inOuts.getSecond());
            e.getBendPoints().add(bendTmp,
                    t.getPosition().y - nodeBendpointPadding * inOuts.getFirst() - edgeEndTexturePadding);
        }
        
        return new Pair<Integer, Integer>(sideOneEdges, sideTwoEdges);
    }

    
    /**
     * Handles the start points for AvoidOverlap edge routing.
     * 
     * @param tGraph the graph to work on
     * @param d the current layout direction
     * @param nodeBendpointPadding the Node BendPoint Padding 
     */
    private void avoidOverlapSetStartPoints(final TGraph tGraph, final Direction d, final double nodeBendpointPadding) {
        for (TNode n : tGraph.getNodes()) {
            if (n.getLabel().equals("SUPER_ROOT")) {
               continue;
            }
            
            // Get a list of all outgoing edges from n
            // I specifically do not use n.getOutgoingEdges() here because that list does not contain all outgoing 
            // edges in certain scenarios. For example if a node has two parents who are also both root nodes
            List<TEdge> outs = TreeUtil.getAllOutgoingEdges(n, tGraph); 
            if (d.isHorizontal()) {
                outs.sort((x, y) -> Double.compare(TreeUtil.getFirstPoint(x).y, 
                        TreeUtil.getFirstPoint(y).y));
            } else {
                outs.sort((x, y) -> Double.compare(TreeUtil.getFirstPoint(x).x, 
                        TreeUtil.getFirstPoint(y).x));
            }

            // Set the bendPoints for all outs
            int num = outs.size();
            for (int i = 0; i < num; i++) {
                TNode tar = outs.get(i).getTarget();
                if (tar.getLabel().equals("n11")) {
                    getClass();
                }
                if (n.getProperty(InternalProperties.COMPACT_LEVEL_ASCENSION)
                        && !TreeUtil.isCycleInducing(outs.get(i), tGraph)) {
                    continue;
                }
                
                double interpolation = num == 1 ? oneHalf : (i + 1) / (double) (num + 1);
                if (d == Direction.LEFT) {
                    double levelEndCoord = n.getProperty(InternalProperties.LEVELMIN);
                    double y = n.getPosition().y + n.getSize().y * interpolation;
                    outs.get(i).getBendPoints().addFirst(Math.min(levelEndCoord, 
                            n.getPosition().x - nodeBendpointPadding), y);
                    outs.get(i).getBendPoints().addFirst(n.getPosition().x, y);
                } else if (d == Direction.RIGHT) {
                    double levelEndCoord = n.getProperty(InternalProperties.LEVELMAX) + nodeBendpointPadding;
                    double y = n.getPosition().y + n.getSize().y * interpolation;
                    outs.get(i).getBendPoints().addFirst(levelEndCoord, y);
                    outs.get(i).getBendPoints().addFirst(n.getPosition().x + n.getSize().x, y);
                } else if (d == Direction.UP) {
                    double levelEndCoord = n.getProperty(InternalProperties.LEVELMIN);
                    double x = n.getPosition().x + n.getSize().x * interpolation;
                    outs.get(i).getBendPoints().addFirst(x, 
                            Math.min(n.getPosition().y - nodeBendpointPadding, levelEndCoord));
                    outs.get(i).getBendPoints().addFirst(x, n.getPosition().y);
                } else {
                    double levelEndCoord = n.getProperty(InternalProperties.LEVELMAX) + nodeBendpointPadding;
                    double x = n.getPosition().x + n.getSize().x * interpolation;
                    outs.get(i).getBendPoints().addFirst(x, levelEndCoord);
                    outs.get(i).getBendPoints().addFirst(x, n.getPosition().y + n.getSize().y);
                }
            }
        }
    }
    
    /**
     * Handles the end points for AvoidOverlap edge routing.
     * 
     * @param tGraph the graph to work on
     * @param d the current layout direction
     * @param nodeBendpointPadding the Node BendPoint Padding
     * @param edgeEndTexturePadding the length of the texture at the end of edges set in the MrTreeOptions 
     */
    private void avoidOverlapSetEndPoints(final TGraph tGraph, final Direction d, final double nodeBendpointPadding, 
            final double edgeEndTexturePadding) {
        for (TNode n : tGraph.getNodes()) {
            
            if (n.getLabel().equals("SUPER_ROOT")) {
                continue;
            }
            
            // Get the incoming edges and sort them by their current bendPoints
            // Similarly to setStartPoints I do not use the n.getIncomingEdges() here for the same reasons
            List<TEdge> ins = TreeUtil.getAllIncomingEdges(n, tGraph).stream().collect(Collectors.toList());
            if (d.isHorizontal()) {
                ins.sort((x, y) -> Double.compare(TreeUtil.getLastPoint(x).y, 
                        TreeUtil.getLastPoint(y).y));
            } else {
                ins.sort((x, y) -> Double.compare(TreeUtil.getLastPoint(x).x, 
                        TreeUtil.getLastPoint(y).x));
            }
            
            // Set the bendPoints
            int num = ins.size();
            for (int i = 0; i < num; i++) {
                double interpolation = num == 1 ? oneHalf : (1 + i) / (double) (num + 1);
                if (d == Direction.LEFT) {
                    double levelStartCoord = n.getProperty(InternalProperties.LEVELMAX);
                    // Only add level bendPoint if the distance is great enough
                    if (n.getPosition().x + n.getSize().x + edgeEndTexturePadding < levelStartCoord) {
                        ins.get(i).getBendPoints().add(levelStartCoord + nodeBendpointPadding, 
                                n.getPosition().y + n.getSize().y * interpolation);
                    // If the angle of the end piece is too steep, add another bendpoint
                    } else if (ins.get(i).getBendPoints().size() > 0) {
                        double lastX = ins.get(i).getBendPoints().getLast().x;
                        double nextX = n.getPosition().x + n.getSize().x / 2;
                        double lastY = ins.get(i).getBendPoints().getLast().y;
                        double nextY = n.getPosition().y + n.getSize().y / 2;
                        if (edgeEndTexturePadding > 0 && Math.abs(lastY - nextY) / (Math.abs(lastX - nextX) 
                                / steepEndEdgeSampleHeight) > steepEndEdgeTheresholdDistance) {
                            if (nextY > lastY) { 
                                // Place it to the left
                                ins.get(i).getBendPoints().add(
                                        n.getPosition().x + n.getSize().x + edgeEndTexturePadding / steepEndEdgeRatio, 
                                        n.getPosition().y + n.getSize().y * interpolation - edgeEndTexturePadding / 2);
                            } else {
                                // Place it to the right
                                ins.get(i).getBendPoints().add(
                                        n.getPosition().x + n.getSize().x + edgeEndTexturePadding / steepEndEdgeRatio, 
                                        n.getPosition().y + n.getSize().y * interpolation + edgeEndTexturePadding / 2);
                            }
                        }
                    }
                    ins.get(i).getBendPoints().add(n.getPosition().x + n.getSize().x, 
                            n.getPosition().y + n.getSize().y * interpolation);
                } else if (d == Direction.RIGHT) {
                    double levelStartCoord = n.getProperty(InternalProperties.LEVELMIN);
                    // Only add level bendPoint if the distance is great enough
                    if (n.getPosition().x - edgeEndTexturePadding > levelStartCoord) {
                        ins.get(i).getBendPoints().add(levelStartCoord - nodeBendpointPadding, 
                                n.getPosition().y + n.getSize().y * interpolation);
                    // If the angle of the end piece is too steep, add another bendpoint
                    } else if (ins.get(i).getBendPoints().size() > 0) {
                        double lastX = ins.get(i).getBendPoints().getLast().x;
                        double nextX = n.getPosition().x + n.getSize().x / 2;
                        double lastY = ins.get(i).getBendPoints().getLast().y;
                        double nextY = n.getPosition().y + n.getSize().y / 2;
                        if (edgeEndTexturePadding > 0 && Math.abs(lastY - nextY) / (Math.abs(lastX - nextX) 
                                / steepEndEdgeSampleHeight) > steepEndEdgeTheresholdDistance) {
                            if (nextY > lastY) {
                                ins.get(i).getBendPoints().add(
                                        n.getPosition().x - edgeEndTexturePadding / steepEndEdgeRatio, 
                                        n.getPosition().y + n.getSize().y * interpolation - edgeEndTexturePadding / 2);
                            } else {
                                ins.get(i).getBendPoints().add(
                                        n.getPosition().x - edgeEndTexturePadding / steepEndEdgeRatio, 
                                        n.getPosition().y + n.getSize().y * interpolation + edgeEndTexturePadding / 2);
                            }
                        }
                    }
                    ins.get(i).getBendPoints().add(n.getPosition().x, 
                            n.getPosition().y + n.getSize().y * interpolation);
                } else if (d == Direction.UP) {
                    double levelStartCoord = n.getProperty(InternalProperties.LEVELMAX);
                    // Only add level bendPoint if the distance is great enough
                    if (n.getPosition().y + n.getSize().y + edgeEndTexturePadding < levelStartCoord) {
                        ins.get(i).getBendPoints().add(n.getPosition().x + n.getSize().x * interpolation, 
                                levelStartCoord + nodeBendpointPadding);
                    // If the angle of the end piece is too steep, add another bendpoint
                    } else if (ins.get(i).getBendPoints().size() > 0) {
                        double lastX = ins.get(i).getBendPoints().getLast().x;
                        double nextX = n.getPosition().x + n.getSize().x / 2;
                        double lastY = ins.get(i).getBendPoints().getLast().y;
                        double nextY = n.getPosition().y + n.getSize().y / 2;
                        if (edgeEndTexturePadding > 0 && Math.abs(lastX - nextX) / (Math.abs(lastY - nextY) 
                                / steepEndEdgeSampleHeight) > steepEndEdgeTheresholdDistance) {
                            if (nextX > lastX) {
                                ins.get(i).getBendPoints().add(
                                        n.getPosition().x + n.getSize().x * interpolation - edgeEndTexturePadding / 2, 
                                        n.getPosition().y + edgeEndTexturePadding / steepEndEdgeRatio + n.getSize().y);
                            } else {
                                ins.get(i).getBendPoints().add(
                                        n.getPosition().x + n.getSize().x * interpolation + edgeEndTexturePadding / 2, 
                                        n.getPosition().y + edgeEndTexturePadding / steepEndEdgeRatio + n.getSize().y);
                            }
                        }
                    }
                    ins.get(i).getBendPoints().add(n.getPosition().x + n.getSize().x * interpolation, 
                            n.getPosition().y + n.getSize().y);
                } else {
                    double levelStartCoord = n.getProperty(InternalProperties.LEVELMIN);
                    // For cycle inducing edges add the end piece
                    if (TreeUtil.isCycleInducing(ins.get(i), tGraph)) {
                        ins.get(i).getBendPoints().add(n.getPosition().x + n.getSize().x * interpolation, 
                                ins.get(i).getBendPoints().getLast().y);
                    // Only add level bendPoint if the distance is great enough
                    } else if (n.getPosition().y - edgeEndTexturePadding > levelStartCoord) {
                        ins.get(i).getBendPoints().add(n.getPosition().x + n.getSize().x * interpolation, 
                                levelStartCoord - nodeBendpointPadding);
                    // If the angle of the end piece is too steep, add another bend point
                    } else if (ins.get(i).getBendPoints().size() > 0) {
                        double lastX = ins.get(i).getBendPoints().getLast().x;
                        double nextX = n.getPosition().x + n.getSize().x / 2;
                        double lastY = ins.get(i).getBendPoints().getLast().y;
                        double nextY = n.getPosition().y + n.getSize().y / 2;
                        if (edgeEndTexturePadding > 0 && Math.abs(lastX - nextX) / (Math.abs(lastY - nextY) 
                                / steepEndEdgeSampleHeight) > steepEndEdgeTheresholdDistance) {
                            if (nextX > lastX) {
                                ins.get(i).getBendPoints().add(
                                        n.getPosition().x + n.getSize().x * interpolation - edgeEndTexturePadding / 2, 
                                        n.getPosition().y - edgeEndTexturePadding / steepEndEdgeRatio);
                            } else {
                                ins.get(i).getBendPoints().add(
                                        n.getPosition().x + n.getSize().x * interpolation + edgeEndTexturePadding / 2, 
                                        n.getPosition().y - edgeEndTexturePadding / steepEndEdgeRatio);
                            }
                        }
                    }
                    
                    ins.get(i).getBendPoints().add(n.getPosition().x + n.getSize().x * interpolation, 
                            n.getPosition().y);
                }
            }
        }
    }

    /**
     * Handles CyrusBeck Edge Routing.
     * 
     * @param tGraph
     */
    private void cyrusBeck(final TGraph tGraph) {
        TNode root = TreeUtil.getRoot(tGraph);
        double nodeBendpointPadding = tGraph.getProperty(MrTreeOptions.SPACING_EDGE_NODE);
        
        for (TEdge edge : tGraph.getEdges()) {
            middleToMiddleEdgeRoute(edge);
        }
        
        for (TEdge edge : tGraph.getEdges()) {
            KVectorChain bends = edge.getBendPoints();
            for (int i = 0; i < bends.size() - 1; i++) {
                KVector start = bends.get(i);
                KVector end = bends.get(i + 1);
                KVector delta = start.clone().sub(end.clone());
                
                if (delta.x == 0 && delta.y == 0) {
                    continue;
                }
                
                for (TNode node : tGraph.getNodes()) {
                    
                    if (node.getLabel().equals(edge.getSource().getLabel()) 
                            || node.getLabel().equals(edge.getTarget().getLabel())) {
                        continue;
                    }
                    
                    ArrayList<KVector> clip = new ArrayList<KVector>();
                    double xmin = node.getPosition().x, xmax = node.getPosition().x + node.getSize().x, 
                            ymin = node.getPosition().y, ymax = node.getPosition().y + node.getSize().y;
                    
                    // Cohen-Sutherland Test
                    boolean so = start.y < ymin,
                            su = start.y > ymax,
                            sr = start.x > xmax,
                            sl = start.x < xmin,
                            eo = end.y < ymin,
                            eu = end.y > ymax,
                            er = end.x > xmax,
                            el = end.x < xmin;
                    
                    if (edge.getSource().getLabel().equals("n8") && node.getLabel().equals("n6")) {
                        getClass();
                    }
                    
                    // None of the points are inside the node but they may still overlap according to the test
                    if (!(so & eo) & !(su & eu) & !(sr & er) & !(sl & el)) {
                        
                        if (delta.x == 0) { // Special case: vert
                            clip.add(new KVector(start.x, ymin));
                            clip.add(new KVector(start.x, ymax));
                        } else if (delta.y == 0) { // Special case: horz
                            clip.add(new KVector(xmin, start.y));
                            clip.add(new KVector(xmax, start.y));
                        } else {
                            ArrayList<Triple<Double, KVector, Double>> pclips = 
                                    new ArrayList<Triple<Double, KVector, Double>>();
                            
                            if (so | eo) { // over
                                double t = (ymin - start.y) / (end.y - start.y);
                                KVector s = start.clone().add(delta.clone().scale(t));
                                pclips.add(new Triple<Double, KVector, Double>(t, s, -delta.y));
                            }
                            if (su | eu) { // under
                                double t = (ymax - start.y) / (end.y - start.y);
                                KVector s = start.clone().add(delta.clone().scale(t));
                                pclips.add(new Triple<Double, KVector, Double>(t, s, delta.y));
                            }
                            if (sl | el) { // left
                                double t = (xmin - start.x) / (end.x - start.x);
                                KVector s = start.clone().add(delta.clone().scale(t));
                                pclips.add(new Triple<Double, KVector, Double>(t, s, -delta.x));
                            }
                            if (sr | er) { // right
                                double t = (xmax - start.x) / (end.x - start.x);
                                KVector s = start.clone().add(delta.clone().scale(t));
                                pclips.add(new Triple<Double, KVector, Double>(t, s, delta.x));
                            }
                            
                            pclips.sort((x, y) -> Double.compare(x.getThird(), y.getThird()));
                            for (int j = 0; j < pclips.size(); j++) {
                                if (pclips.get(j).getThird() > 0) {
                                    clip.add(pclips.get(j).getSecond());
                                    clip.add(pclips.get(j - 1).getSecond());
                                    break;
                                }
                            }
                        }
                        
                        clip.sort((x, y) -> Double.compare(
                                delta.clone().dotProduct(x.clone()), 
                                delta.clone().dotProduct(y.clone())));
                        
                        KVector middle = clip.get(0).clone().add(clip.get(1).clone()).scale(1 / 2);
                        
                        bends.addAll(0, clip);
                        i += 2 * 2 * 2 * 2 * 2;
                        break;
                    }
                }
            }
        }
    }
}
