/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.alg.layered.intermediate.compaction;

import java.util.Collections;
import java.util.EnumSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Objects;
import java.util.Set;
import java.util.function.Function;
import java.util.stream.StreamSupport;

import org.eclipse.elk.alg.common.compaction.oned.CGraph;
import org.eclipse.elk.alg.common.compaction.oned.CGroup;
import org.eclipse.elk.alg.common.compaction.oned.CNode;
import org.eclipse.elk.alg.common.compaction.oned.CompareFuzzy;
import org.eclipse.elk.alg.common.compaction.oned.Quadruplet;
import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LNode.NodeType;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.layered.p5edges.splines.SplineSegment;
import org.eclipse.elk.core.math.ElkRectangle;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.options.Direction;
import org.eclipse.elk.core.options.EdgeRouting;
import org.eclipse.elk.core.options.PortSide;
import org.eclipse.elk.core.util.Pair;

import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;

/**
 * Manages the transformation of {@link LNode}s and {@link LEdge}s from an {@link LGraph} into
 * compactable {@link CNode}s and {@link CGroup}s.
 */
public final class LGraphToCGraphTransformer {
    
    /** the output CGraph. */
    private CGraph cGraph;
    /** the layered graph. */
    private LGraph layeredGraph;
    /** Current style of edge routing. */
    private EdgeRouting edgeRouting;
    /** remember comment boxes as we neglect them during compaction and offset them afterwards. */
    private Map<LNode, Pair<LNode, KVector>> commentOffsets = Maps.newHashMap();
    
    /* ------------------ Internal Mappings ------------------ */
    private Map<LNode, CNode> nodesMap = Maps.newHashMap();
    private Map<VerticalSegment, CNode> verticalSegmentsMap = Maps.newHashMap();
    private Map<CNode, Quadruplet> lockMap = Maps.newHashMap();
    
    // for debugging
    private static final Function<CNode, String> NODE_TO_STRING_DELEGATE =
            n -> ((LNode) n.origin).getProperty(InternalProperties.ORIGIN).toString();
    private static final Function<CNode, String> VS_TO_STRING_DELEGATE = n -> ((VerticalSegment) n.origin).toString();
    
    /**
     * @return the input graph transformed into a corresponding constraint graph.
     */
    public CGraph transform(final LGraph inputGraph) {
        this.layeredGraph = inputGraph;
        this.edgeRouting = layeredGraph.getProperty(LayeredOptions.EDGE_ROUTING);
        
        init();
        transformNodes();
        transformEdges();

        return cGraph;
    }

    /**
     * @return the lockMap
     */
    public Map<CNode, Quadruplet> getLockMap() {
        return lockMap;
    }
    
    private void init() {
        // checking if the graph has edges and possibly prohibiting vertical compaction
        boolean hasEdges = false;

        int index = 0;
        for (Layer l : layeredGraph) {
            l.id = index++;
            for (LNode n : l) {
                // avoid calling Iterables.isEmpty unnecessarily
                if (!hasEdges && !Iterables.isEmpty(n.getConnectedEdges())) {
                    hasEdges = true;
                }
            }
        }

        EnumSet<Direction> supportedDirections = EnumSet.of(Direction.UNDEFINED, Direction.LEFT, Direction.RIGHT);
        if (!hasEdges) {
            supportedDirections.add(Direction.UP);
            supportedDirections.add(Direction.DOWN);
        }

        // initializing members
        cGraph = new CGraph(supportedDirections);
        nodesMap.clear();
        commentOffsets.clear();
        lockMap.clear();
        verticalSegmentsMap.clear();
    }

    private void transformNodes() {
        for (Layer layer : layeredGraph) {
            for (LNode node : layer) {

                // comment boxes are part of a node's margins
                // hence we can neglect them here without the risk
                // of other nodes overlapping them after compaction
                if (node.getProperty(LayeredOptions.COMMENT_BOX)) {
                    if (!Iterables.isEmpty(node.getConnectedEdges())) {
                        LEdge e = Iterables.get(node.getConnectedEdges(), 0);
                        LNode other = e.getSource().getNode();
                        if (other == node) {
                            other = e.getTarget().getNode();
                        }
                        Pair<LNode, KVector> p = Pair.of(other, node.getPosition().clone().sub(other.getPosition()));
                        commentOffsets.put(node, p);
                        continue;
                    }
                }
                
                ElkRectangle hitbox =
                      new ElkRectangle(
                              node.getPosition().x - node.getMargin().left, 
                              node.getPosition().y - node.getMargin().top, 
                              node.getSize().x + node.getMargin().left + node.getMargin().right, 
                              node.getSize().y + node.getMargin().top + node.getMargin().bottom);

                // create the node in the compaction graph
                CNode cNode = CNode.of()
                                   .origin(node)
                                   .hitbox(hitbox)
                                   .toStringDelegate(NODE_TO_STRING_DELEGATE)
                                   .create(cGraph);
                // the node lives in its own group
                CGroup.of()
                      .nodes(cNode)
                      .master(cNode)
                      .create(cGraph);
                
                Quadruplet nodeLock = new Quadruplet();
                lockMap.put(cNode, nodeLock);

                // locking the node for directions that fewer edges are connected in
                // (only used if LEFT_RIGHT_CONNECTION_LOCKING is used)
                int difference = Iterables.size(node.getIncomingEdges())
                                 - Iterables.size(node.getOutgoingEdges());
                if (difference < 0) {
                    nodeLock.set(true, Direction.LEFT);
                } else if (difference > 0) {
                    nodeLock.set(true, Direction.RIGHT);
                }
                
                // excluding external port dummies
                if (node.getType() == NodeType.EXTERNAL_PORT) {
                    nodeLock.set(false, false, false, false);
                }
                
                nodesMap.put(node, cNode);
            }
        }
    }

    private void transformEdges() {
        EdgeRouting style = layeredGraph.getProperty(LayeredOptions.EDGE_ROUTING);
        List<VerticalSegment> verticalSegments;
        switch (style) {
            case ORTHOGONAL:
                verticalSegments = collectVerticalSegmentsOrthogonal();
                break;
            case SPLINES: 
                verticalSegments = collectVerticalSegmentsSplines();
                break;
            default:
                throw new IllegalStateException("Compaction not supported for " + style + " edges.");
        }

        // merge them
        mergeVerticalSegments(verticalSegments);
        
        // create precomputed constraints
        verticalSegmentsMap.keySet()
            .forEach(vs -> {
                CNode vsNode = verticalSegmentsMap.get(vs);
                vs.constraints.forEach(other -> {
                   CNode otherNode = verticalSegmentsMap.get(other);
                   cGraph.predefinedHorizontalConstraints.add(Pair.of(vsNode, otherNode));
                });
            });
    }
    
    private List<VerticalSegment> collectVerticalSegmentsOrthogonal() {
        List<VerticalSegment> verticalSegments = Lists.newArrayList();
        
        for (Layer layer : layeredGraph) {
            for (LNode node : layer) {
                final CNode cNode = nodesMap.get(node);
                
                // add vertical edge segments
                for (LEdge edge : node.getOutgoingEdges()) {

                    Iterator<KVector> bends = edge.getBendPoints().iterator();

                    boolean first = true;
                    VerticalSegment lastSegment = null;
                    // infer vertical segments from positions of bendpoints
                    if (bends.hasNext()) {
                        KVector bend1 = bends.next();
                        KVector bend2 = null;

                        // get segment of source n/s port of outgoing segment
                        if (edge.getSource().getSide() == PortSide.NORTH) {
                            VerticalSegment vs =
                                    new VerticalSegment(bend1, new KVector(bend1.x, cNode.hitbox.y), cNode, edge);
                            vs.ignoreSpacing.down = true;
                            vs.aPort = edge.getSource();
                            verticalSegments.add(vs);
                        }
                        if (edge.getSource().getSide() == PortSide.SOUTH) {
                            VerticalSegment vs = new VerticalSegment(bend1, new KVector(bend1.x,
                                    cNode.hitbox.y + cNode.hitbox.height), cNode, edge);
                            vs.ignoreSpacing.up = true;
                            vs.aPort = edge.getSource();
                            verticalSegments.add(vs);
                        }

                        // get regular segments
                        while (bends.hasNext()) {
                            bend2 = bends.next();
                            if (!CompareFuzzy.eq(bend1.y, bend2.y)) {
                                lastSegment = new VerticalSegment(bend1, bend2, null, edge);
                                verticalSegments.add(lastSegment);

                                // the first vertical segment of an outgoing edge
                                if (first) {
                                   first = false;
                                   
                                   // make sure that "invalid" spacing that sometimes is required 
                                   // to route the edges away from the nodes, does not become an issue
                                   if (bend2.y < cNode.hitbox.y) {
                                       lastSegment.ignoreSpacing.down = true;
                                   } else if (bend2.y > cNode.hitbox.y + cNode.hitbox.height) {
                                       lastSegment.ignoreSpacing.up = true;
                                   } else {
                                       // completely surrounded
                                       lastSegment.ignoreSpacing.up = true;
                                       lastSegment.ignoreSpacing.down = true;
                                   }
                                }
                            }

                            if (bends.hasNext()) {
                                bend1 = bend2;
                            }
                        }
                        
                        // handle last vertical segment
                        if (lastSegment != null) {
                            CNode cTargetNode = nodesMap.get(edge.getTarget().getNode());
                            if (bend1.y < cTargetNode.hitbox.y) {
                                lastSegment.ignoreSpacing.down = true;
                            } else if (bend1.y > cTargetNode.hitbox.y + cTargetNode.hitbox.height) {
                                lastSegment.ignoreSpacing.up = true;
                            } else {
                                // completely surrounded
                                lastSegment.ignoreSpacing.up = true;
                                lastSegment.ignoreSpacing.down = true;
                            }
                        }
                    }
                }

                // same for incoming edges to get north/south segments on target side
                for (LEdge edge : node.getIncomingEdges()) {
                    if (!edge.getBendPoints().isEmpty()) {

                        // get segment of target n/s port
                        KVector bend1 = edge.getBendPoints().getLast();
                        if (edge.getTarget().getSide() == PortSide.NORTH) {
                            VerticalSegment vs = new VerticalSegment(bend1, new KVector(bend1.x,
                                    cNode.hitbox.y), cNode, edge);
                            vs.ignoreSpacing.down = true;
                            vs.aPort = edge.getTarget();
                            verticalSegments.add(vs);
                        }
                        if (edge.getTarget().getSide() == PortSide.SOUTH) {
                            VerticalSegment vs = new VerticalSegment(bend1, new KVector(bend1.x,
                                    cNode.hitbox.y + cNode.hitbox.height), cNode, edge);
                            vs.ignoreSpacing.up = true;
                            vs.aPort = edge.getTarget();
                            verticalSegments.add(vs);
                        }
                        
                    }
                }
            }
        }
        
        return verticalSegments;
    }
    
    private List<VerticalSegment> collectVerticalSegmentsSplines() {
        List<VerticalSegment> verticalSegments = Lists.newArrayList();
        
        layeredGraph.getLayers().stream()
            .flatMap(l -> l.getNodes().stream())
            .flatMap(n -> StreamSupport.stream(n.getOutgoingEdges().spliterator(), false))
            .map(out -> out.getProperty(InternalProperties.SPLINE_ROUTE_START))
            .filter(Objects::nonNull) 
            .forEach(spline -> {
                VerticalSegment lastVs = null;
                for (SplineSegment s : spline) {
                    if (s.isStraight) {
                        continue;
                    }
                    KVector leftTop = s.boundingBox.getTopLeft();
                    KVector rightBottom = s.boundingBox.getBottomRight();

                    VerticalSegment vs = new VerticalSegment(leftTop, rightBottom, null, s.edges.iterator().next());
                    vs.affectedBoundingBoxes.add(s.boundingBox);

                    verticalSegments.add(vs);
                    
                    // remember that there has to be a constraint between the two (non-straight) segments
                    if (lastVs != null) {
                        lastVs.constraints.add(vs);
                    } 
                    lastVs = vs;
                }
            });
        
        return verticalSegments;
    }
    
    private void mergeVerticalSegments(final List<VerticalSegment> verticalSegments) {
        if (verticalSegments.isEmpty()) {
            return;
        }
        
        // sort the segments such that segments with the same 
        //  x coordinate are consecutive within the list, 
        //  for equal x coordinates, sort with ascending y coordinate 
        Collections.sort(verticalSegments);

        // merge segments that intersect
        Iterator<VerticalSegment> vsIt = verticalSegments.iterator();
        VerticalSegment survivor = vsIt.next();

        while (vsIt.hasNext()) {
            VerticalSegment next = vsIt.next();
            
            if (survivor.intersects(next)) {
                survivor = survivor.joinWith(next);
            } else {
                verticalSegmentToCNode(survivor);
                // start with a new segment
                survivor = next;
            }
        }
        
        // transform the last segment as well
        verticalSegmentToCNode(survivor);
    }
    
    private void verticalSegmentToCNode(final VerticalSegment verticalSegment) {
        // create CNode representation for the last segment
        CNode cNode = CNode.of()
                           .origin(verticalSegment)
                           .type("vs")
                           .hitbox(new ElkRectangle(verticalSegment.hitbox))
                           .toStringDelegate(VS_TO_STRING_DELEGATE)
                           .create(cGraph);
        
        // group the node if requested, currently this is only the case for north/south segments of orthogonal edges 
        // assumption: during creation of LNode representing CNodes, these CNodes have been fitted with their own group
        if (!verticalSegment.potentialGroupParents.isEmpty()) {
            verticalSegment.potentialGroupParents.get(0).cGroup.addCNode(cNode);
        }
        
        Quadruplet vsLock = new Quadruplet();
        lockMap.put(cNode, vsLock);

        // segments belonging to multiple edges should be locked 
        // in the direction that fewer different ports are connected in
        // (only used if LEFT_RIGHT_CONNECTION_LOCKING is active)
        Set<LPort> inc = Sets.newHashSet();
        Set<LPort> out = Sets.newHashSet();
        for (LEdge e : verticalSegment.representedLEdges) {
            inc.add(e.getSource());
            out.add(e.getTarget());
        }
        int difference = inc.size() - out.size();
        if (difference < 0) {
            vsLock.set(true, Direction.LEFT);
            vsLock.set(false, Direction.RIGHT);
        } else if (difference > 0) {
            vsLock.set(false, Direction.LEFT);
            vsLock.set(true, Direction.RIGHT);
        }
        
        verticalSegment.joined.forEach(other -> verticalSegmentsMap.put(other, cNode));
        verticalSegmentsMap.put(verticalSegment, cNode);
    }

    /**
     * Apply the layout from the internal constraint graph back to the original lgraph.
     */
    public void applyLayout() {
        
        // apply the compacted positions to nodes
        cGraph.cNodes.stream()
            .filter(cNode -> cNode.origin instanceof LNode)
            .forEach(cNode -> {
                LNode lNode = (LNode) cNode.origin;
                lNode.getPosition().x = cNode.hitbox.x + lNode.getMargin().left;
            });
        
        // adjust comment boxes
        applyCommentPositions();

        // apply new positions to vertical segments
        cGraph.cNodes.stream()
            .filter(cNode -> cNode.origin instanceof VerticalSegment)
            .forEach(cNode -> {
                double deltaX = cNode.hitbox.x - cNode.hitboxPreCompaction.x;
                VerticalSegment vs = (VerticalSegment) cNode.origin;
                vs.affectedBends.forEach(b -> b.x += deltaX);
                vs.affectedBoundingBoxes.forEach(bb -> bb.x += deltaX);
                vs.junctionPoints.forEach(jp -> jp.x += deltaX);
            });
        
        // special treatment of spline edge routes
        //  1) self loops have already been routed and are not part of the compaction graph, 
        //     their positions have to be adjusted
        //  2) control points of straight parts of a spline have to be offset to 
        //     be properly located in between non-straight parts
        if (edgeRouting == EdgeRouting.SPLINES) {
            
            // offset selfloops of splines
            nodesMap.keySet().stream()
                .flatMap(n -> StreamSupport.stream(n.getOutgoingEdges().spliterator(), false))
                .filter(e -> e.isSelfLoop())
                .forEach(sl -> {
                   LNode lNode = sl.getSource().getNode();
                   CNode cNode = nodesMap.get(lNode);
                   double deltaX = cNode.hitbox.x - cNode.hitboxPreCompaction.x;
                   sl.getBendPoints().offset(deltaX, 0);
                });
            
            // offset straight segments
            layeredGraph.getLayers().stream()
                .flatMap(l -> l.getNodes().stream())
                .flatMap(n -> StreamSupport.stream(n.getOutgoingEdges().spliterator(), false))
                .map(e -> e.getProperty(InternalProperties.SPLINE_ROUTE_START))
                .filter(chain -> chain != null && !chain.isEmpty())
                .forEach(spline -> adjustSplineControlPoints(spline));
            
        }
        
        // offset selfloop labels 
        nodesMap.keySet().stream()
        .flatMap(n -> StreamSupport.stream(n.getOutgoingEdges().spliterator(), false))
        .filter(e -> e.isSelfLoop())
        .forEach(sl -> {
           LNode lNode = sl.getSource().getNode();
           CNode cNode = nodesMap.get(lNode);
           double deltaX = cNode.hitbox.x - cNode.hitboxPreCompaction.x;
           sl.getLabels().forEach(l -> {
              l.getPosition().x += deltaX;
           });
        });
        
        // calculating new graph size and offset
        KVector topLeft = new KVector(Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY);
        KVector bottomRight = new KVector(Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY);
        for (CNode cNode : cGraph.cNodes) {
            topLeft.x = Math.min(topLeft.x, cNode.hitbox.x);
            topLeft.y = Math.min(topLeft.y, cNode.hitbox.y);
            bottomRight.x = Math.max(bottomRight.x, cNode.hitbox.x + cNode.hitbox.width);
            bottomRight.y = Math.max(bottomRight.y, cNode.hitbox.y + cNode.hitbox.height);
        }
        layeredGraph.getOffset().reset().add(topLeft.clone().negate());
        layeredGraph.getSize().reset().add(bottomRight.clone().sub(topLeft));
        
        
        // external port dummies may have been moved ... put them back
        applyExternalPortPositions(topLeft, bottomRight);
        
        cleanup();
    }
    
    /**
     * Adjust the control points of the straight segments of the passed spline. The control points of non-straight
     * segments are directly represented by {@link CNode}s during compaction and thus adjusted during the compaction
     * process.
     */
    private void adjustSplineControlPoints(final List<SplineSegment> spline) {
        if (spline.isEmpty()) {
            return;
        }
        
        SplineSegment lastSeg = spline.get(0);

        // first case: a signle segment
        if (spline.size() == 1) {
            // the two indexes are selected such that the method's loop runs exactly once
            adjustControlPointBetweenSegments(lastSeg, lastSeg, 1, 0, spline);
            return;
        }

        // ... more than one segment
        int i = 1;
        while (i < spline.size()) {
            if (lastSeg.initialSegment || !lastSeg.isStraight) {
                
                // find the next non-straight segment (or the very last segment of this spline)
                Pair<Integer, SplineSegment> needle = firstNonStraightSegment(spline, i);
                if (needle != null) {
                    int j = needle.getFirst();
                    SplineSegment nextSeg = needle.getSecond();
                    
                    adjustControlPointBetweenSegments(lastSeg, nextSeg, i, j, spline);
                    
                    i = j + 1; 
                    lastSeg = nextSeg;
                }
            }
        }
    }

    /**
     * Finds the next non-straight segment of {@code spline}, starting at {@code index} (inclusive), or the very last
     * segment of the spline.
     * 
     * @return a {@link Pair} of the found segment's index and the segment. {@code null} if
     *         {@code index < 0 || index >= spline.size()}.
     */
    private Pair<Integer, SplineSegment> firstNonStraightSegment(final List<SplineSegment> spline, final int index) {
        if (index < 0 || index >= spline.size()) {
            return null;
        }
        for (int i = index; i < spline.size(); ++i) {
            SplineSegment seg = spline.get(i);
            if ((i == spline.size() - 1) || !seg.isStraight) {
                return Pair.of(i, seg);
            }
        }
        return null;
    }
    
    /**
     * Adjust the control points of the sequence of spline segments {@code left, ..., right}. The {@code left} and
     * {@code right} segments do not have to be straight in which case they are dropped from the sequence and their
     * positions remain unaltered. For the remaining sequence of length {@code n}, the left-most and right-most
     * coordinate of the straight path are determined. The center control points of the sequence's segments are then
     * distributed equidistantly on this path.
     */
    private void adjustControlPointBetweenSegments(final SplineSegment left, final SplineSegment right, 
            final int leftIdx, final int rightIdx,
            final List<SplineSegment> spline) {
        
        // check if the initial segment of the spline is a straight one
        double startX;
        int idx1 = leftIdx;
        if (left.initialSegment && left.isStraight) {
            CNode n = nodesMap.get(left.sourceNode);
            startX = n.hitbox.x + n.hitbox.width;
            idx1--;
        } else {
            startX = left.boundingBox.x + left.boundingBox.width;
        }
        
        // ... the same for the last segment
        double endX;
        int idx2 = rightIdx;
        if (right.lastSegment && right.isStraight) {
            CNode n = nodesMap.get(right.targetNode);
            endX = n.hitbox.x;
            idx2++;
        } else {
            endX = right.boundingBox.x;
        }
        
        // divide the available space into equidistant chunks 
        double strip = endX - startX;
        int chunks = Math.max(2, idx2 - idx1);
        double chunk = strip / chunks;

        // apply new positions to the control points
        double newPos = startX + chunk;
        for (int k = idx1; k < idx2; ++k) {
            SplineSegment adjust = spline.get(k);
            double width = adjust.boundingBox.width;
            adjust.boundingBox.x = newPos - width / 2;
            
            newPos += chunk;
        }
    }

    private void applyCommentPositions() {
        for (Entry<LNode, Pair<LNode, KVector>> e : commentOffsets.entrySet()) {
            LNode comment = e.getKey();
            LNode other = e.getValue().getFirst();
            KVector offset = e.getValue().getSecond();
            comment.getPosition().reset().add(other.getPosition().clone().add(offset));
        }
    }
    
    private void applyExternalPortPositions(final KVector topLeft, final KVector bottomRight) {
        for (CNode cNode : cGraph.cNodes) {
            LNode lNode = HorizontalGraphCompactor.getLNodeOrNull(cNode);
            if (lNode != null) {
                if (lNode.getType() == NodeType.EXTERNAL_PORT) {
                    switch (lNode.getProperty(InternalProperties.EXT_PORT_SIDE)) {
                    case WEST:
                        lNode.getPosition().x = topLeft.x;
                        break;
                    case EAST:
                        lNode.getPosition().x = bottomRight.x
                                - (lNode.getSize().x + lNode.getMargin().right); 
                        break;
                    case NORTH: 
                        lNode.getPosition().y = topLeft.y;
                        break;
                    case SOUTH:
                        lNode.getPosition().y = bottomRight.y
                                - (lNode.getSize().y + lNode.getMargin().bottom);
                        break;
                    }
                }
            }
        }
    }
    
    private void cleanup() {
        nodesMap.clear();
        commentOffsets.clear();
        verticalSegmentsMap.clear();
        lockMap.clear();
        cGraph.cGroups.clear();
        cGraph.cNodes.clear();
        cGraph = null;
        layeredGraph = null;
    }

}
