/*******************************************************************************
 * Copyright (c) 2016 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.layered.components;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import org.eclipse.elk.alg.layered.compaction.components.IComponent;
import org.eclipse.elk.alg.layered.compaction.components.IConnectedComponents;
import org.eclipse.elk.alg.layered.compaction.components.IExternalExtension;
import org.eclipse.elk.alg.layered.compaction.components.OneDimensionalComponentsCompaction;
import org.eclipse.elk.alg.layered.compaction.recthull.Point;
import org.eclipse.elk.alg.layered.compaction.recthull.RectilinearConvexHull;
import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LGraphUtil;
import org.eclipse.elk.alg.layered.graph.LLabel;
import org.eclipse.elk.alg.layered.graph.LMargin;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LNode.NodeType;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.graph.LShape;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.core.math.ElkRectangle;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.math.KVectorChain;
import org.eclipse.elk.core.options.Direction;
import org.eclipse.elk.core.options.PortSide;
import org.eclipse.elk.core.util.BasicProgressMonitor;
import org.eclipse.elk.core.util.Pair;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;
import com.google.common.collect.Sets;
import com.google.common.math.DoubleMath;

/**
 * Contains implementations of the {@link IConnectedComponents} and {@link IComponent} interfaces
 * for {@link LGraph}s. As a consequence allows to compact a set of components using one dimensional
 * compaction.
 * 
 * @see OneDimensionalComponentsCompaction
 * @see org.eclipse.elk.alg.layered.compaction.components.ComponentsToCGraphTransformer ComponentsToCGraphTransformer
 */
public class ComponentsCompactor {

    /** Instance of the compactor we use. */
    private OneDimensionalComponentsCompaction<LNode, Set<LEdge>> compactor;
    
    /** The offset after compaction has been applied. */
    private KVector yetAnotherOffset;
    /** The graph size after compaciton. */
    private KVector compactedGraphSize;
    
    private KVector graphTopLeft;
    private KVector graphBottomRight;
    
    /** Epsilon for double comparisons. */
    private static final double EPSILON = 0.0001;

    // ------------------------------------------------------------------------------------------------
    //                                          public API
    // ------------------------------------------------------------------------------------------------
    
    /**
     * @param graphs
     *            the components to be compacted
     * @param originalGraphsSize
     *            the size of the overall graph as it is currently
     * @param spacing
     *            the desired spacing to be preserved between any pair of components
     */
    public void compact(final List<LGraph> graphs, final KVector originalGraphsSize,
            final double spacing) {

        // determine the extreme points of the current diagram, 
        // we will reuse this 'frame' to cut external extensions at appropriate lengths
        graphTopLeft = new KVector(Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY);
        graphBottomRight = new KVector(Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY);
        for (LGraph graph : graphs) {
            for (LNode node : graph.getLayerlessNodes()) {
                graphTopLeft.x = Math.min(graphTopLeft.x, node.getPosition().x - node.getMargin().left);
                graphTopLeft.y = Math.min(graphTopLeft.y, node.getPosition().y - node.getMargin().top);
                graphBottomRight.x = Math.max(graphBottomRight.x,
                                node.getPosition().x + node.getSize().x + node.getMargin().right);
                graphBottomRight.y = Math.max(graphBottomRight.y,
                                node.getPosition().y + node.getSize().y + node.getMargin().bottom);
            }
        }
        
        // from the lgraphs, create connected components
        IConnectedComponents<LNode, Set<LEdge>> ccs = new InternalConnectedComponents();
        for (LGraph graph : graphs) {
            IComponent<LNode, Set<LEdge>> c = transformLGraph(graph);
            ccs.getComponents().add(c);
            ((InternalComponent) c).containsRegularNodes |=
                    !c.getExternalExtensionSides().isEmpty();
        }
        
        // for every component we create an element in the compactor
        compactor = OneDimensionalComponentsCompaction.init(ccs, spacing);
        
        // execute compaction
        compactor.compact(new BasicProgressMonitor());

        yetAnotherOffset = new KVector();
        compactedGraphSize = compactor.getGraphSize();
        
        // apply the positions
        for (IComponent<LNode, Set<LEdge>> cc : ccs.getComponents()) {
            
            // retrieve the common offset for the currently handled connected component
            KVector offset = compactor.getOffset(cc);

            // move it
            LGraphUtil.offsetGraph(((InternalComponent) cc).graph, offset.x, offset.y);

            // adjust positions of external ports
            for (LNode n : ((InternalComponent) cc).getNodes()) {
  
                if (n.getType() == NodeType.EXTERNAL_PORT) {

                    KVector newPos = getExternalPortPosition(n.getPosition(), 
                            n.getProperty(InternalProperties.EXT_PORT_SIDE));
                    
                    n.getPosition().reset().add(newPos);
                }
            }
        }
        
        // external edges contribute to the graph's size ... however, only certain segments do.
        for (IComponent<LNode, Set<LEdge>> cc : ccs.getComponents()) {
            for (LEdge e : ((InternalComponent) cc).getExternalEdges()) {
                KVectorChain vc = new KVectorChain(e.getBendPoints());
                vc.add(0, e.getSource().getAbsoluteAnchor());
                vc.add(e.getTarget().getAbsoluteAnchor());
                
                KVector last = null;
                for (KVector v : vc) {
                    if (last == null) {
                        last = v;
                        continue;
                    }
                    if (DoubleMath.fuzzyEquals(last.x, v.x, EPSILON)) {
                        yetAnotherOffset.x = Math.min(yetAnotherOffset.x, last.x);
                        compactedGraphSize.x = Math.max(compactedGraphSize.x, last.x);
                    } else if (DoubleMath.fuzzyEquals(last.y, v.y, EPSILON)) {
                        yetAnotherOffset.y = Math.min(yetAnotherOffset.y, last.y);
                        compactedGraphSize.y = Math.max(compactedGraphSize.y, last.y);
                    }
                    last = v;
                }
            }
        }
        
        yetAnotherOffset.negate();
        compactedGraphSize.add(yetAnotherOffset);
    }
    
    /**
     * @return the offset by which each component has to be shifted after compaction such that the
     *         top-left-most point is (0, 0).
     */
    public KVector getOffset() {
        return yetAnotherOffset;
    }
    
    /**
     * @return the new graph size.
     */
    public KVector getGraphSize() {
        return compactedGraphSize;
    }
    
    // ------------------------------------------------------------------------------------------------
    //                                          private API
    // ------------------------------------------------------------------------------------------------

    
    /**
     * Converts a {@link LGraph} into an component, i.e. a set of rectangles 
     * describing the space occupied by the component and a set of external extensions.
     */
    private IComponent<LNode, Set<LEdge>> transformLGraph(final LGraph graph) {
        
        InternalComponent component = new InternalComponent(graph);
        
        if (!component.containsRegularNodes) {
            createDummyNode(graph);
        }
        
        // #1 convert the nodes to hull points
        Hullpoints hullPoints = componentHullPoints(graph);
       
        // #2 convert external edges,
        //    while doing this, add further hull points contributed by inner segments
        //    and remember extreme points of external segments
        Multimap<Direction, LEdge> externalExtensions = HashMultimap.create();
        OuterSegments outerSegments = new OuterSegments();
        
        for (LNode node : graph.getLayerlessNodes()) {
            for (LEdge edge : node.getOutgoingEdges()) {
                if (isExternalEdge(edge)) {
                    IExternalExtension<LEdge> iee = transformLEdge(edge, hullPoints, outerSegments);
                    externalExtensions.put(iee.getDirection(), iee.getRepresentative());
                }
            }
        }

        // #3 create the common external extensions for this component
        //    (there can be 4 per component, one in each direction)
        List<InternalUnionExternalExtension> extensions = Lists.newArrayList();
        for (PortSide ps : component.getExternalExtensionSides()) {
            double min = outerSegments.min[ps.ordinal()];
            double max = outerSegments.max[ps.ordinal()];
            double extent = outerSegments.extent[ps.ordinal()];
            
            ElkRectangle extension = null;
            ElkRectangle placeholder = null;
            switch (ps) {
            case WEST:
                extension = new ElkRectangle(graphTopLeft.x,
                                          min,
                                          hullPoints.topLeft.x - graphTopLeft.x,
                                          max - min);
                placeholder = new ElkRectangle(graphTopLeft.x,
                                            min,
                                            extent,
                                            max - min);
                hullPoints.add(extension.getTopRight());
                hullPoints.add(extension.getBottomRight());
                break;
                
            case EAST: 
                extension = new ElkRectangle(hullPoints.bottomRight.x,
                                          min,
                                          graphBottomRight.x - hullPoints.bottomRight.x,
                                          max - min);
                placeholder = new ElkRectangle(graphBottomRight.x - extent,
                                            min,
                                            extent,
                                            max - min);
                hullPoints.add(extension.getTopLeft());
                hullPoints.add(extension.getBottomLeft());
                break;
                
            case NORTH:
                extension = new ElkRectangle(min, 
                                          graphTopLeft.y,
                                          max - min,
                                          hullPoints.topLeft.y - graphTopLeft.y);
                placeholder = new ElkRectangle(min, 
                                            graphTopLeft.y,
                                            max - min,
                                            extent);
                hullPoints.add(extension.getBottomLeft());
                hullPoints.add(extension.getBottomRight());
                break;
                
            case SOUTH: 
                extension = new ElkRectangle(min, 
                                          hullPoints.bottomRight.y, 
                                          max - min, 
                                          graphBottomRight.y - hullPoints.bottomRight.y);
                placeholder = new ElkRectangle(min, 
                                            graphBottomRight.y - extent, 
                                            max - min, 
                                            extent);
                hullPoints.add(extension.getTopLeft());
                hullPoints.add(extension.getTopRight());
                break;
            }
            // instantiate the external extension object
            if (extension != null) {
                InternalUnionExternalExtension iuee = new InternalUnionExternalExtension();
                iuee.side = ps;
                iuee.extension = extension;
                iuee.placeholder = placeholder;
                iuee.edges = Sets.newHashSet(externalExtensions.get(portSideToDirection(ps)));
                extensions.add(iuee);
            }
            
        }
        
        // #4 calculate the hull for the component
        component.externalExtensions.addAll(extensions);
        component.rectilinearConvexHull = RectilinearConvexHull.of(hullPoints)
                                                               .splitIntoRectangles();
        
        return component;
    }
    
    private LNode createDummyNode(final LGraph graph) {

        // this can only be the case if there's a single external port
        //  to handle it properly, we give it a tiny little dummy rectangle
        assert graph.getLayerlessNodes().size() == 1;
        
        LNode extPortDummy = graph.getLayerlessNodes().get(0);

        // create a small dummy node for the external port dummy
        LNode dummy = new LNode(graph);
        graph.getLayerlessNodes().add(dummy);
        // reassemble the size of the external dummy but assure 
        //  that no dimension is zero. Otherwise we may end up 
        //  with dimension-less hull rectangles which causes
        //  the scanline to fail
        dummy.getSize().x = Math.max(1, extPortDummy.getSize().x);
        dummy.getSize().y = Math.max(1, extPortDummy.getSize().y);
        
        // set the position such that it is on the proper side 
        // with respect to the dummy port
        dummy.getPosition().x = extPortDummy.getPosition().x;
        dummy.getPosition().y = extPortDummy.getPosition().y;
        switch (extPortDummy.getProperty(InternalProperties.EXT_PORT_SIDE)) {
            case WEST:
                dummy.getPosition().x += 2;
                break;
            case NORTH:
                dummy.getPosition().y += 2;
                break;
            case EAST:
                dummy.getPosition().x -= 2;
                break;
            case SOUTH:
                dummy.getPosition().y -= 2;
                break;
        }
        
        // give the dummy a port and create an edge for it
        LPort dummyPort = new LPort();
        dummyPort.setNode(dummy);
        LEdge dummyEdge = new LEdge();
        LPort extPortDummyPort = extPortDummy.getPorts().get(0);
        dummyEdge.setSource(extPortDummyPort);
        dummyEdge.setTarget(dummyPort);
        // position the dummy port such that the edge will be straight
        dummyPort.getPosition().reset().add(extPortDummyPort.getPosition());
        dummyPort.getAnchor().reset().add(extPortDummyPort.getAnchor());
        
        return dummy;
    }
    
    private Hullpoints componentHullPoints(final LGraph graph) {
        
        final Hullpoints pts = new Hullpoints();
        
        for (LNode n : graph.getLayerlessNodes()) {
           
            if (n.getType() == NodeType.EXTERNAL_PORT) {
                continue;   
            }

            // Note that labels of nodes are already part of a node's margins
            // the same is true for ports and their labels
            addLGraphElementBounds(pts, n, new KVector());

            // add bend points of the edges
            for (LEdge edge : n.getOutgoingEdges()) {
                if (isExternalEdge(edge)) {
                    continue;
                }
                
                for (KVector bp : edge.getBendPoints()) {
                    KVector absolute = bp;
                    pts.add(new Point(absolute.x, absolute.y));   
                }
            }
        }
        
        return pts;
    }
    
    private void addLGraphElementBounds(final List<Point> pts, final LShape element,
            final KVector offset) {
        
        // extract the relevant margins object.
        //  there's LayoutOptions.MARGINS as well,
        //  however, this is only used outside of elk layered.
        LMargin margins = null;
        if (element instanceof LNode) {
            margins = ((LNode) element).getMargin();
        } else if (element instanceof LPort) {
            margins = ((LPort) element).getMargin();
        } else if (element instanceof LLabel) {
            margins = new LMargin();
        }
        // add bounding box of the node
        pts.add(new Point(element.getPosition().x - margins.left + offset.x, 
                element.getPosition().y - margins.top + offset.y));
        pts.add(new Point(element.getPosition().x - margins.left + offset.x, 
                element.getPosition().y + element.getSize().y + margins.bottom + offset.y));
        pts.add(new Point(
                element.getPosition().x + element.getSize().x + margins.right + offset.x, 
                element.getPosition().y - margins.top + offset.y));
        pts.add(new Point(
                element.getPosition().x + element.getSize().x + margins.right + offset.x, 
                element.getPosition().y + element.getSize().y + margins.bottom + offset.y));
    }
    
    private boolean isExternalEdge(final LEdge edge) {
        return edge.getSource().getNode().getType() == NodeType.EXTERNAL_PORT
                || edge.getTarget().getNode().getType() == NodeType.EXTERNAL_PORT;
    }
    
    /**
     * Converts an external edge to an external extension. While doing so, the the original edge
     * contributes three things as explained in the code.
     */
    private IExternalExtension<LEdge> transformLEdge(final LEdge externalEdge,
            final Hullpoints hullPoints, 
            final OuterSegments outerSegments) {
        
        InternalExternalExtension externalExtension = new InternalExternalExtension(externalEdge);
        
        // #1 convert the edge's path into a set of segments
        Segments segments = edgeToSegments(externalEdge, externalExtension);
        
        // #2 all 'inner' segments contribute to the hull (consider the edge's thickness)
        double thickness = Math.max(externalEdge.getProperty(LayeredOptions.EDGE_THICKNESS).doubleValue(), 1);
        for (Pair<KVector, KVector> segment : segments.innerSegments) {
            ElkRectangle rect = segmentToRectangle(segment.getFirst(), segment.getSecond(), thickness);
            hullPoints.add(rect);
        }
        
        // #3 the 'outer' segment, being the segment that actually connects to the external port,
        //    contributes to the 'union external segment' that we create 
        //    for one direction of the component
        PortSide side = externalExtension.externalPortSide;
        ElkRectangle outerSegmentRect =
                segmentToRectangle(segments.outerSegment.getFirst(),
                        segments.outerSegment.getSecond(), thickness);
        if (side == PortSide.WEST || side == PortSide.EAST) {
            outerSegments.min[side.ordinal()] =
              Math.min(outerSegments.min[side.ordinal()], outerSegmentRect.y);
            outerSegments.max[side.ordinal()] =
              Math.max(outerSegments.max[side.ordinal()], outerSegmentRect.y + outerSegmentRect.height);
        } else {
            outerSegments.min[side.ordinal()] =
              Math.min(outerSegments.min[side.ordinal()], outerSegmentRect.x);
            outerSegments.max[side.ordinal()] =
              Math.max(outerSegments.max[side.ordinal()], outerSegmentRect.x + outerSegmentRect.width);
        }
        
        // extent
        double extent = Double.NEGATIVE_INFINITY;
        LMargin margins = externalExtension.externalPort.getNode().getMargin();
        switch (side) {
            case WEST:
                extent = margins.right;
                break;
            case EAST:
                extent = margins.left;
                break;
            case NORTH:
                extent = margins.bottom;
                break;
            case SOUTH:
                extent = margins.top;
                break;
        }
        outerSegments.extent[side.ordinal()] = Math.max(outerSegments.extent[side.ordinal()], extent);
        
        return externalExtension;
    }
    
    // ------------------------------------------------------------------------------------------------
    //                                   Auxiliary stuff
    // ------------------------------------------------------------------------------------------------
    
    private ElkRectangle segmentToRectangle(final KVector p1, final KVector p2, final double extent) {
        return new ElkRectangle(Math.min(p1.x, p2.x) - extent / 2f,
                             Math.min(p1.y, p2.y) - extent / 2f,
                             Math.abs(p1.x - p2.x) + extent,
                             Math.abs(p1.y - p2.y) + extent);
    }
    
    private Segments edgeToSegments(final LEdge edge,
            final InternalExternalExtension externalExtension) {

        LPort externalPort = externalExtension.externalPort;
        PortSide externalPortSide = externalExtension.externalPortSide;
        
        // extract the correct segment that is to be used for the external extension
        KVector p1 = edge.getSource().getAbsoluteAnchor();
        KVector p2 = edge.getTarget().getAbsoluteAnchor();
        
        // external ports of components may not be routed until the very border
        // of the diagram yet, hence we adjust their position here
        if (externalPort == edge.getSource()) {
            p1 = getExternalPortPosition(p1, externalPortSide);
            p2 = getPortPositionOnMargin(edge.getTarget());
        } else {
            p1 = getPortPositionOnMargin(edge.getSource());
            p2 = getExternalPortPosition(p2, externalPortSide);
        }

        KVectorChain points = new KVectorChain(edge.getBendPoints());
        points.addFirst(p1);
        points.addLast(p2);
        
        boolean outerSegmentIsFirst = edge.getSource() == externalPort;
        
        // for easy processing, make it a list of segments
        Segments segments = new Segments();
        for (int i = 0; i < points.size() - 1; ++i) {
            Pair<KVector, KVector> segment = Pair.of(points.get(i), points.get(i + 1));

            if ((outerSegmentIsFirst && i == 0) || (!outerSegmentIsFirst && i == points.size() - 2)) {
                segments.outerSegment = segment;
            } else {
                segments.innerSegments.add(segment);
            }
        }
        
        return segments;
    }
    
    private KVector getExternalPortPosition(final KVector pos, final PortSide ps) {
        switch (ps) {
        case NORTH:
            return new KVector(pos.x, Math.min(graphTopLeft.y, pos.y));
        case EAST:
            return new KVector(Math.max(graphBottomRight.x, pos.x), pos.y);
        case SOUTH:
            return new KVector(pos.x, Math.max(graphBottomRight.y, pos.y));
        case WEST:
            return new KVector(Math.min(pos.x, graphTopLeft.x), pos.y);
        }
        return pos.clone();
    }
    
    private KVector getPortPositionOnMargin(final LPort port) {
        KVector pos = port.getAbsoluteAnchor().clone();
        LMargin margins = port.getNode().getMargin();
        
        switch (port.getSide()) {
        case NORTH:
            pos.y -= margins.top;
            break;
        case EAST:
            pos.x += margins.right;
            break;
        case SOUTH:
            pos.y += margins.bottom;
            break;
        case WEST:
            pos.x -= margins.left;
            break;
        }
        return pos;
    }
    
    private Direction portSideToDirection(final PortSide side) {
        switch (side) {
        case NORTH:
            return Direction.UP;
        case WEST:
            return Direction.LEFT;
        case EAST:
            return Direction.RIGHT;
        case SOUTH:
            return Direction.DOWN;
        }
        return Direction.UNDEFINED;
    }
    
    // ------------------------------------------------------------------------------------------------
    //                                     Interface Implementations
    // ------------------------------------------------------------------------------------------------

    /**
     * Implementation of the {@link IConnectedComponents} interface. Holds the components 
     * and whether the overall graph contains external edges
     */
    private class InternalConnectedComponents implements IConnectedComponents<LNode, Set<LEdge>> {

        private List<IComponent<LNode, Set<LEdge>>> components = Lists.newArrayList();
        private boolean containsExternalPorts = false;
        
        InternalConnectedComponents() { }
        
        @Override
        public Iterator<IComponent<LNode, Set<LEdge>>> iterator() {
            return components.iterator();
        }

        @Override
        public List<IComponent<LNode, Set<LEdge>>> getComponents() {
            return components;
        }

        @Override
        public boolean isContainsExternalExtensions() {
            return containsExternalPorts;
        }
    }
    
    /**
     * Basically, this is a wrapper around an {@link LGraph} since the {@link LGraph} represents a
     * connected component, if the user desires to (i.e. {@link LayeredOptions#SEPARATE_CC} is true).
     */
    private class InternalComponent implements IComponent<LNode, Set<LEdge>> {
        
        private LGraph graph;
        private boolean containsRegularNodes;

        private List<ElkRectangle> rectilinearConvexHull;
        private List<IExternalExtension<Set<LEdge>>> externalExtensions = Lists.newArrayList();
        
        InternalComponent(final LGraph graph) {
            this.graph = graph;
            
            containsRegularNodes = false;
            for (LNode n : graph.getLayerlessNodes()) {
                containsRegularNodes |= n.getType() == NodeType.NORMAL;
            }
        }
        
        @Override
        public Set<PortSide> getExternalExtensionSides() {
            return graph.getProperty(InternalProperties.EXT_PORT_CONNECTIONS);
        }
        
        @Override
        public List<ElkRectangle> getHull() {
            return rectilinearConvexHull;
        }
        
        @Override
        public List<IExternalExtension<Set<LEdge>>> getExternalExtensions() {
            return externalExtensions;
        }
        
        private List<LNode> getNodes() {
            return graph.getLayerlessNodes();
        }
        
        private List<LEdge> getExternalEdges() {
            List<LEdge> edges = Lists.newArrayList();
            for (IExternalExtension<Set<LEdge>> ee : externalExtensions) {
                edges.addAll(ee.getRepresentative());
            }
            return edges;
        }
    }
    
    /**
     * Represents the union of multiple external extensions.
     */
    private final class InternalUnionExternalExtension implements IExternalExtension<Set<LEdge>> {

        private Set<LEdge> edges = Sets.newHashSet();
        private PortSide side;
        private ElkRectangle extension;
        private ElkRectangle placeholder;
        
        @Override
        public Set<LEdge> getRepresentative() {
            return edges;
        }

        @Override
        public ElkRectangle getRepresentor() {
            return extension;
        }

        @Override
        public ElkRectangle getParent() {
            throw new UnsupportedOperationException();
        }

        @Override
        public Direction getDirection() {
            return portSideToDirection(side);
        }

        @Override
        public ElkRectangle getPlaceholder() {
            return placeholder;
        }
        
    }
    
    /**
     * Internal implementation of the {@link IExternalExtension} interface. Contains the
     * {@link LEdge} which represents the external edge alongside which of the two ports is the
     * external port and on which side of the parent node the port is located.
     * 
     * The external edge is split into two parts, the last segment incident to the external port
     * will be represented by an external extension, all other segments will be added to the hull of
     * the component itself.
     * 
     * Care has to be taken when generating the segment representing rectangles. To be uses with the
     * {@link org.eclipse.elk.alg.layered.compaction.oned.algs.ScanlineConstraintCalculator}, they
     * are not allowed to overlap.
     */
    private final class InternalExternalExtension implements IExternalExtension<LEdge> {

        private LEdge edge;
        private LPort externalPort;
        private PortSide externalPortSide;
        
        // the following are calculated 
        private ElkRectangle externalExtension;
        private ElkRectangle parent;
        
        InternalExternalExtension(final LEdge edge) {
            // housekeeping
            this.edge = edge;
            if (edge.getSource().getNode().getType() == NodeType.EXTERNAL_PORT) {
                externalPort = edge.getSource();
                externalPortSide =
                        edge.getSource().getNode().getProperty(InternalProperties.EXT_PORT_SIDE);
            } else if (edge.getTarget().getNode().getType() == NodeType.EXTERNAL_PORT) {
                externalPort = edge.getTarget();
                externalPortSide =
                        edge.getTarget().getNode().getProperty(InternalProperties.EXT_PORT_SIDE);
            } else {
                throw new IllegalArgumentException("Edge " + edge + " is not an external edge.");
            }
        }
        
        @Override
        public Direction getDirection() {
            return portSideToDirection(externalPortSide);
        }
        
        @Override
        public LEdge getRepresentative() {
            return edge;
        }

        @Override
        public ElkRectangle getRepresentor() {
            return externalExtension;
        }
        
        @Override
        public ElkRectangle getParent() {
            return parent;
        }
    }
    
    // ------------------------------------------------------------------------------------------------
    //                                     Further classess
    // ------------------------------------------------------------------------------------------------

    
    /**
     * {@link ArrayList} extension, keeping track of minimal and maximal values.
     */
    private static final class Hullpoints extends ArrayList<Point> {
        
        private static final long serialVersionUID = -8667344007490805272L;

        private KVector topLeft = new KVector(Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY);
        private KVector bottomRight = new KVector(Double.NEGATIVE_INFINITY,
                Double.NEGATIVE_INFINITY);

        @Override
        public boolean add(final Point e) {
            topLeft.x = Math.min(topLeft.x, e.x);
            topLeft.y = Math.min(topLeft.y, e.y);
            bottomRight.x = Math.max(bottomRight.x, e.x);
            bottomRight.y = Math.max(bottomRight.y, e.y);
            return super.add(e);
        }
        
        public boolean add(final KVector e) {
            return add(Point.from(e));
        }
        
        public boolean add(final ElkRectangle rect) {
            boolean returnVal = true;
            returnVal &= add(rect.getPosition());
            returnVal &= add(rect.getPosition().add(rect.width, 0));
            returnVal &= add(rect.getPosition().add(0, rect.height));
            returnVal &= add(rect.getPosition().add(rect.width, rect.height));
            return returnVal;
        }
    }
    
    /**
     * Internal collector for inner and outer segments.
     */
    private static final class Segments {
        private List<Pair<KVector, KVector>> innerSegments = Lists.newArrayList();
        private Pair<KVector, KVector> outerSegment;
    }

    /** Internal collector for extreme points of unions of outer segments. */
    private static class OuterSegments {
        private double[] min = new double[PortSide.values().length];
        private double[] max = new double[PortSide.values().length];
        private double[] extent = new double[PortSide.values().length];

        OuterSegments() {
            // initialize with proper values
            Arrays.fill(min, Double.POSITIVE_INFINITY);
            Arrays.fill(max, Double.NEGATIVE_INFINITY);
            
            // Arrays.fill(extent, Double.POSITIVE_INFINITY);
            Arrays.fill(extent, Double.NEGATIVE_INFINITY);
        }
    }
    
}
