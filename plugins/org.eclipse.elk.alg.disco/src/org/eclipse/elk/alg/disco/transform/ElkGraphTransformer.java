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
package org.eclipse.elk.alg.disco.transform;

import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map.Entry;
import java.util.Set;
import java.util.function.BiConsumer;
import java.util.function.Consumer;

import org.eclipse.elk.alg.disco.graph.DCDirection;
import org.eclipse.elk.alg.disco.graph.DCElement;
import org.eclipse.elk.alg.disco.graph.DCExtension;
import org.eclipse.elk.alg.disco.graph.DCGraph;
import org.eclipse.elk.alg.disco.options.DisCoOptions;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.math.KVectorChain;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.util.ElkUtil;
import org.eclipse.elk.graph.ElkEdge;
import org.eclipse.elk.graph.ElkEdgeSection;
import org.eclipse.elk.graph.ElkGraphElement;
import org.eclipse.elk.graph.ElkLabel;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.ElkPort;
import org.eclipse.elk.graph.ElkShape;
import org.eclipse.elk.graph.util.ElkGraphUtil;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;

/**
 * Class for transforming any graph structure into a {@link DCGraph}, in order to use layout algorithms provided by
 * DisCo for compacting disconnected components of the given graph.
 */
public class ElkGraphTransformer implements IGraphTransformer<ElkNode> {

    ///////////////////////////////////////////////////////////////////////////////
    // Nifty constant

    private static final double HALF_PI = Math.PI / 2;

    ///////////////////////////////////////////////////////////////////////////////
    // Variables

    /** the graph represented by its hierarchical parent {@link ElkNode}. */
    private ElkNode parent;
    /**
     * Maps Elements of original graph to DCGraph elements, (but only if their coordinates have to be changed, when
     * components move. This doesn't concern ElkPorts (their coordinates are relative their ElkNode) and ElkLabels of
     * ElkNodes and ElkPorts).
     */
    private HashMap<ElkGraphElement, DCElement> elementMapping = Maps.newHashMap();

    private HashMap<ElkEdge, DCExtension> incomingExtensionsMapping = Maps.newHashMap();
    private HashMap<ElkEdge, DCExtension> outgoingExtensionsMapping = Maps.newHashMap();

    /** the graph after transformation into the simpler {@link DCGraph} used for layout out connected components. */
    private DCGraph transformedGraph;

    private double componentSpacing = 0.0d;

    ///////////////////////////////////////////////////////////////////////////////
    // Constructor

    /**
     * Creates a new graph transformer with a minimum spacing of 0 between each two components.
     */
    public ElkGraphTransformer() {
        super();
    }

    /**
     * Creates a new graph transformer with a minimum given spacing between each two components.
     * 
     * @param componentSpacing
     *            minimum space among components
     */
    public ElkGraphTransformer(final double componentSpacing) {
        this.componentSpacing = componentSpacing;
    }

    ///////////////////////////////////////////////////////////////////////////////
    // Implementation of interface methods

    @Override
    public DCGraph importGraph(final ElkNode graph) {
        parent = graph;
        // Split graph into components
        List<List<ElkNode>> components = ElkGraphComponentsProcessor.split(graph);
        // This list of lists is used to construct the DCGraph from its smallest units, so called DCElements.
        List<List<DCElement>> result = Lists.newArrayList();

        for (List<ElkNode> component : components) {
            // Each of the subResult lists will represent a connected component in the DCGraph.
            List<DCElement> subResult = Lists.newArrayList();
            result.add(subResult);
            Set<ElkEdge> edgeSet = Sets.newHashSet();
            for (ElkNode node : component) {

                // "true" indicates - we need to consider the position of ElkNodes when applying the new layout of the
                // DCGraph back to the KGraph.
                DCElement componentNode = importElkShape(node, true, 0.0f, 0.0f);
                subResult.add(componentNode);

                // Compute offset for labels (used for computing absolute position)
                double nodeX = node.getX();
                double nodeY = node.getY();

                // For use in the Debugview only
                componentNode.setParentCoords(new KVector(nodeX, nodeY));

                // next look at all labels of this ElkNode ...
                List<ElkLabel> labels = node.getLabels();
                for (ElkLabel label : labels) {
                    // "false" indicates - ElkLabels belonging to nodes have coordinates relative to their ElkNode. So
                    // don't
                    // adjust their position when applying the layout of the DCGraph back to the KGraph.
                    DCElement componentLabel = importElkShape(label, false, nodeX, nodeY);
                    subResult.add(componentLabel);
                }

                // ... then the ports of the ElkNode
                List<ElkPort> ports = node.getPorts();
                for (ElkPort port : ports) {
                    // "false" - ElkPorts have coordinates relative to their ElkNode
                    DCElement componentPort = importElkShape(port, false, nodeX, nodeY);
                    subResult.add(componentPort);

                    // ElkPorts can have labels, too!
                    // Compute offset for labels of ElkPorts
                    double portX = port.getX() + nodeX;
                    double portY = port.getY() + nodeY;

                    labels = port.getLabels();
                    for (ElkLabel label : labels) {
                        // false - ElkLabels have coordinates relative to their ElkPort
                        DCElement componentLabel = importElkShape(label, false, portX, portY);
                        subResult.add(componentLabel);
                    }

                }

                edgeSet.addAll(Sets.newHashSet(ElkGraphUtil.allIncidentEdges(node)));
            }
            importElkEdges(edgeSet, subResult);
        }
        // Finally create the DCGraph
        transformedGraph = new DCGraph(result, componentSpacing / 2);

        // Don't forget to copy properties of parent ElkNode to the DCGraph
        transformedGraph.copyProperties(graph);

        return transformedGraph;
    }

    @Override
    public void applyLayout() {
        KVector graphDimensions = transformedGraph.getDimensions();
        double newWidth = graphDimensions.x;
        double newHeight = graphDimensions.y;

        double oldWidth = parent.getWidth();
        double oldHeight = parent.getHeight();

        // Adjust size of layout
        parent.setDimensions(graphDimensions.x, graphDimensions.y);

        double xFactor = newWidth / oldWidth;
        double yFactor = newHeight / oldHeight;

        for (ElkLabel label : parent.getLabels()) {
            label.setX(label.getX() * xFactor);
            label.setY(label.getY() * yFactor);
        }

        for (ElkPort port : parent.getPorts()) {
            double px = port.getX();
            double py = port.getY();

            if (px > 0) {
                port.setX(px * xFactor);
            }
            if (py > 0) {
                port.setY(py * yFactor);

            }

        }

        // Apply offsets, whenever necessary.
        elementMapping.forEach(new OffsetApplier());

        List<ElkPort> adjustedPorts = Lists.newArrayList();

        ElkPort portToAdjust;
        for (Entry<ElkEdge, DCExtension> inEntry : incomingExtensionsMapping.entrySet()) {
            ElkEdge edge = inEntry.getKey();
            DCDirection dir = inEntry.getValue().getDirection();
            ElkEdgeSection edgeSection = ElkGraphUtil.firstEdgeSection(edge, false, false);
            KVectorChain newPoints =
                    adjustFirstSegment(ElkGraphUtil.getSourceNode(edge), ElkUtil.createVectorChain(edgeSection), dir);
            ElkUtil.applyVectorChain(newPoints, edgeSection);

            portToAdjust = ElkGraphUtil.getSourcePort(edge);
            if (portToAdjust != null && !adjustedPorts.contains(portToAdjust)) {
                adjustedPorts.add(portToAdjust);
                adjustRelatedPort(portToAdjust, newPoints.getFirst(), dir);
            }
        }
        for (Entry<ElkEdge, DCExtension> outEntry : outgoingExtensionsMapping.entrySet()) {
            ElkEdge edge = outEntry.getKey();
            DCDirection dir = outEntry.getValue().getDirection();
            ElkEdgeSection edgeSection = ElkGraphUtil.firstEdgeSection(edge, false, false);
            KVectorChain newPoints = adjustFirstSegment(ElkGraphUtil.getTargetNode(edge),
                    KVectorChain.reverse(ElkUtil.createVectorChain(edgeSection)), dir);
            newPoints = KVectorChain.reverse(newPoints);
            ElkUtil.applyVectorChain(newPoints, edgeSection);

            portToAdjust = ElkGraphUtil.getTargetPort(edge);
            if (portToAdjust != null && !adjustedPorts.contains(portToAdjust)) {
                adjustedPorts.add(portToAdjust);
                adjustRelatedPort(portToAdjust, newPoints.getLast(), dir);
            }
        }

    }

    ///////////////////////////////////////////////////////////////////////////////
    // Private helper methods

    /**
     * Adjusts the position of a port after compaction. Note that DisCo currently only supports the fixed side port
     * constraint.
     * 
     * @param port
     *            Port to be adjusted
     * @param edgePoint
     *            The end point of an edge (extension) connected to the port after compaction took place
     * @param dir
     *            Cardinal direction the extension is facing, which is at the same time the side of the parent node the
     *            port is placed on.
     */
    private void adjustRelatedPort(final ElkPort port, final KVector edgePoint, final DCDirection dir) {
        if (dir.isHorizontal()) {
            port.setY(edgePoint.y - (port.getHeight() / 2));
        } else { // dir is vertical
            port.setX(edgePoint.x - (port.getWidth() / 2));
        }
    }

    /**
     * Transforms a {@link ElkNode}, {@link ElkLabel} or {@link ElkPort} into a {@link DCElement} without destroying it.
     * An offset may be applied, as some of these elements may have coordinates relative to another element, but the
     * {@link DCGraph} only works with absolute coordinates.
     * 
     * @param element
     *            Elemment to be transformed
     * @param considerWhenApplyingOffset
     *            true - a key-value pair consisting of the original element and its {@link DCElement} representation
     *            will be saved and used later by {@link ElkGraphTransformer#applyLayout()} to set the original element
     *            to its new position after the layouting process has finished. false - the original element has
     *            coordinates relative to another element already moved and doesn't need new coordinates after
     *            layouting.
     * @param offsetX
     *            X-Coordinate of the offset needed to make the given element's coordinates absolute (only needed, when
     *            the element has coordinates not relative to its parent {@link ElkNode}).
     * @param offsetY
     *            Y-Coordinate of the offset needed to make the given element's coordinates absolute (only needed, when
     *            the element has coordinates not relative to its parent {@link ElkNode}).
     * @return result of transformation
     * @throws IllegalArgumentException
     *             A {@link ElkGraphElement} other than {@link ElkNode}, {@link ElkLabel} or {@link ElkPort} has been
     *             given as an argument.
     */
    private <E extends ElkShape> DCElement importElkShape(final E element, final boolean considerWhenApplyingOffset,
            final double offsetX, final double offsetY) throws IllegalArgumentException {
        if (!(element instanceof ElkNode || element instanceof ElkLabel || element instanceof ElkPort)) {
            throw new IllegalArgumentException("Method only works for ElkNode-, ElkLabel and ElkPort-objects.");
        }

        double halfComponentSpacing = componentSpacing / 2;
        double x0 = element.getX() + offsetX - halfComponentSpacing;
        double y0 = element.getY() + offsetY - halfComponentSpacing;
        double x1 = x0 + element.getWidth() + componentSpacing;
        double y1 = y0 + element.getHeight() + componentSpacing;

        KVectorChain coords = new KVectorChain();
        coords.add(newPoint(x0, y0));
        coords.add(newPoint(x0, y1));
        coords.add(newPoint(x1, y1));
        coords.add(newPoint(x1, y0));

        DCElement shape = new DCElement(coords);
        // copy all properties of original layout
        shape.copyProperties(element);

        if (considerWhenApplyingOffset) {
            elementMapping.put(element, shape);
        }

        return shape;
    }

    /**
     * Transforms a {@link ElkEdge} into a {@link DCElement} without destroying it. Edges can have their own
     * {@link ElkLabel ElkLabels}, so they will be transformed, too.
     * 
     * @param edge
     *            Edge to be transformed into a {@link DCElement}
     * @param newComponent
     *            Collection representing the component the edge and its associated labels (if any) belong to. Newly
     *            generated {@link DCElement DCElements} will be added to it
     * @return {@link DCElement} resulting from the transformation
     */
    private DCElement importElkEdge(final ElkEdge edge, final Collection<DCElement> newComponent) {
        ElkEdgeSection edgeSection = ElkGraphUtil.firstEdgeSection(edge, false, false);
        List<KVector> points = ElkUtil.createVectorChain(edgeSection);

        double thickness = edge.getProperty(DisCoOptions.EDGE_THICKNESS);
        KVectorChain contour = getContour(points, thickness + componentSpacing);

        DCElement shape = new DCElement(contour);
        shape.copyProperties(edge);
        elementMapping.put(edge, shape);
        newComponent.add(shape);

        // ElkEdges can have labels, too!
        List<ElkLabel> labels = edge.getLabels();
        for (ElkLabel label : labels) {
            // "true" - ElkLabels belonging to an ElkEdge have absolute coordinates and have to be considered when
            // applying
            // changes to the DCGraph back to the original graph.
            DCElement componentLabel = importElkShape(label, true, 0.0f, 0.0f);
            newComponent.add(componentLabel);
        }

        return shape;
    }

    /**
     * Transforms a {@link ElkEdge} into a {@link DCElement} without destroying it. Edges can have their own
     * {@link ElkLabel ElkLabels}, so they will be transformed, too.
     * 
     * @param edge
     *            Edge to be transformed into a {@link DCElement}
     * @param newComponent
     *            Collection representing the component the edge and its associated labels (if any) belong to. Newly
     *            generated {@link DCElement DCElements} will be added to it
     * @return {@link DCElement} resulting from the transformation
     */
    /**
     * Transforms a short hierarchical {@link ElkEdge} into a {@link DCExtension} without destroying it. Edges can have
     * their own {@link ElkLabel ElkLabels}, so they will be transformed, too.
     * 
     * @param edge
     *            Edge to be transformed into a {@link DCExtension}
     * @param newComponent
     *            component the extension will belong to
     * @param outgoingExtension
     *            true, if the edge should be handled as an outgoing extension; false, otherwise
     */
    private void importExtension(final ElkEdge edge, final Collection<DCElement> newComponent,
            final boolean outgoingExtension) {
        ElkEdgeSection edgeSection = ElkGraphUtil.firstEdgeSection(edge, false, false);
        KVectorChain points = ElkUtil.createVectorChain(edgeSection);

        if (outgoingExtension) {
            points = KVectorChain.reverse(points);
        }

        // Is going to hold the extension
        DCElement shape;

        double thickness = edge.getProperty(DisCoOptions.EDGE_THICKNESS);

        KVector outerPoint = points.getFirst();
        KVector innerPoint = points.get(1);

        if (points.size() > 2) {
            List<KVector> fixedEdgePoints = Lists.newArrayList();
            fixedEdgePoints.addAll(points.subList(1, points.size()));

            KVectorChain contour = getContour(fixedEdgePoints, thickness + componentSpacing);

            shape = new DCElement(contour);
            shape.copyProperties(edge);
            newComponent.add(shape);
        } else {
            if (outgoingExtension) {
                shape = elementMapping.get(ElkGraphUtil.getSourceNode(edge));
            } else {
                shape = elementMapping.get(ElkGraphUtil.getTargetNode(edge));
            }
        }

        // Construct the extension and add to mapping
        ElkNode extParent = ElkGraphUtil.getSourceNode(edge);
        if (outgoingExtension) {
            extParent = ElkGraphUtil.getTargetNode(edge);
        }
        DCDirection dir = nearestSide(outerPoint, extParent);
        double extensionWidth = thickness + componentSpacing;
        KVector middlePos;

        if (dir.isHorizontal()) { // West or east extension
            extensionWidth += Math.abs(outerPoint.y - innerPoint.y);
            middlePos = new KVector(innerPoint.x, (innerPoint.y + outerPoint.y) / 2);
        } else {
            extensionWidth += Math.abs(outerPoint.x - innerPoint.x);
            middlePos = new KVector((innerPoint.x + outerPoint.x) / 2, innerPoint.y);
        }

        if (outgoingExtension) {
            outgoingExtensionsMapping.put(edge, new DCExtension(shape, dir, middlePos, extensionWidth));
        } else {
            incomingExtensionsMapping.put(edge, new DCExtension(shape, dir, middlePos, extensionWidth));
        }
        elementMapping.put(edge, shape);

        // ElkEdges can have labels, too!
        List<ElkLabel> labels = edge.getLabels();
        for (ElkLabel label : labels) {
            // "true" - ElkLabels belonging to an ElkEdge have absolute coordinates and have to be considered when
            // applying
            // changes to the DCGraph back to the original graph.
            DCElement componentLabel = importElkShape(label, true, 0.0f, 0.0f);
            newComponent.add(componentLabel);
        }

    }

    /**
     * Computes for a given {@link KVector} relative to the upper left corner of a given node, which of its four sides
     * is closest to the point the vector describes.
     * 
     * @param point
     *            vector to compare
     * @param node
     *            node as reference to the vector
     * @return side of the node which is closest to the point
     */
    private DCDirection nearestSide(final KVector point, final ElkNode node) {
        double distance;
        double shortestDistance = Double.MAX_VALUE;
        DCDirection result = DCDirection.NORTH;

        // NORTHVALUE
        shortestDistance = Math.abs(point.y);
        // SOUTHVALUE
        distance = Math.abs(node.getHeight() - point.y);
        if (distance < shortestDistance) {
            shortestDistance = distance;
            result = DCDirection.SOUTH;
        }
        // WESTVALUE
        distance = Math.abs(point.x);
        if (distance < shortestDistance) {
            shortestDistance = distance;
            result = DCDirection.WEST;
        }
        // EASTVALUE
        distance = Math.abs(node.getWidth() - point.x);
        if (distance < shortestDistance) {
            shortestDistance = distance;
            result = DCDirection.EAST;
        }
        return result;

    }

    /**
     * Transforms multiple {@link ElkEdge ElkEdges} (and their labels) into {@link DCElement DCElements} and adds them
     * to a given collection, which is representing a connected component in the final {@link DCGraph}.
     * 
     * @param edges
     *            Edges to be transformed, all from the same connected component
     * @param newComponent
     *            Collection representing one connected component of the graph
     */
    private void importElkEdges(final Collection<ElkEdge> edges, final Collection<DCElement> newComponent) {
        for (ElkEdge edge : edges) {
            // For Edges, you always have to check, if the corresponding DCElement already exists.
            DCElement componentEdge = elementMapping.get(edge);
            if (componentEdge == null) {
                if (ElkGraphUtil.getSourceNode(edge).getParent().equals(ElkGraphUtil.getTargetNode(edge).getParent())) {
                    importElkEdge(edge, newComponent);
                } else {
                    if (ElkGraphUtil.getSourceNode(edge).equals(ElkGraphUtil.getTargetNode(edge).getParent())) {
                        if (incomingExtensionsMapping.get(edge) == null
                                && elementMapping.get(ElkGraphUtil.getTargetNode(edge)) != null) { // otherwise the
                                                                                                   // extension has
                            // already been taken care of
                            importExtension(edge, newComponent, false);
                        }
                    } else {
                        // outGoingExtensions
                        if (outgoingExtensionsMapping.get(edge) == null
                                && elementMapping.get(ElkGraphUtil.getSourceNode(edge)) != null) { // otherwise the
                                                                                                   // extension has
                            // already been taken care of
                            importExtension(edge, newComponent, true);
                        }
                    }
                }
            }
        }
    }

    /**
     * Given all source, bend and target points of an edge in correct order, this method computes all vertices needed to
     * construct a polygonal {@link DCElement} for it.
     * 
     * @param edgePoints
     *            All relevant points of an Edge as described above
     * @param thickness
     *            Desired thickness of the edge
     * @return List of vertices of the edge as a polygon considering the given thickness
     */
    private KVectorChain getContour(final List<KVector> edgePoints, final double thickness) {
        List<KVector> ccwPoints = Lists.newArrayList();
        List<KVector> cwPoints = Lists.newArrayList();
        double radius = thickness / 2;

        int numberOfPoints = edgePoints.size();

        // special case: first edge segment. Assumption an edge has at least two points
        KVector predecessor;
        KVector current = edgePoints.get(0);
        KVector successor = edgePoints.get(1);
        List<KVector> orthPoints = getOrthogonalPoints(current.x, current.y, successor.x, successor.y, radius);
        ccwPoints.add(orthPoints.get(0));
        cwPoints.add(orthPoints.get(1));

        // normal case: bendpoints have a preceding and a succeeding neighbor
        for (int i = 2; i < numberOfPoints; i++) {
            predecessor = current;
            current = successor;
            successor = edgePoints.get(i);
            orthPoints = getOrthogonalPoints(current.x, current.y, predecessor.x, predecessor.y, radius);
            ccwPoints.add(orthPoints.get(1));
            cwPoints.add(orthPoints.get(0));

            orthPoints = getOrthogonalPoints(current.x, current.y, successor.x, successor.y, radius);
            ccwPoints.add(orthPoints.get(0));
            cwPoints.add(orthPoints.get(1));
        }
        // For the last point, consider the line connecting back to the previous point (the resulting orthogonal points
        // will be in reversed order in the resulting list, though.).
        orthPoints = getOrthogonalPoints(successor.x, successor.y, current.x, current.y, radius);
        ccwPoints.add(orthPoints.get(1));
        cwPoints.add(orthPoints.get(0));

        // Compute the intersections of the line segments of the orthogonal points
        KVectorChain ccwMerged = new KVectorChain();
        List<KVector> cwMerged = Lists.newArrayList();

        ccwMerged.add(ccwPoints.get(0));
        for (int i = 1; i < ccwPoints.size() - 2; i += 2) {
            KVector currentPoint = ccwPoints.get(i);
            KVector intersectionPoint =
                    computeIntersection(ccwPoints.get(i - 1), currentPoint, ccwPoints.get(i + 1), ccwPoints.get(i + 2));

            // Keep currentPoint if lines are coincident (or parallel).
            if (!Double.isFinite(intersectionPoint.x) || !Double.isFinite(intersectionPoint.y)) {
                ccwMerged.add(currentPoint);
            } else {
                ccwMerged.add(intersectionPoint);
            }
        }
        ccwMerged.add(ccwPoints.get(ccwPoints.size() - 1));

        cwMerged.add(cwPoints.get(0));
        for (int i = 1; i < cwPoints.size() - 2; i += 2) {
            KVector currentPoint = cwPoints.get(i);
            KVector intersectionPoint =
                    computeIntersection(cwPoints.get(i - 1), currentPoint, cwPoints.get(i + 1), cwPoints.get(i + 2));

            // Keep currentPoint if lines are coincident (or parallel).
            if (!Double.isFinite(intersectionPoint.x) || !Double.isFinite(intersectionPoint.y)) {
                cwMerged.add(currentPoint);
            } else {
                cwMerged.add(intersectionPoint);
            }
        }
        cwMerged.add(cwPoints.get(cwPoints.size() - 1));

        // merge lists (one of them has to be added in reverse order)
        for (int i = cwMerged.size() - 1; i >= 0; i--) {
            ccwMerged.add(cwMerged.get(i));
        }

        return ccwMerged;
    }

    /**
     * Creates two points from one point (of an edge) orthogonal to the edges orientation with a distance of
     * {@code radius} on each side of the edge.
     * 
     * @param curX
     *            X-coordinate of current source, bend or target point on the edge
     * @param curY
     *            Y-coordinate of current source, bend or target point on the edge
     * @param nxtX
     *            X-coordinate of next source, bend or target point on the edge
     * @param nxtY
     *            Y-coordinate of next source, bend or target point on the edge
     * @param radius
     *            Distance of the generated point from the current point
     * @return Generated points as a list
     */
    private List<KVector> getOrthogonalPoints(final double curX, final double curY, final double nxtX,
            final double nxtY, final double radius) {
        double difX = nxtX - curX;
        double difY = nxtY - curY;

        double angleRadians = Math.atan2(difX, difY);
        double orthAngleCCW = angleRadians + HALF_PI;
        double orthAngleCW = angleRadians - HALF_PI;

        double xCCW = radius * Math.sin(orthAngleCCW) + curX;
        double yCCW = radius * Math.cos(orthAngleCCW) + curY;
        double xCW = radius * Math.sin(orthAngleCW) + curX;
        double yCW = radius * Math.cos(orthAngleCW) + curY;

        return Lists.newArrayList(newPoint(xCCW, yCCW), newPoint(xCW, yCW));
    }

    // might result in +inf or -inf
    private KVector computeIntersection(final KVector p1, final KVector p2,
            final KVector p3, final KVector p4) {
        double x1 = p1.x;
        double y1 = p1.y;
        double x2 = p2.x;
        double y2 = p2.y;
        double x3 = p3.x;
        double y3 = p3.y;
        double x4 = p4.x;
        double y4 = p4.y;

        double factor1 = x1 * y2 - y1 * x2;
        double factor2 = x3 * y4 - y3 * x4;
        double denominator = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);

        double x = (factor1 * (x3 - x4) - factor2 * (x1 - x2)) / denominator;
        double y = (factor1 * (y3 - y4) - factor2 * (y1 - y2)) / denominator;

        return newPoint(x, y);
    }

    /**
     * Creates a new instance of {@link KVector}.
     * 
     * @param x
     *            X-coordinate
     * @param y
     *            Y-coordinate
     * @return Freshly instantiated point
     */
    private KVector newPoint(final double x, final double y) {
        return new KVector(x, y);
    }

    /**
     * Applies all changes to an edge segment which has been treated as an extension during compaction after the
     * compaction process has finished.
     * 
     * @param source
     *            node the edge belongs to
     * @param chain
     *            edge in its vector chain form
     * @param dir
     *            direction the extension points to
     * @return adjusted vector chain
     */
    private KVectorChain adjustFirstSegment(final ElkNode source, final KVectorChain chain, final DCDirection dir) {
        KVector firstPoint = chain.remove();
        switch (dir) {
        case NORTH:
            firstPoint.y = 0;
            break;
        case SOUTH:
            firstPoint.y = source.getHeight();
            break;
        case WEST:
            firstPoint.x = 0;
            break;
        default: // EAST
            firstPoint.x = source.getWidth();
        }

        chain.add(0, firstPoint);

        return chain;

    }

    ///////////////////////////////////////////////////////////////////////////////
    // Inner nested helper classes

    /**
     * Class implementing a function for mapping over a {@link HashMap}, in this case applying the layout of a
     * {@link DCElement} to the original KGraph-element it represents.
     */
    private class OffsetApplier implements BiConsumer<ElkGraphElement, DCElement> {
        private KVector offset;

        /*
         * (non-Javadoc)
         *
         * @see java.util.function.BiConsumer#accept(java.lang.Object, java.lang.Object)
         */
        @Override
        public void accept(final ElkGraphElement elem, final DCElement poly) {
            offset = poly.getOffset();
            if (elem instanceof ElkEdge) {
                ElkEdgeSection edgeSection = ElkGraphUtil.firstEdgeSection((ElkEdge) elem, false, false);
                KVectorChain points = ElkUtil.createVectorChain(edgeSection);
                OffSetToChainApplier applier = new OffSetToChainApplier();
                points.forEach(applier);
                ElkUtil.applyVectorChain(points, edgeSection);
                if (elem.getProperty(CoreOptions.JUNCTION_POINTS) != null) {
                    elem.getProperty(CoreOptions.JUNCTION_POINTS).forEach(applier);
                }
            } else {
                ElkShape shape = (ElkShape) elem;
                shape.setX(shape.getX() + offset.x);
                shape.setY(shape.getY() + offset.y);
            }
        }

        /**
         * Class implementing a function for mapping over an {@link Iterable}, in this case adding a specific
         * {@link KVector} to a lot of (different) {@link KVector KVectors}.
         */
        private class OffSetToChainApplier implements Consumer<KVector> {

            /*
             * (non-Javadoc)
             *
             * @see java.util.function.Consumer#accept(java.lang.Object)
             */
            @Override
            public void accept(final KVector point) {
                point.add(offset.x, offset.y);
            }

        }

    }
}
