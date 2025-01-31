/*******************************************************************************
 * Copyright (c) 2008, 2020 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.core.util;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.ListIterator;
import java.util.Map;
import java.util.Set;
import java.util.regex.Pattern;

import org.eclipse.elk.core.math.ElkMargin;
import org.eclipse.elk.core.math.ElkRectangle;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.math.KVectorChain;
import org.eclipse.elk.core.options.ContentAlignment;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.options.Direction;
import org.eclipse.elk.core.options.EdgeLabelPlacement;
import org.eclipse.elk.core.options.NodeLabelPlacement;
import org.eclipse.elk.core.options.PortConstraints;
import org.eclipse.elk.core.options.PortLabelPlacement;
import org.eclipse.elk.core.options.PortSide;
import org.eclipse.elk.core.options.SizeConstraint;
import org.eclipse.elk.core.options.SizeOptions;
import org.eclipse.elk.core.util.adapters.GraphAdapters.LabelAdapter;
import org.eclipse.elk.core.util.adapters.GraphAdapters.PortAdapter;
import org.eclipse.elk.core.validation.GraphIssue;
import org.eclipse.elk.core.validation.GraphValidationException;
import org.eclipse.elk.core.validation.IValidatingGraphElementVisitor;
import org.eclipse.elk.graph.ElkBendPoint;
import org.eclipse.elk.graph.ElkConnectableShape;
import org.eclipse.elk.graph.ElkEdge;
import org.eclipse.elk.graph.ElkEdgeSection;
import org.eclipse.elk.graph.ElkGraphElement;
import org.eclipse.elk.graph.ElkGraphFactory;
import org.eclipse.elk.graph.ElkLabel;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.ElkPort;
import org.eclipse.elk.graph.ElkShape;
import org.eclipse.elk.graph.util.ElkGraphUtil;
import org.eclipse.emf.ecore.EObject;

import com.google.common.base.Strings;
import com.google.common.collect.Iterables;
import com.google.common.collect.Iterators;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

/**
 * Utility methods for layout-related things.
 */
public final class ElkUtil {

    /**
     * Default minimal width for nodes.
     */
    public static final float DEFAULT_MIN_WIDTH = 20.0f;

    /**
     * Default minimal height for nodes.
     */
    public static final float DEFAULT_MIN_HEIGHT = 20.0f;


    /**
     * Hidden constructor to avoid instantiation.
     */
    private ElkUtil() {
    }


    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // LAYOUT-RELATED METHODS

    /**
     * Determines the port side for the given port from its relative position at
     * its corresponding node.
     *
     * @param port port to analyze
     * @param direction the overall layout direction
     * @return the port side relative to its containing node
     * @throws IllegalStateException if the port does not have a parent node.
     */
    public static PortSide calcPortSide(final ElkPort port, final Direction direction) {
        if (port.getParent() == null) {
            throw new IllegalStateException("port must have a parent node to calculate the port side");
        }

        // if the node has zero size, we cannot decide anything
        ElkNode node = port.getParent();
        double nodeWidth = node.getWidth(), nodeHeight = node.getHeight();
        if (nodeWidth <= 0 && nodeHeight <= 0) {
            return PortSide.UNDEFINED;
        }

        // check direction-dependent criterion
        double xpos = port.getX(), ypos = port.getY();
        switch (direction) {
        case LEFT:
        case RIGHT:
            if (xpos < 0) {
                return PortSide.WEST;
            } else if (xpos + port.getWidth() > nodeWidth) {
                return PortSide.EAST;
            }
            break;
        case UP:
        case DOWN:
            if (ypos < 0) {
                return PortSide.NORTH;
            } else if (ypos + port.getHeight() > nodeHeight) {
                return PortSide.SOUTH;
            }
        }

        // check general criterion
        double widthPercent = (xpos + port.getWidth() / 2) / nodeWidth;
        double heightPercent = (ypos + port.getHeight() / 2) / nodeHeight;
        if (widthPercent + heightPercent <= 1 && widthPercent - heightPercent <= 0) {
            // port is on the left
            return PortSide.WEST;
        } else if (widthPercent + heightPercent >= 1 && widthPercent - heightPercent >= 0) {
            // port is on the right
            return PortSide.EAST;
        } else if (heightPercent < 1.0f / 2) {
            // port is on the top
            return PortSide.NORTH;
        } else {
            // port is on the bottom
            return PortSide.SOUTH;
        }
    }

    /**
     * Calculate the offset for a port, that is the amount by which it is moved outside of the node.
     * An offset value of 0 means the port has no intersection with the node and touches the outside
     * border of the node.
     *
     * @param port a port
     * @param side the side on the node for the given port
     * @return the offset on the side
     * @throws IllegalStateException if the port does not have a parent node.
     */
    public static double calcPortOffset(final ElkPort port, final PortSide side) {
        if (port.getParent() == null) {
            throw new IllegalStateException("port must have a parent node to calculate the port side");
        }

        ElkNode node = port.getParent();
        switch (side) {
        case NORTH:
            return -(port.getY() + port.getHeight());
        case EAST:
            return port.getX() - node.getWidth();
        case SOUTH:
            return port.getY() - node.getHeight();
        case WEST:
            return -(port.getX() + port.getWidth());
        }
        return 0;
    }

    /**
     * Sets the size of a given node, depending on the minimal size, the number of ports
     * on each side and the label.
     *
     * @param node the node that shall be resized
     * @return a vector holding the width and height resizing ratio, or {@code null} if the size
     *     constraint is set to {@code FIXED}
     */
    public static KVector resizeNode(final ElkNode node) {
        Set<SizeConstraint> sizeConstraint = node.getProperty(CoreOptions.NODE_SIZE_CONSTRAINTS);
        if (sizeConstraint.isEmpty()) {
            return null;
        }

        double newWidth = 0, newHeight = 0;

        if (sizeConstraint.contains(SizeConstraint.PORTS)) {
            PortConstraints portConstraints = node.getProperty(CoreOptions.PORT_CONSTRAINTS);
            double minNorth = 2, minEast = 2, minSouth = 2, minWest = 2;
            Direction direction = node.getParent() == null
                    ? node.getProperty(CoreOptions.DIRECTION)
                    : node.getParent().getProperty(CoreOptions.DIRECTION);
            for (ElkPort port : node.getPorts()) {
                PortSide portSide = port.getProperty(CoreOptions.PORT_SIDE);
                if (portSide == PortSide.UNDEFINED) {
                    portSide = calcPortSide(port, direction);
                    port.setProperty(CoreOptions.PORT_SIDE, portSide);
                }
                if (portConstraints == PortConstraints.FIXED_POS) {
                    switch (portSide) {
                    case NORTH:
                        minNorth = Math.max(minNorth, port.getX() + port.getWidth());
                        break;
                    case EAST:
                        minEast = Math.max(minEast, port.getY() + port.getHeight());
                        break;
                    case SOUTH:
                        minSouth = Math.max(minSouth, port.getX() + port.getWidth());
                        break;
                    case WEST:
                        minWest = Math.max(minWest, port.getY() + port.getHeight());
                        break;
                    }
                } else {
                    switch (portSide) {
                    case NORTH:
                        minNorth += port.getWidth() + 2;
                        break;
                    case EAST:
                        minEast += port.getHeight() + 2;
                        break;
                    case SOUTH:
                        minSouth += port.getWidth() + 2;
                        break;
                    case WEST:
                        minWest += port.getHeight() + 2;
                        break;
                    }
                }
            }

            newWidth = Math.max(minNorth, minSouth);
            newHeight = Math.max(minEast, minWest);
        }

        return resizeNode(node, newWidth, newHeight, true, true);
    }

    /**
     * Resize a node to the given width and height, adjusting port and label positions if needed.
     *
     * @param node a node
     * @param newWidth the new width to set
     * @param newHeight the new height to set
     * @param movePorts whether port positions should be adjusted
     * @param moveLabels whether label positions should be adjusted
     * @return a vector holding the width and height resizing ratio
     */
    public static KVector resizeNode(final ElkNode node, final double newWidth, final double newHeight,
            final boolean movePorts, final boolean moveLabels) {
        
        KVector oldSize = new KVector(node.getWidth(), node.getHeight());
        
        KVector newSize = effectiveMinSizeConstraintFor(node);
        newSize.x = Math.max(newSize.x, newWidth);
        newSize.y = Math.max(newSize.y, newHeight);

        double widthRatio = newSize.x / oldSize.x;
        double heightRatio = newSize.y / oldSize.y;
        double widthDiff = newSize.x - oldSize.x;
        double heightDiff = newSize.y - oldSize.y;

        // update port positions
        if (movePorts) {
            Direction direction = node.getParent() == null
                    ? node.getProperty(CoreOptions.DIRECTION)
                    : node.getParent().getProperty(CoreOptions.DIRECTION);
            boolean fixedPorts = node.getProperty(CoreOptions.PORT_CONSTRAINTS) == PortConstraints.FIXED_POS;

            for (ElkPort port : node.getPorts()) {
                PortSide portSide = port.getProperty(CoreOptions.PORT_SIDE);

                if (portSide == PortSide.UNDEFINED) {
                    portSide = calcPortSide(port, direction);
                    port.setProperty(CoreOptions.PORT_SIDE, portSide);
                }

                switch (portSide) {
                case NORTH:
                    if (!fixedPorts) {
                        port.setX(port.getX() * widthRatio);
                    }
                    break;
                case EAST:
                    port.setX(port.getX() + widthDiff);
                    if (!fixedPorts) {
                        port.setY(port.getY() * heightRatio);
                    }
                    break;
                case SOUTH:
                    if (!fixedPorts) {
                        port.setX(port.getX() * widthRatio);
                    }
                    port.setY(port.getY() + heightDiff);
                    break;
                case WEST:
                    if (!fixedPorts) {
                        port.setY(port.getY() * heightRatio);
                    }
                    break;
                }
            }
        }

        // resize the node AFTER ports have been placed, since calcPortSide needs the old size
        node.setDimensions(newSize.x, newSize.y);

        // update label positions
        if (moveLabels) {
            for (ElkLabel label : node.getLabels()) {
                double midx = label.getX() + label.getWidth() / 2;
                double midy = label.getY() + label.getHeight() / 2;
                double widthPercent = midx / oldSize.x;
                double heightPercent = midy / oldSize.y;

                if (widthPercent + heightPercent >= 1) {
                    if (widthPercent - heightPercent > 0 && midy >= 0) {
                        // label is on the right
                        label.setX(label.getX() + widthDiff);
                        label.setY(label.getY() + heightDiff * heightPercent);
                    } else if (widthPercent - heightPercent < 0 && midx >= 0) {
                        // label is on the bottom
                        label.setX(label.getX() + widthDiff * widthPercent);
                        label.setY(label.getY() + heightDiff);
                    }
                }
            }
        }

        // set fixed size option for the node: now the size is assumed to stay as determined here
        node.setProperty(CoreOptions.NODE_SIZE_CONSTRAINTS, SizeConstraint.fixed());

        return new KVector(widthRatio, heightRatio);
    }
    
    /**
     * Returns the minimum size of the node according to the {@link CoreOptions#NODE_SIZE_MINIMUM} constraint. If that
     * constraint is not set, the size returned by this method will be {@code (0, 0)}.
     * 
     * @param node the node whose minimum size to compute.
     * @return the minimum size.
     */
    public static KVector effectiveMinSizeConstraintFor(final ElkNode node) {
        Set<SizeConstraint> sizeConstraint = node.getProperty(CoreOptions.NODE_SIZE_CONSTRAINTS);
        
        if (sizeConstraint.contains(SizeConstraint.MINIMUM_SIZE)) {
            Set<SizeOptions> sizeOptions = node.getProperty(CoreOptions.NODE_SIZE_OPTIONS);
            KVector minSize = new KVector(node.getProperty(CoreOptions.NODE_SIZE_MINIMUM));

            // If minimum width or height are not set, maybe default to default values
            if (sizeOptions.contains(SizeOptions.DEFAULT_MINIMUM_SIZE)) {
                if (minSize.x <= 0) {
                    minSize.x = DEFAULT_MIN_WIDTH;
                }
                if (minSize.y <= 0) {
                    minSize.y = DEFAULT_MIN_HEIGHT;
                }
            }

            return minSize;
            
        } else {
            return new KVector();
        }
    }

    /**
     * Applies the scaling factor configured in terms of {@link CoreOptions#SCALE_FACTOR} to {@code node}'s
     * size data, and updates the layout data of {@code node}'s ports and labels accordingly.<br>
     * <b>Note:</b> The scaled layout data won't be reverted during the layout process, see
     * {@link CoreOptions#SCALE_FACTOR}.
     *
     * @param node
     *            the node to be scaled
     */
    public static void applyConfiguredNodeScaling(final ElkNode node) {
        final double scalingFactor = node.getProperty(CoreOptions.SCALE_FACTOR);
        if (scalingFactor == 1) {
            return;
        }

        node.setDimensions(scalingFactor * node.getWidth(), scalingFactor * node.getHeight());

        final Iterable<ElkLabel> portLabels = Iterables.concat(
                Iterables.transform(node.getPorts(), p -> p.getLabels()));
        for (ElkShape shape : Iterables.concat(node.getLabels(), node.getPorts(), portLabels)) {
            shape.setLocation(scalingFactor * shape.getX(), scalingFactor * shape.getY());
            shape.setDimensions(scalingFactor * shape.getWidth(), scalingFactor * shape.getHeight());

            final KVector anchor = shape.getProperty(CoreOptions.PORT_ANCHOR);
            if (anchor != null) {
                anchor.x *= scalingFactor;
                anchor.y *= scalingFactor;
            }
        }
    }
    
    /**
     * Computes the area occupied by this node's layout and stores the values in {@link CoreOptions#CHILD_AREA_WIDTH} 
     * and {@link CoreOptions#CHILD_AREA_HEIGHT}.
     * @param node the node whose child area should be computed
     */
    public static void computeChildAreaDimensions(final ElkNode node) {

        double minX = Double.POSITIVE_INFINITY;
        double minY = Double.POSITIVE_INFINITY;
        double maxX = 0.0;
        double maxY = 0.0;
        
        // iterate over all nodes and labels and get their coordinate bounds
        Iterable<ElkLabel> edgeLabels = new ArrayList<>();
        for (ElkEdge edge : node.getContainedEdges()) {
            edgeLabels = Iterables.concat(edgeLabels, edge.getLabels());
        }
        for (ElkShape shape : Iterables.concat(node.getLabels(), node.getChildren(), edgeLabels)) {

            ElkMargin margins = shape.getProperty(CoreOptions.MARGINS);
            if (minX > shape.getX() - margins.left) {
                minX = shape.getX() - margins.left;
            }
            if (minY > shape.getY() - margins.top) {
                minY = shape.getY() - margins.top;
            }
            if (maxX < shape.getX() + shape.getWidth() + margins.right) {
                maxX = shape.getX() + shape.getWidth() + margins.right; 
            }
            if (maxY < shape.getY() + shape.getHeight() + margins.bottom) {
                maxY = shape.getY() + shape.getHeight() + margins.bottom; 
            }
        }
        
        // iterate over all contained edges and check their bounds
        for (ElkEdge edge : node.getContainedEdges()) {
            for (ElkEdgeSection section : edge.getSections()) {
                double sX = section.getStartX();
                double eX = section.getEndX();
                double sY = section.getStartY();
                double eY = section.getEndY();

                minX = Math.min(minX, sX);
                minX = Math.min(minX, eX);

                maxX = Math.max(maxX, sX);
                maxX = Math.max(maxX, eX);

                minY = Math.min(minY, sY);
                minY = Math.min(minY, eY);

                maxY = Math.max(maxY, sY);
                maxY = Math.max(maxY, eY);
                
                for (ElkBendPoint bendpoint : section.getBendPoints()) {

                    minX = Math.min(minX, bendpoint.getX());
                    maxX = Math.max(maxX, bendpoint.getX());
                    minY = Math.min(minY, bendpoint.getY());
                    maxY = Math.max(maxY, bendpoint.getY());
                    
                }
            }
        }
        // max and min value represent outermost bounds of the layout
        node.setProperty(CoreOptions.CHILD_AREA_WIDTH, maxX - minX);
        node.setProperty(CoreOptions.CHILD_AREA_HEIGHT, maxY - minY);
    }

    /**
     * Determine the junction points of the given edge. This is done by comparing the bend points
     * of the given edge with the bend points of all other edges that are connected to the same
     * source port or the same target port. Note that this method requires all edges to have exactly one
     * edge section.
     * 
     * @param edge an edge
     * @return a list of junction points
     * @throws IllegalArgumentException if the edge has no or more than one edge section.
     */
    public static KVectorChain determineJunctionPoints(final ElkEdge edge) {
        if (edge.getSections().size() != 1) {
            throw new IllegalArgumentException(
                    "The edge needs to have exactly one edge section. Found: " + edge.getSections().size());
        }

        KVectorChain junctionPoints = new KVectorChain();

        if (ElkGraphUtil.connectableShapeToPort(edge.getSources().get(0)) != null) {
            junctionPoints.addAll(determineJunctionPoints(edge,
                    ElkGraphUtil.connectableShapeToPort(edge.getSources().get(0)), false));
        }
        if (ElkGraphUtil.connectableShapeToPort(edge.getTargets().get(0)) != null) {
            junctionPoints.addAll(
                    determineJunctionPoints(edge, ElkGraphUtil.connectableShapeToPort(edge.getTargets().get(0)), true));
        }

        return junctionPoints;
    }

    /**
     * Determine the junction points of the given edge with any edge connected to the given port.
     * This is done by comparing the bend points of the given edge with the bend points of all other edges
     * that are connected to the given port. Note that this method requires all edges to have exactly one
     * edge section.
     * 
     * @param edge an edge
     * @param port one of the ports the edge is connected to
     * @param reverse flag to indicate whether the points are traversed forward or reverse
     * @return a list of junction points
     * @throws IllegalArgumentException if a connected edge has no or more than one edge section.
     */
    private static KVectorChain determineJunctionPoints(final ElkEdge edge, final ElkPort port, final boolean reverse) {
        // Ensure exactly one section
        assert edge.getSections().size() == 1;

        // Grab the edge section
        final ElkEdgeSection section = edge.getSections().get(0);

        // Collection for the junction points of the current edge
        KVectorChain junctionPoints = new KVectorChain();

        // Store the points of the edge in a map for efficiency
        Map<ElkEdgeSection, ArrayList<KVector>> pointsMap = Maps.newHashMap();
        ArrayList<KVector> sectionPoints = getPoints(section);
        pointsMap.put(section, sectionPoints);

        // Store the offset of the other edges
        Map<ElkEdgeSection, KVector> offsetMap = Maps.newHashMap();

        // let allConnectedEdges be the set of edge sections connected to port without the main edge
        List<ElkEdgeSection> allConnectedSections = Lists.newLinkedList();
        for (ElkEdge otherEdge : ElkGraphUtil.allIncidentEdges(port)) {
            if (edge.getSections().size() != 1) {
                throw new IllegalArgumentException(
                        "The edge needs to have exactly one edge section. Found: " + edge.getSections().size());
            }
            if (otherEdge != edge) {
                ElkEdgeSection otherSection = otherEdge.getSections().get(0);
                allConnectedSections.add(otherSection);

                // Edges might have different containments leading to different coordinate systems
                // We can calculate the offset between the edges by comparing the shared port
                ArrayList<KVector> otherPoints = pointsMap.get(otherSection);
                if (otherPoints == null) {
                    otherPoints = getPoints(otherSection);
                    pointsMap.put(otherSection, otherPoints);
                }

                KVector offset = reverse
                        ? new KVector(sectionPoints.get(sectionPoints.size() - 1))
                                .sub(otherPoints.get(otherPoints.size() - 1))
                        : new KVector(sectionPoints.get(0)).sub(otherPoints.get(0));

                offsetMap.put(otherSection, offset);
            }
        }

        if (!allConnectedSections.isEmpty()) {

            // let p1 be the start point
            KVector p1 = sectionPoints.get(reverse ? sectionPoints.size() - 1 : 0);

            // for all bend points of this connection
            for (int i = 1; i < sectionPoints.size(); i++) {
                // let p2 be the next bend point on this connection
                KVector p2 = sectionPoints.get(reverse ? sectionPoints.size() - 1 - i : i);

                // for all other connections that are still on the same track as this one
                Iterator<ElkEdgeSection> allSectIter = allConnectedSections.iterator();
                while (allSectIter.hasNext()) {
                    ElkEdgeSection otherSection = allSectIter.next();
                    ArrayList<KVector> otherPoints = pointsMap.get(otherSection);
                    if (otherPoints.size() <= i) {
                        allSectIter.remove();
                    } else {
                        // let p3 be the next bend point of the other connection
                        KVector p3 = new KVector(otherPoints.get(reverse ? otherPoints.size() - 1 - i : i))
                                .add(offsetMap.get(otherSection));
                        if (p2.x != p3.x || p2.y != p3.y) {
                            // the next point of this and the other connection differ
                            double dx2 = p2.x - p1.x;
                            double dy2 = p2.y - p1.y;
                            double dx3 = p3.x - p1.x;
                            double dy3 = p3.y - p1.y;
                            if ((dx3 * dy2) == (dy3 * dx2) && Math.signum(dx2) == Math.signum(dx3)
                                    && Math.signum(dy2) == Math.signum(dy3)) {

                                // the points p1, p2, p3 form a straight line,
                                // now check whether p2 is between p1 and p3
                                if (Math.abs(dx2) < Math.abs(dx3) || Math.abs(dy2) < Math.abs(dy3)) {
                                    junctionPoints.add(p2);
                                }

                            } else if (i > 1) {
                                // p2 and p3 have diverged, so the last common point is p1
                                junctionPoints.add(p1);
                            }

                            // do not consider the other connection in the next iterations
                            allSectIter.remove();
                        }
                    }
                }
                // for the next iteration p2 is taken as reference point
                p1 = p2;
            }
        }

        return junctionPoints;
    }

    /**
     * Compute the bounding box of the labels of the port corresponding to the {@code port}. This bounding box
     * coordinates are relative to the port (as the location of each port).
     * 
     * @param port
     *            The port of the labels
     * @return the bounds of the labels
     */
    public static ElkRectangle getLabelsBounds(final PortAdapter<?> port) {
        ElkRectangle bounds = null;
        for (LabelAdapter<?> label : port.getLabels()) {
            ElkRectangle currentLabelBounds = new ElkRectangle(label.getPosition().x, label.getPosition().y,
                    label.getSize().x, label.getSize().y);
            if (bounds == null) {
                // First label
                bounds = currentLabelBounds;
            } else {
                bounds.union(currentLabelBounds);
            }
        }
        if (bounds == null) {
            // No label case, return an empty rectangle
            bounds = new ElkRectangle();
        }
        return bounds;
    }

    /**
     * Compute the part of the labels of the port that is inside the node. This method must be used for label with
     * "fixed" port label placement ({@link PortLabelPlacement#isFixed(Set)}).
     * 
     * @param port
     *            The port
     * @param portBorderOffset
     *            The border offset of the port
     * @return the part of the labels that is inside the node
     */
    public static double computeInsidePart(final PortAdapter<?> port, final double portBorderOffset) {
        ElkRectangle labelBounds = getLabelsBounds(port);
        return computeInsidePart(new KVector(labelBounds.x, labelBounds.y),
                new KVector(labelBounds.width, labelBounds.height), port.getSize(), portBorderOffset, port.getSide());
    }

    /**
     * Compute the part of the label of the port that is inside the node. This method must be used for label with
     * "fixed" port label placement ({@link PortLabelPlacement#isFixed(Set)}).
     * 
     * @param labelPosition
     *            The coordinates of the label (relative to the port)
     * @param labelSize
     *            The size of the label
     * @param portSize
     *            The size of the port
     * @param portBorderOffset
     *            The border offset of the port
     * @param portSide
     *            The side of the port
     * @return the part of the label that is inside the node
     */
    public static double computeInsidePart(final KVector labelPosition, final KVector labelSize, final KVector portSize,
            final double portBorderOffset, final PortSide portSide) {
        double insidePart = 0;
        switch (portSide) {
        case NORTH:
            insidePart = Math.max(0, labelSize.y + labelPosition.y - (portSize.y + portBorderOffset));
            break;

        case SOUTH:
            insidePart = Math.max(0, -labelPosition.y - portBorderOffset);
            break;

        case EAST:
            insidePart = Math.max(0, -labelPosition.x - portBorderOffset);
            break;

        case WEST:
            insidePart = Math.max(0, labelSize.x + labelPosition.x - (portSize.x + portBorderOffset));
            break;
        }

        return insidePart;
    }
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // COORDINATE TRANSLATION

    /**
     * Translates the contents of the given node by an offset.
     *
     * @param parent parent node whose children shall be translated
     * @param xoffset x coordinate offset
     * @param yoffset y coordinate offset
     */
    public static void translate(final ElkNode parent, final double xoffset, final double yoffset) {
        for (ElkNode child : parent.getChildren()) {
            // Translate node position
            child.setLocation(child.getX() + xoffset, child.getY() + yoffset);
        }
        // Translates all edges contained in the parent. This includes edges connecting the parent to its
        // children. For these edges the start or end point might get separated from the node boundary.
        parent.getContainedEdges().forEach(edge -> translate(edge, xoffset, yoffset));
    }

    /**
     * Translates the given edge by an offset. This includes all routing information, junction points (if any), and
     * edge labels.
     *
     * @param edge edge that shall be translated
     * @param xoffset x coordinate offset
     * @param yoffset y coordinate offset
     */
    public static void translate(final ElkEdge edge, final double xoffset, final double yoffset) {
        // Edge sections
        edge.getSections().stream().forEach(
                s -> translate(s, xoffset, yoffset));

        // Edge labels
        edge.getLabels().stream().forEach(
                label -> label.setLocation(label.getX() + xoffset, label.getY() + yoffset));

        // Junction points
        KVectorChain junctionPoints = edge.getProperty(CoreOptions.JUNCTION_POINTS);
        if (junctionPoints != null) {
            junctionPoints.offset(xoffset, yoffset);
        }
    }

    /**
     * Translates the given edge section by an offset.
     *
     * @param section edge section that shall be translated
     * @param xoffset x coordinate offset
     * @param yoffset y coordinate offset
     */
    public static void translate(final ElkEdgeSection section, final double xoffset, final double yoffset) {
        // Translate source point
        section.setStartLocation(section.getStartX() + xoffset, section.getStartY() + yoffset);

        // Translate bend points
        for (ElkBendPoint bendPoint : section.getBendPoints()) {
            bendPoint.set(bendPoint.getX() + xoffset, bendPoint.getY() + yoffset);
        }

        // Translate target point
        section.setEndLocation(section.getEndX() + xoffset, section.getEndY() + yoffset);
    }
    
    /**
     * Translates the contents of the given node based on the content alignment property without resizing the node itself.
     * 
     * @param parent The parent node.
     * @param newSize The new size.
     * @param oldSize The old size.
     */
    public static void translate(final ElkNode parent, final KVector newSize, final KVector oldSize) {
        Set<ContentAlignment> contentAlignment =
                parent.getProperty(CoreOptions.CONTENT_ALIGNMENT);
        double xTranslate = 0;
        double yTranslate = 0;
        
        // Horizontal alignment
        if (newSize.x > oldSize.x) {
            if (contentAlignment.contains(ContentAlignment.H_CENTER)) {
                xTranslate = (newSize.x - oldSize.x) / 2f;
            } else if (contentAlignment.contains(ContentAlignment.H_RIGHT)) {
                xTranslate = newSize.x - oldSize.x;
            }
        }

        // Vertical alignment
        if (newSize.y > oldSize.y) {
            if (contentAlignment.contains(ContentAlignment.V_CENTER)) {
                yTranslate = (newSize.y - oldSize.y) / 2f;
            } else if (contentAlignment.contains(ContentAlignment.V_BOTTOM)) {
                yTranslate = newSize.y - oldSize.y;
            }
        }

        translate(parent, xTranslate, yTranslate);
    }


    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // COORDINATE-SYSTEM AND TYPE "CONVERSION"
    
    /**
     * Returns the absolute position of the given element. For nodes and ports, this is exactly what you would expect.
     * Edges don't exactly have an absolute position, so we simply return the absolute position of their containign
     * node. For labels, we walk up the parent relationship and keep adding up positions.
     */
    public static KVector absolutePosition(final ElkGraphElement element) {
        if (element instanceof ElkNode) {
            ElkNode node = (ElkNode) element;
            return toAbsolute(new KVector(node.getX(), node.getY()), node.getParent());
        
        } else if (element instanceof ElkPort) {
            ElkPort port = (ElkPort) element;
            return toAbsolute(new KVector(port.getX(), port.getY()), port.getParent());
            
        } else if (element instanceof ElkEdge) {
            ElkEdge edge = (ElkEdge) element;
            return absolutePosition(edge.getContainingNode());
            
        } else if (element instanceof ElkLabel) {
            ElkLabel label = (ElkLabel) element;
            KVector absoluteParentPosition = absolutePosition(label.getParent());
            return new KVector(absoluteParentPosition.x + label.getX(), absoluteParentPosition.y + label.getY());
        
        } else {
            return null;
        }
    }

    /**
     * Converts the given relative point to an absolute location.
     *
     * @param point a relative point
     * @param parent the parent node to which the point is relative to
     * @return {@code point} for convenience
     */
    public static KVector toAbsolute(final KVector point, final ElkNode parent) {
        ElkNode node = parent;
        while (node != null) {
            point.add(node.getX(), node.getY());
            node = node.getParent();
        }
        return point;
    }

    /**
     * Converts the given absolute point to a relative location.
     *
     * @param point an absolute point
     * @param parent the parent node to which the point shall be made relative to
     * @return {@code point} for convenience
     */
    public static KVector toRelative(final KVector point, final ElkNode parent) {
        ElkNode node = parent;
        while (node != null) {
            point.add(-node.getX(), -node.getY());
            node = node.getParent();
        }
        return point;
    }

    /**
     * Get the edge section points as an array of vectors. 
     * Unnecessary bend points are filtered out automatically. 
     * 
     * @param section an edge section
     * @return an array with all needed edge section points
     */
    private static ArrayList<KVector> getPoints(final ElkEdgeSection section) {
        int n = section.getBendPoints().size() + 2;
        ArrayList<KVector> points = new ArrayList<KVector>(n);

        // Source point
        points.add(new KVector(section.getStartX(), section.getStartY()));

        // Bend points
        section.getBendPoints().stream()
                .forEach(bendPoint -> points.add(new KVector(bendPoint.getX(), bendPoint.getY())));

        // Target point
        points.add(new KVector(section.getEndX(), section.getEndY()));

        // Filter unnecessary bend points from the list
        int i = 1;
        while (i < points.size() - 1) {
            // Unnecessary bend points are given if three points in a row have the same x or y coordinate
            KVector p1 = points.get(i - 1);
            KVector p2 = points.get(i);
            KVector p3 = points.get(i + 1);
            
            if ((p1.x == p2.x && p2.x == p3.x) || (p1.y == p2.y && p2.y == p3.y)) {
                // Found a straight segment, drop p2 and re-check with the same i
                points.remove(i);
            } else {
                // Points are not a straight line, advance
                ++i;
            }
        }
        
        return points;
    }

    /**
     * Creates a vector chain containing the start point, bend points, and end point of the given edge section. Note
     * that modifying the vector chain will be of no consequence to the edge section.
     *
     * @param edgeSection
     *            the edge section to initialize the vector chain with.
     * @return the vector chain.
     */
    public static KVectorChain createVectorChain(final ElkEdgeSection edgeSection) {
        KVectorChain chain = new KVectorChain();

        chain.add(new KVector(edgeSection.getStartX(), edgeSection.getStartY()));

        for (ElkBendPoint bendPoint : edgeSection.getBendPoints()) {
            chain.add(new KVector(bendPoint.getX(), bendPoint.getY()));
        }

        chain.add(new KVector(edgeSection.getEndX(), edgeSection.getEndY()));

        return chain;
    }

    /**
     * Applies the vector chain's vectors to the given edge section. The first and the last point of the vector chain
     * are used as the section's new source and start point, respectively. The remaining points become the section's
     * new bend points. The method tries to reuse as many bend points as possible instead of wiping all bend points
     * out and creating new ones.
     *
     * @param vectorChain the vector chain to apply.
     * @param section the edge section to apply the chain to.
     * @throws IllegalArgumentException if the vector chain contains less than two vectors.
     */
    public static void applyVectorChain(final KVectorChain vectorChain, final ElkEdgeSection section) {
        // We need at least a start and an end point
        if (vectorChain.size() < 2) {
            throw new IllegalArgumentException("The vector chain must contain at least a source and a target point.");
        }

        // Start point
        KVector firstPoint = vectorChain.getFirst();
        section.setStartLocation(firstPoint.x, firstPoint.y);

        // Reuse as many existing bend points as possible
        ListIterator<ElkBendPoint> oldPointIter = section.getBendPoints().listIterator();
        ListIterator<KVector> newPointIter = vectorChain.listIterator(1);

        while (newPointIter.nextIndex() < vectorChain.size() - 1) {
            KVector nextPoint = newPointIter.next();
            ElkBendPoint bendpoint;
            if (oldPointIter.hasNext()) {
                bendpoint = oldPointIter.next();
            } else {
                bendpoint = ElkGraphFactory.eINSTANCE.createElkBendPoint();
                oldPointIter.add(bendpoint);
            }

            bendpoint.set(nextPoint.x, nextPoint.y);
        }

        // Remove existing bend points that we did not use
        while (oldPointIter.hasNext()) {
            oldPointIter.next();
            oldPointIter.remove();
        }

        // End point
        KVector lastPoint = vectorChain.getLast();
        section.setEndLocation(lastPoint.x, lastPoint.y);
    }


    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // DEFAULT LAYOUT SETTINGS

    /**
     * Recursively configures default values for all child elements of the passed graph. This
     * includes node, ports, and edges.
     *
     * Note that it is <b>not</b> checked whether the defaults flag is set on the graph.
     *
     * @param graph
     *            the graph to configure.
     */
    public static void configureDefaultsRecursively(final ElkNode graph) {
        // note that Iterators#filter(Iterator, Class) isn't used on purpose since
        // it's incompatible with gwt
        Iterator<EObject> kgeIt = Iterators.filter(graph.eAllContents(),
            e -> e instanceof ElkGraphElement);
        while (kgeIt.hasNext()) {
            EObject kge = kgeIt.next();
            if (kge instanceof ElkNode) {
                configureWithDefaultValues((ElkNode) kge);
            } else if (kge instanceof ElkPort) {
                configureWithDefaultValues((ElkPort) kge);
            } else if (kge instanceof ElkEdge) {
                configureWithDefaultValues((ElkEdge) kge);
            }
        }
    }

    /**
     * Adds some default values to the passed node. This includes a reasonable size, a label based
     * on the node's identifier and a inside center node label placement.
     *
     * Such default values are useful for fast test case generation.
     *
     * @param node
     *            a node of a graph
     */
    public static void configureWithDefaultValues(final ElkNode node) {
        // Make sure the node has a size if the size constraints are fixed
        Set<SizeConstraint> sc = node.getProperty(CoreOptions.NODE_SIZE_CONSTRAINTS);

        if (sc.equals(SizeConstraint.fixed()) && node.getWidth() == 0f && node.getHeight() == 0f) {
            node.setWidth(DEFAULT_MIN_WIDTH * 2 * 2);
            node.setHeight(DEFAULT_MIN_HEIGHT * 2 * 2);
        }

        // label
        ensureLabel(node);

        Set<NodeLabelPlacement> nlp = node.getProperty(CoreOptions.NODE_LABELS_PLACEMENT);
        if (nlp.equals(NodeLabelPlacement.fixed())) {
            node.setProperty(CoreOptions.NODE_LABELS_PLACEMENT, NodeLabelPlacement.insideCenter());
        }

    }

    /**
     * Adds some default values to the passed port. This includes a reasonable size
     * and a label based on the port's identifier.
     *
     * Such default values are useful for fast test case generation.
     *
     * @param port
     *            a port of a node of a graph
     */
    public static void configureWithDefaultValues(final ElkPort port) {
        if (port.getWidth() == 0 && port.getHeight() == 0) {
            port.setWidth(DEFAULT_MIN_WIDTH / 2 / 2);
            port.setHeight(DEFAULT_MIN_HEIGHT / 2 / 2);
        }

        ensureLabel(port);
    }

    /**
     * Configures the {@link EdgeLabelPlacement} of the passed edge to be center of the edge.
     *
     * @param edge
     *            an edge of a graph
     */
    public static void configureWithDefaultValues(final ElkEdge edge) {
        if (!edge.hasProperty(CoreOptions.EDGE_LABELS_PLACEMENT)) {
            edge.setProperty(CoreOptions.EDGE_LABELS_PLACEMENT, EdgeLabelPlacement.CENTER);
        }
    }

    /**
     * If the element does not already own a label, a label is created based on the element's
     * identifier.
     */
    private static void ensureLabel(final ElkGraphElement klge) {
        if (klge.getLabels().isEmpty()) {
            if (!Strings.isNullOrEmpty(klge.getIdentifier())) {
                ElkLabel label = ElkGraphUtil.createLabel(klge);
                label.setText(klge.getIdentifier());
            }
        }
    }


    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // VISITORS

    /**
     * Apply the given graph element visitors to the content of the given graph. If validators are involved
     * they are not queried about errors.
     *
     * @param graph the graph the visitors shall be applied to.
     * @param visitors the visitors to apply.
     */
    public static void applyVisitors(final ElkNode graph, final IGraphElementVisitor... visitors) {
        Iterator<EObject> allElements = ElkGraphUtil.propertiesSkippingIteratorFor(graph, true);
        while (allElements.hasNext()) {
            EObject nextElement = allElements.next();

            if (nextElement instanceof ElkGraphElement) {
                ElkGraphElement graphElement = (ElkGraphElement) nextElement;
                for (int i = 0; i < visitors.length; i++) {
                    visitors[i].visit(graphElement);
                }
            }
        }
    }

    /**
     * Apply the given graph element visitors to the content of the given graph. If validators are involved
     * and at least one error is found, a {@link GraphValidationException} is thrown.
     *
     * @param graph the graph the visitors shall be applied to.
     * @param visitors the visitors to apply.
     * @throws GraphValidationException if an error is found while validating the graph
     */
    public static void applyVisitorsWithValidation(final ElkNode graph, final IGraphElementVisitor... visitors)
                throws GraphValidationException {

        applyVisitors(graph, visitors);

        // Gather validator results and generate an error message
        List<GraphIssue> allIssues = null;
        for (int i = 0; i < visitors.length; i++) {
            if (visitors[i] instanceof IValidatingGraphElementVisitor) {
                Collection<GraphIssue> issues = ((IValidatingGraphElementVisitor) visitors[i]).getIssues();
                if (!issues.isEmpty()) {
                    if (allIssues == null) {
                        allIssues = new ArrayList<GraphIssue>(issues);
                    } else {
                        allIssues.addAll(issues);
                    }
                }
            }
        }

        // TODO print warnings to log if there are no errors
        if (allIssues != null
                && allIssues.stream().anyMatch(issue -> issue.getSeverity() == GraphIssue.Severity.ERROR)) {
            StringBuilder message = new StringBuilder();
            for (GraphIssue issue : allIssues) {
                if (message.length() > 0) {
                    message.append("\n");
                }
                message.append(issue.getSeverity())
                    .append(": ")
                    .append(issue.getMessage())
                    .append("\n\tat ");
                ElkUtil.printElementPath(issue.getElement(), message);
            }
            throw new GraphValidationException(message.toString(), allIssues);
        }
    }


    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // DEBUGGING

    /**
     * Print information on the given graph element to the given string builder.
     */
    public static void printElementPath(final ElkGraphElement element, final StringBuilder builder) {
        // Print the containing element
        if (element.eContainer() instanceof ElkGraphElement) {
            printElementPath((ElkGraphElement) element.eContainer(), builder);
            builder.append(" > ");
        } else {
            builder.append("Root ");
        }

        // Print the class name
        String className = element.eClass().getName();
        if (className.startsWith("Elk")) {
            // CHECKSTYLEOFF MagicNumber
            builder.append(className.substring(3));
            // CHECKSTYLEON MagicNumber
        } else {
            builder.append(className);
        }

        // Print the identifier if present
        String identifier = element.getIdentifier();
        if (!Strings.isNullOrEmpty(identifier)) {
            builder.append(' ').append(identifier);
            return;
        }
        // Print the label if present
        if (element instanceof ElkLabel) {
            String text = ((ElkLabel) element).getText();
            if (!Strings.isNullOrEmpty(text)) {
                builder.append(' ').append(text);
                return;
            }
        }
        for (ElkLabel label : element.getLabels()) {
            String text = label.getText();
            if (!Strings.isNullOrEmpty(text)) {
                builder.append(' ').append(text);
                return;
            }
        }
        // If it's an edge and no identifier nor label is present, print source and target
        if (element instanceof ElkEdge) {
            ElkEdge edge = (ElkEdge) element;
            if (edge.isConnected()) {
                builder.append(" (");
                ListIterator<ElkConnectableShape> sourceIter = edge.getSources().listIterator();
                while (sourceIter.hasNext()) {
                    if (sourceIter.nextIndex() > 0) {
                        builder.append(", ");
                    }
                    printElementPath(sourceIter.next(), builder);
                }
                builder.append(" -> ");
                ListIterator<ElkConnectableShape> targetIter = edge.getTargets().listIterator();
                while (targetIter.hasNext()) {
                    if (targetIter.nextIndex() > 0) {
                        builder.append(", ");
                    }
                    printElementPath(targetIter.next(), builder);
                }
                builder.append(")");
            }
        }
    }

    /**
     * Returns a path where debug files can be stored. All debug folders end up in an ELK-specific folder placed
     * inside the user's home folder. The returned path can either be that ELK-specific folder itself or a
     * subfolder.
     * <p>
     * If the returned path does not exist, it is not automatically created.
     *
     * @param subfolders optional subfolder names. Can be empty, in which case the ELK-specific subfolder of the user's
     *                   home folder is returned.
     * @return debug folder path, including a trailing separator character. Can return {@code null} if the user's
     *         home folder is not defined.
     */
    public static String debugFolderPath(final String... subfolders) {
        // elkjs-exclude-start
        String userHome = System.getProperty("user.home");
        if (userHome != null) {
            StringBuilder path = new StringBuilder(userHome);

            // Make sure we end the path with a separator character
            if (path.charAt(path.length() - 1) != File.separatorChar) {
                path.append(File.separatorChar);
            }

            // The ELK debug directory
            path.append("elk").append(File.separatorChar);

            // Append the subfolder names, if any
            if (subfolders != null) {
                for (String s : subfolders) {
                    path.append(s).append(File.separatorChar);
                }
            }

            return path.toString();
        }
        // elkjs-exclude-end

        return null;
    }
    
    /**
     * Takes the given name and makes it safe to be used as a file or folder name. To do so, we replace all spaces by
     * underscores and everything that is neither digit not standard character by hyphens.
     * 
     * @param name the name to convert to a proper path name.
     * @return the proper path name.
     */
    // elkjs-exclude-start
    public static String toSafePathName(final String name) {
        // Replace whitespace by _
        Pattern whitespace = Pattern.compile("\\s");
        String nameWithoutWhitespace = whitespace.matcher(name).replaceAll("_");
        
        // Replace everything which isn't a-z, A-Z, 0-9 or _ with -
        Pattern allButAllowedCharacters = Pattern.compile("[^a-zA-Z0-9_]");
        return allButAllowedCharacters.matcher(nameWithoutWhitespace).replaceAll("-");
    }
    // elkjs-exclude-end


    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // MISCELLANEOUS



    /**
     * Create a unique identifier for the given graph element. Note that this identifier
     * is not necessarily universally unique, since it uses the hash code, which
     * usually covers only the range of heap space addresses.
     *
     * @param element a graph element
     */
    public static void createIdentifier(final ElkGraphElement element) {
        element.setIdentifier(Integer.toString(element.hashCode()));
    }

}
