/*******************************************************************************
 * Copyright (c) 2014, 2020 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.layered.graph;

import java.util.Collection;
import java.util.Set;

import org.eclipse.elk.alg.layered.graph.LNode.NodeType;
import org.eclipse.elk.alg.layered.options.EdgeConstraint;
import org.eclipse.elk.alg.layered.options.GraphProperties;
import org.eclipse.elk.alg.layered.options.InLayerConstraint;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.LayerConstraint;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.layered.options.PortType;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.math.KVectorChain;
import org.eclipse.elk.core.options.Alignment;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.options.Direction;
import org.eclipse.elk.core.options.PortConstraints;
import org.eclipse.elk.core.options.PortSide;
import org.eclipse.elk.core.options.SizeConstraint;
import org.eclipse.elk.core.util.IndividualSpacings;
import org.eclipse.elk.graph.properties.IProperty;
import org.eclipse.elk.graph.properties.IPropertyHolder;

/**
 * Utility class for importing graphs into the {@link LGraph} format.
 */
public final class LGraphUtil {
    
    /**
     * Hidden constructor to avoid instantiation.
     */
    private LGraphUtil() { }
    
    /**
     * Create a new array by copying the content of the given node collection.
     * 
     * @param nodes a collection of nodes
     * @return an array of nodes
     */
    public static LNode[] toNodeArray(final Collection<LNode> nodes) {
        return nodes.toArray(new LNode[nodes.size()]);
    }
    
    /**
     * Create a new array by copying the content of the given edge collection.
     * 
     * @param edges a collection of edges
     * @return an array of edges
     */
    public static LEdge[] toEdgeArray(final Collection<LEdge> edges) {
        return edges.toArray(new LEdge[edges.size()]);
    }
    
    /**
     * Create a new array by copying the content of the given port collection.
     * 
     * @param ports a collection of ports
     * @return an array of ports
     */
    public static LPort[] toPortArray(final Collection<LPort> ports) {
        return ports.toArray(new LPort[ports.size()]);
    }
    
    ///////////////////////////////////////////////////////////////////////////////
    // Node Resizing

    /**
     * Resize a node to the given width and height, adjusting port and label positions if needed.
     * 
     * @param node a node
     * @param newSize the new size for the node
     * @param movePorts whether port positions should be adjusted
     * @param moveLabels whether label positions should be adjusted
     */
    public static void resizeNode(final LNode node, final KVector newSize, final boolean movePorts,
            final boolean moveLabels) {
        KVector oldSize = new KVector(node.getSize());
        
        float widthRatio = (float) (newSize.x / oldSize.x);
        float heightRatio = (float) (newSize.y / oldSize.y);
        float widthDiff = (float) (newSize.x - oldSize.x);
        float heightDiff = (float) (newSize.y - oldSize.y);

        // Update port positions
        if (movePorts) {
            boolean fixedPorts =
                    node.getProperty(LayeredOptions.PORT_CONSTRAINTS) == PortConstraints.FIXED_POS;
            
            for (LPort port : node.getPorts()) {
                switch (port.getSide()) {
                case NORTH:
                    if (!fixedPorts) {
                        port.getPosition().x *= widthRatio;
                    }
                    break;
                case EAST:
                    port.getPosition().x += widthDiff;
                    if (!fixedPorts) {
                        port.getPosition().y *= heightRatio;
                    }
                    break;
                case SOUTH:
                    if (!fixedPorts) {
                        port.getPosition().x *= widthRatio;
                    }
                    port.getPosition().y += heightDiff;
                    break;
                case WEST:
                    if (!fixedPorts) {
                        port.getPosition().y *= heightRatio;
                    }
                    break;
                }
            }
        }
        
        // Update label positions
        if (moveLabels) {
            for (LLabel label : node.getLabels()) {
                double midx = label.getPosition().x + label.getSize().x / 2;
                double midy = label.getPosition().y + label.getSize().y / 2;
                double widthPercent = midx / oldSize.x;
                double heightPercent = midy / oldSize.y;
                
                if (widthPercent + heightPercent >= 1) {
                    if (widthPercent - heightPercent > 0 && midy >= 0) {
                        // label is on the right
                        label.getPosition().x += widthDiff;
                        label.getPosition().y += heightDiff * heightPercent;
                    } else if (widthPercent - heightPercent < 0 && midx >= 0) {
                        // label is on the bottom
                        label.getPosition().x += widthDiff * widthPercent;
                        label.getPosition().y += heightDiff;
                    }
                }
            }
        }
        
        // Set the new node size
        node.getSize().x = newSize.x;
        node.getSize().y = newSize.y;
        
        // Set fixed size option for the node: now the size is assumed to stay as determined here
        node.setProperty(LayeredOptions.NODE_SIZE_CONSTRAINTS, SizeConstraint.fixed());
    }
    
    ///////////////////////////////////////////////////////////////////////////////
    // Graph Offsetting

    /**
     * Offsets the given graphs by a given offset without moving their nodes to another graph. The order
     * of the graphs in the collection must not depend on their hash values. Otherwise, subsequent
     * layout calls will most likely produce different results.
     *
     * @param graphs the graph to offset.
     * @param offsetx x coordinate offset.
     * @param offsety y coordinate offset.
     */
    public static void offsetGraphs(final Collection<LGraph> graphs, final double offsetx,
            final double offsety) {

        for (LGraph graph : graphs) {
            offsetGraph(graph, offsetx, offsety);
        }
    }
    
    /**
     * Offsets the given graph by a given offset without moving its nodes to another graph. This
     * method can be called as many times as required on a given graph: it does not take the graph's
     * offset into account.
     *
     * @param graph the graph to offset.
     * @param offsetx x coordinate offset.
     * @param offsety y coordinate offset.
     */
    public static void offsetGraph(final LGraph graph, final double offsetx, final double offsety) {
        KVector graphOffset = new KVector(offsetx, offsety);

        for (LNode node : graph.getLayerlessNodes()) {
            node.getPosition().add(graphOffset);
            for (LPort port : node.getPorts()) {
                for (LEdge edge : port.getOutgoingEdges()) {
                    edge.getBendPoints().offset(graphOffset);
                    KVectorChain junctionPoints = edge.getProperty(LayeredOptions.JUNCTION_POINTS);
                    if (junctionPoints != null) {
                        junctionPoints.offset(graphOffset);
                    }
                    for (LLabel label : edge.getLabels()) {
                        label.getPosition().add(graphOffset);
                    }
                }
            }
        }
    }
    
    ///////////////////////////////////////////////////////////////////////////////
    // Layer Things
    
    /**
     * Determines a horizontal placement for all nodes of a layer. The size of the layer is assumed
     * to be already set to the maximal width of the contained nodes (usually done during node
     * placement).
     * 
     * @param layer the layer in which to place the nodes
     * @param xoffset horizontal offset for layer placement
     */
    public static void placeNodesHorizontally(final Layer layer, final double xoffset) {
        // determine maximal left and right margin
        double maxLeftMargin = 0, maxRightMargin = 0;
        for (LNode node : layer.getNodes()) {
            maxLeftMargin = Math.max(maxLeftMargin, node.getMargin().left);
            maxRightMargin = Math.max(maxRightMargin, node.getMargin().right);
        }

        // CHECKSTYLEOFF MagicNumber
        for (LNode node : layer.getNodes()) {
            Alignment alignment = node.getProperty(LayeredOptions.ALIGNMENT);
            double ratio;
            switch (alignment) {
            case LEFT:
                ratio = 0.0;
                break;
            case RIGHT:
                ratio = 1.0;
                break;
            case CENTER:
                ratio = 0.5;
                break;
            default:
                // determine the number of input and output ports for the node
                int inports = 0, outports = 0;
                for (LPort port : node.getPorts()) {
                    if (!port.getIncomingEdges().isEmpty()) {
                        inports++;
                    }
                    
                    if (!port.getOutgoingEdges().isEmpty()) {
                        outports++;
                    }
                }
                
                // calculate node placement based on the port numbers
                if (inports + outports == 0) {
                    ratio = 0.5;
                } else {
                    ratio = (double) outports / (inports + outports);
                }
            }
            
            // align nodes to the layer's maximal margin
            KVector size = layer.getSize();
            double nodeSize = node.getSize().x;
            double xpos = (size.x - nodeSize) * ratio;
            if (ratio > 0.5) {
                xpos -= maxRightMargin * 2 * (ratio - 0.5);
            } else if (ratio < 0.5) {
                xpos += maxLeftMargin * 2 * (0.5 - ratio);
            }
            
            // consider the node's individual margin
            double leftMargin = node.getMargin().left;
            if (xpos < leftMargin) {
                xpos = leftMargin;
            }
            double rightMargin = node.getMargin().right;
            if (xpos > size.x - rightMargin - nodeSize) {
                xpos = size.x - rightMargin - nodeSize;
            }
            
            node.getPosition().x = xoffset + xpos;
        }
    }
    
    /**
     * Finds the maximum width of non-dummy nodes in the given layer.
     * <p>
     * If the graph is laid out in a vertical direction, the maximum non-dummy node width doesn't
     * mean anything since the labels are narrow, but very high. So in that case, {@code 0.0} is
     * returned.
     * </p>
     * <p>
     * When calling the function prior to node margin calculation, the node margins are
     * unknown/invalid. One may set the {@code respectNodeMargins} flag accordingly.It would seem
     * that the solution is to execute label management after node margin calculation. But there's a
     * chicken-and-egg problem right there: if we execute node margin calculation first, it won't
     * know about shortened edge end labels and port labels. If we execute label management first,
     * it won't include node margins in the space usable for center edge labels. Damn!
     * </p>
     *
     * @param layer
     *            the layer to iterate over.
     * @param respectNodeMargins
     *            whether to include node margins in width calculation.
     * @return the maximum width of non-dummy nodes. If there are none or the layout direction is
     *         vertical, {@code 0.0} is returned.
     */
    public static double findMaxNonDummyNodeWidth(final Layer layer,
            final boolean respectNodeMargins) {

        if (layer.getGraph().getProperty(LayeredOptions.DIRECTION).isVertical()) {
            return 0.0;
        }

        double maxWidth = 0.0;

        for (LNode node : layer) {
            if (node.getType() == NodeType.NORMAL) {
                double width = node.getSize().x;
                if (respectNodeMargins) {
                    width += node.getMargin().left + node.getMargin().right;
                }
                maxWidth = Math.max(maxWidth, width);
            }
        }

        return maxWidth;
    }
    
    ///////////////////////////////////////////////////////////////////////////////
    // Graph Properties
    
    /**
     * Compute the graph properties of the given layered graph. These properties are important
     * to determine which intermediate processors are included in the layout run.
     * Ideally the properties are computed during the import of the source format into {@link LGraph},
     * e.g. as done in {@code KGraphImporter}. This method is offered only for convenience.
     * <p>
     * The nodes are expected to be in the {@link LGraph#getLayerlessNodes()} list.
     * </p>
     * 
     * @param layeredGraph a layered graph
     */
    public static void computeGraphProperties(final LGraph layeredGraph) {
        Set<GraphProperties> props = layeredGraph.getProperty(InternalProperties.GRAPH_PROPERTIES);
        if (!props.isEmpty()) {
            props.clear();
        }
        
        Direction direction = getDirection(layeredGraph);
        for (LNode node : layeredGraph.getLayerlessNodes()) {
            if (node.getProperty(LayeredOptions.COMMENT_BOX)) {
                props.add(GraphProperties.COMMENTS);
            } else if (node.getProperty(LayeredOptions.HYPERNODE)) {
                props.add(GraphProperties.HYPERNODES);
                props.add(GraphProperties.HYPEREDGES);
            } else if (node.getType() == NodeType.EXTERNAL_PORT) {
                props.add(GraphProperties.EXTERNAL_PORTS);
            }
            
            PortConstraints portConstraints = node.getProperty(LayeredOptions.PORT_CONSTRAINTS);
            if (portConstraints == PortConstraints.UNDEFINED) {
                // correct the port constraints value
                node.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FREE);
            } else if (portConstraints != PortConstraints.FREE) {
                props.add(GraphProperties.NON_FREE_PORTS);
            }
            
            for (LPort port : node.getPorts()) {
                if (port.getIncomingEdges().size() + port.getOutgoingEdges().size() > 1) {
                    props.add(GraphProperties.HYPEREDGES);
                }
                
                PortSide portSide = port.getSide();
                switch (direction) {
                case UP:
                case DOWN:
                    if (portSide == PortSide.EAST || portSide == PortSide.WEST) {
                        props.add(GraphProperties.NORTH_SOUTH_PORTS);
                    }
                    break;
                default:
                    if (portSide == PortSide.NORTH || portSide == PortSide.SOUTH) {
                        props.add(GraphProperties.NORTH_SOUTH_PORTS);
                    }
                }
                
                for (LEdge edge : port.getOutgoingEdges()) {
                    if (edge.getTarget().getNode() == node) {
                        props.add(GraphProperties.SELF_LOOPS);
                    }
                    
                    for (LLabel label : edge.getLabels()) {
                        switch (label.getProperty(LayeredOptions.EDGE_LABELS_PLACEMENT)) {
                        case CENTER:
                            props.add(GraphProperties.CENTER_LABELS);
                            break;
                        case HEAD:
                        case TAIL:
                            props.add(GraphProperties.END_LABELS);
                            break;
                        }
                    }
                }
            }
        }
    }
    
    
    ///////////////////////////////////////////////////////////////////////////////
    // Handling of Ports

    /**
     * Create a port for an edge that is not connected to a port. This is necessary because ELK
     * Layered wants all edges to have a source port and a target port. The port side is computed
     * from the given absolute end point position of the edge.
     * 
     * @param node
     *            the node at which the edge is incident
     * @param endPoint
     *            the absolute point where the edge ends, or {@code null} if unknown
     * @param type
     *            the port type
     * @param layeredGraph
     *            the layered graph
     * @return a new port
     */
    public static LPort createPort(final LNode node, final KVector endPoint, final PortType type,
            final LGraph layeredGraph) {
        LPort port;
        Direction direction = getDirection(layeredGraph);
        boolean mergePorts = layeredGraph.getProperty(LayeredOptions.MERGE_EDGES);
        
        if ((mergePorts || node.getProperty(LayeredOptions.HYPERNODE))
                && !node.getProperty(LayeredOptions.PORT_CONSTRAINTS).isSideFixed()) {
            
            // Hypernodes have one output port and one input port
            PortSide defaultSide = PortSide.fromDirection(direction);
            port = provideCollectorPort(layeredGraph, node, type,
                    type == PortType.OUTPUT ? defaultSide : defaultSide.opposed());
        } else {
            port = new LPort();
            port.setNode(node);
            
            if (endPoint != null) {
                KVector pos = port.getPosition();
                pos.x = endPoint.x - node.getPosition().x;
                pos.y = endPoint.y - node.getPosition().y;
                pos.bound(0, 0, node.getSize().x, node.getSize().y);
                port.setSide(calcPortSide(port, direction));
            } else {
                PortSide defaultSide = PortSide.fromDirection(direction);
                port.setSide(type == PortType.OUTPUT ? defaultSide : defaultSide.opposed());
            }
            
            Set<GraphProperties> graphProperties = layeredGraph.getProperty(
                    InternalProperties.GRAPH_PROPERTIES);
            PortSide portSide = port.getSide();
            switch (direction) {
            case LEFT:
            case RIGHT:
                if (portSide == PortSide.NORTH || portSide == PortSide.SOUTH) {
                    graphProperties.add(GraphProperties.NORTH_SOUTH_PORTS);
                }
                break;
            case UP:
            case DOWN:
                if (portSide == PortSide.EAST || portSide == PortSide.WEST) {
                    graphProperties.add(GraphProperties.NORTH_SOUTH_PORTS);
                }
                break;
            }
        }
        
        return port;
    }
    
    /**
     * Determine the port side for the given port from its relative position at
     * its corresponding node.
     * 
     * @param port port to analyze
     * @param direction the overall layout direction
     * @return the port side relative to its containing node
     */
    static PortSide calcPortSide(final LPort port, final Direction direction) {
        LNode node = port.getNode();
        // if the node has zero size, we cannot decide anything
        double nodeWidth = node.getSize().x;
        double nodeHeight = node.getSize().y;
        if (nodeWidth <= 0 && nodeHeight <= 0) {
            return PortSide.UNDEFINED;
        }

        // check direction-dependent criterion
        double xpos = port.getPosition().x;
        double ypos = port.getPosition().y;
        double width = port.getSize().x;
        double height = port.getSize().y;
        switch (direction) {
        case LEFT:
        case RIGHT:
            if (xpos < 0) {
                return PortSide.WEST;
            } else if (xpos + width > nodeWidth) {
                return PortSide.EAST;
            }
            break;
        case UP:
        case DOWN:
            if (ypos < 0) {
                return PortSide.NORTH;
            } else if (ypos + height > nodeHeight) {
                return PortSide.SOUTH;
            }
        }
        
        // check general criterion
        double widthPercent = (xpos + width / 2) / nodeWidth;
        double heightPercent = (ypos + height / 2) / nodeHeight;
        if (widthPercent + heightPercent <= 1
                && widthPercent - heightPercent <= 0) {
            // port is on the left
            return PortSide.WEST;
        } else if (widthPercent + heightPercent >= 1
                && widthPercent - heightPercent >= 0) {
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
     * Compute the offset for a port, that is the amount by which it is moved outside of the node.
     * An offset value of 0 means the port has no intersection with the node and touches the outside
     * border of the node.
     * 
     * @param port a port
     * @param side the side on the node for the given port
     * @return the offset on the side
     */
    static double calcPortOffset(final LPort port, final PortSide side) {
        LNode node = port.getNode();
        switch (side) {
        case NORTH:
            return -(port.getPosition().y + port.getSize().y);
        case EAST:
            return port.getPosition().x - node.getSize().x;
        case SOUTH:
            return port.getPosition().y - node.getSize().y;
        case WEST:
            return -(port.getPosition().x + port.getSize().x);
        }
        return 0;
    }

    /**
     * Center the given point on one side of a boundary.
     * 
     * @param point
     *            a point to change
     * @param boundary
     *            the boundary to use for centering
     * @param side
     *            the side of the boundary
     */
    static void centerPoint(final KVector point, final KVector boundary,
            final PortSide side) {
        
        switch (side) {
        case NORTH:
            point.x = boundary.x / 2;
            point.y = 0;
            break;
        case EAST:
            point.x = boundary.x;
            point.y = boundary.y / 2;
            break;
        case SOUTH:
            point.x = boundary.x / 2;
            point.y = boundary.y;
            break;
        case WEST:
            point.x = 0;
            point.y = boundary.y / 2;
            break;
        }
    }

    /**
     * Return a collector port of given type, creating it if necessary. A collector port is used to
     * merge all incident edges that originally had no ports.
     * 
     * @param layeredGraph
     *            the layered graph
     * @param node
     *            a node
     * @param type
     *            if {@code INPUT}, an input collector port is returned; if {@code OUTPUT}, an
     *            output collector port is returned
     * @param side
     *            the side to set for a newly created port
     * @return a collector port
     */
    public static LPort provideCollectorPort(final LGraph layeredGraph,
            final LNode node, final PortType type, final PortSide side) {
        
        LPort port = null;
        switch (type) {
        case INPUT:
            for (LPort inport : node.getPorts()) {
                if (inport.getProperty(InternalProperties.INPUT_COLLECT)) {
                    return inport;
                }
            }
            port = new LPort();
            port.setProperty(InternalProperties.INPUT_COLLECT, true);
            break;
        case OUTPUT:
            for (LPort outport : node.getPorts()) {
                if (outport.getProperty(InternalProperties.OUTPUT_COLLECT)) {
                    return outport;
                }
            }
            port = new LPort();
            port.setProperty(InternalProperties.OUTPUT_COLLECT, true);
            break;
        }
        if (port != null) {
            port.setNode(node);
            port.setSide(side);
            centerPoint(port.getPosition(), node.getSize(), side);
        }
        return port;
    }
    
    /**
     * Initialize the side, offset, and anchor point of the given port. The port is assumed to
     * be also present in the original graph structure. The port's current position and the size
     * of the corresponding node are used to determine missing values.
     * 
     * @param port a port
     * @param portConstraints the port constraints of the containing node
     * @param direction the overall layout direction
     * @param anchorPos the anchor position, or {@code null} if it shall be determined automatically
     */
    public static void initializePort(final LPort port, final PortConstraints portConstraints,
            final Direction direction, final KVector anchorPos) {
        
        PortSide portSide = port.getSide();
        
        if (portSide == PortSide.UNDEFINED && portConstraints.isSideFixed()) {
            // calculate the port side and offset from the port's current position
            portSide = calcPortSide(port, direction);
            port.setSide(portSide);
            
            // if port coordinates are (0,0), we default to port offset 0 to make the common case
            // frustration-free
            if (!port.getAllProperties().containsKey(LayeredOptions.PORT_BORDER_OFFSET)
                    && portSide != PortSide.UNDEFINED
                    && (port.getPosition().x != 0 || port.getPosition().y != 0)) {
                
                port.setProperty(LayeredOptions.PORT_BORDER_OFFSET, calcPortOffset(port, portSide));
            }
        }
        
        // if the port constraints are set to fixed ratio, remember the current ratio
        if (portConstraints.isRatioFixed()) {
            double ratio = 0.0;
            
            switch (portSide) {
            case NORTH:
            case SOUTH:
                double nodeWidth = port.getNode().getSize().x;
                if (nodeWidth > 0) {
                    ratio = port.getPosition().x / nodeWidth;
                }
                break;
            case EAST:
            case WEST:
                double nodeHeight = port.getNode().getSize().y;
                if (nodeHeight > 0) {
                    ratio = port.getPosition().y / nodeHeight;
                }
                break;
            }
            
            port.setProperty(InternalProperties.PORT_RATIO_OR_POSITION, ratio);
        }

        KVector portSize = port.getSize();
        KVector portAnchor = port.getAnchor();
        
        // if the port anchor property is set, use it as anchor point
        if (anchorPos != null) {
            portAnchor.x = anchorPos.x;
            portAnchor.y = anchorPos.y;
            
            // Since we have applied an explicit anchor, assume the user knows what they are doing and fix it
            port.setExplicitlySuppliedPortAnchor(true);
        } else if (portConstraints.isSideFixed() && portSide != PortSide.UNDEFINED) {
            // set the anchor point according to the port side
            switch (portSide) {
            case NORTH:
                portAnchor.x = portSize.x / 2;
                break;
            case EAST:
                portAnchor.x = portSize.x;
                portAnchor.y = portSize.y / 2;
                break;
            case SOUTH:
                portAnchor.x = portSize.x / 2;
                portAnchor.y = portSize.y;
                break;
            case WEST:
                portAnchor.y = portSize.y / 2;
                break;
            }
        } else {
            // the port side will be decided later, so set the anchor to the center point
            portAnchor.x = portSize.x / 2;
            portAnchor.y = portSize.y / 2;
        }
    }
    
    
    ///////////////////////////////////////////////////////////////////////////////
    // External Ports  (Ports on the boundary of the parent node)
    
    /**
     * Creates a dummy for an external port. The dummy will have just one port. The port is on
     * the eastern side for western external ports, on the western side for eastern external ports,
     * on the southern side for northern external ports, and on the northern side for southern
     * external ports.
     * 
     * <p>The returned dummy node is decorated with some properties:</p>
     * <ul>
     *   <li>Its node type is set to {@link LNode.NodeType#EXTERNAL_PORT}.</li>
     *   <li>Its {@link InternalProperties#ORIGIN} is set to the external port object.</li>
     *   <li>The {@link LayeredOptions#PORT_CONSTRAINTS} are set to
     *     {@link PortConstraints#FIXED_POS}.</li>
     *   <li>For western and eastern port dummies, the {@link LayeredOptions#LAYER_CONSTRAINT} is set to
     *     {@link LayerConstraint#FIRST_SEPARATE} and {@link LayerConstraint#LAST_SEPARATE},
     *     respectively.</li>
     *   <li>For northern and southern port dummies, the {@link InternalProperties#IN_LAYER_CONSTRAINT}
     *     is set to {@link InLayerConstraint#TOP} and {@link InLayerConstraint#BOTTOM},
     *     respectively.</li>
     *   <li>For eastern dummies, the {@link InternalProperties#EDGE_CONSTRAINT} is set to
     *     {@link EdgeConstraint#OUTGOING_ONLY}; for western dummies, it is set to
     *     {@link EdgeConstraint#INCOMING_ONLY}; for all other dummies, it is left unset.</li>
     *   <li>{@link LayeredOptions#EXT_PORT_SIDE} is set to the side of the external port represented.</li>
     *   <li>If the port constraints of the original port's node are set to
     *     {@link PortConstraints#FIXED_RATIO} or {@link PortConstraints#FIXED_POS}, the dummy node's
     *     {@link InternalProperties#PORT_RATIO_OR_POSITION} property is set to the port's original
     *     position, defined relative to the original node's origin. (as opposed to relative to the
     *     node's content area)</li>
     *   <li>The {@link LayeredOptions#EXT_PORT_SIZE} property is set to the size of the external port the
     *     the dummy represents, while the size of the dummy itself is set to {@code (0, 0)}.</li>
     * </ul>
     * 
     * <p>The layout direction of a graph has implications on the side external ports are placed at.
     * If port constraints imply fixed sides for ports, the {@link LayeredOptions#EXT_PORT_SIDE} property is
     * set to whatever the external port's port side is. If the port side needs to be determined, it
     * depends on the port type (input port or output port; determined by the number of incoming and
     * outgoing edges) and on the layout direction as follows:</p>
     * <table>
     *   <tr><th></th>     <th>Input port</th><th>Output port</th></tr>
     *   <tr><th>Right</th><td>WEST</td>      <td>EAST</td></tr>
     *   <tr><th>Left</th> <td>EAST</td>      <td>WEST</td></tr>
     *   <tr><th>Down</th> <td>NORTH</td>     <td>SOUTH</td></tr>
     *   <tr><th>Up</th>   <td>SOUTH</td>     <td>NORTH</td></tr>
     * </table>
     * 
     * @param propertyHolder property holder for layout options that are set on the original port
     * @param portConstraints constraints for external ports.
     * @param portSide the side of the external port.
     * @param netFlow the number of incoming minus the number of outgoing edges.
     * @param portNodeSize the size of the node the port belongs to. Only relevant if the port
     *                     constraints are {@code FIXED_RATIO}.
     * @param portPosition the current port position. Only relevant if the port constraints are
     *                     {@code FIXED_ORDER}, {@code FIXED_RATIO} or {@code FIXED_POSITION}.
     * @param portSize size of the port. Depending on the port's side, the created dummy will
     *                 have the same width or height as the port, with the other dimension set
     *                 to zero.
     * @param layoutDirection layout direction of the node that owns the port. Must not be
     *                        {@link Direction#UNDEFINED}.
     * @param layeredGraph the layered graph
     * @return a dummy node representing the external port.
     */
    public static LNode createExternalPortDummy(final IPropertyHolder propertyHolder,
            final PortConstraints portConstraints, final PortSide portSide, final int netFlow,
            final KVector portNodeSize, final KVector portPosition, final KVector portSize,
            final Direction layoutDirection, final LGraph layeredGraph) {
        
        PortSide finalExternalPortSide = portSide;
        
        // Create the dummy with one port
        LNode dummy = new LNode(layeredGraph);
        dummy.setType(NodeType.EXTERNAL_PORT);
        dummy.setProperty(InternalProperties.EXT_PORT_SIZE, portSize);
        dummy.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_POS);
        double portBorderOffset = propertyHolder.getProperty(LayeredOptions.PORT_BORDER_OFFSET);
        dummy.setProperty(LayeredOptions.PORT_BORDER_OFFSET, portBorderOffset);
        
        LPort dummyPort = new LPort();
        dummyPort.setNode(dummy);
        
        // If the port constraints are free, we need to determine where to put the dummy (and its port)
        if (!portConstraints.isSideFixed()) {
            // we need some layout direction here, use RIGHT as default in case it is undefined
            assert layoutDirection != Direction.UNDEFINED;
            if (netFlow >= 0) {
                finalExternalPortSide = PortSide.fromDirection(layoutDirection);
            } else {
                finalExternalPortSide = PortSide.fromDirection(layoutDirection).opposed();
            }
            propertyHolder.setProperty(LayeredOptions.PORT_SIDE, finalExternalPortSide);
        }
        
        // Retrieve the anchor point, possibly to be modified later
        KVector anchor = new KVector();
        boolean explicitAnchor = false;
        
        if (propertyHolder.hasProperty(LayeredOptions.PORT_ANCHOR)) {
            anchor.set(propertyHolder.getProperty(LayeredOptions.PORT_ANCHOR));
            explicitAnchor = true;
            
        } else {
            anchor.set(portSize.x / 2, portSize.y / 2);
        }
        
        // With the port side at hand, set the necessary properties and place the dummy's port
        // at the dummy's center
        switch (finalExternalPortSide) {
        case WEST:
            dummy.setProperty(LayeredOptions.LAYERING_LAYER_CONSTRAINT, LayerConstraint.FIRST_SEPARATE);
            dummy.setProperty(InternalProperties.EDGE_CONSTRAINT, EdgeConstraint.OUTGOING_ONLY);
            dummy.getSize().y = portSize.y;
            if (portBorderOffset < 0) {
                dummy.getSize().x = -portBorderOffset;
            }
            dummyPort.setSide(PortSide.EAST);
            if (!explicitAnchor) {
                anchor.x = portSize.x;
            }
            
            // The port anchors think that there is a difference between the port's left and right border
            // coordinates, which makes sense if the port has a non-zero width. The port dummy, however,
            // will have a width of zero. Thus, the anchor must be relative to -portWidth. This fixes #546.
            anchor.x -= portSize.x;
            
            break;
        
        case EAST:
            dummy.setProperty(LayeredOptions.LAYERING_LAYER_CONSTRAINT, LayerConstraint.LAST_SEPARATE);
            dummy.setProperty(InternalProperties.EDGE_CONSTRAINT, EdgeConstraint.INCOMING_ONLY);
            dummy.getSize().y = portSize.y;
            if (portBorderOffset < 0) {
                dummy.getSize().x = -portBorderOffset;
            }
            dummyPort.setSide(PortSide.WEST);
            if (!explicitAnchor) {
                anchor.x = 0;
            }
            break;
        
        case NORTH:
            dummy.setProperty(InternalProperties.IN_LAYER_CONSTRAINT, InLayerConstraint.TOP);
            dummy.getSize().x = portSize.x;
            if (portBorderOffset < 0) {
                dummy.getSize().y = -portBorderOffset;
            }
            dummyPort.setSide(PortSide.SOUTH);
            if (!explicitAnchor) {
                anchor.y = portSize.y;
            }
            
            // See comments in case WEST. This partly fixes #680.
            anchor.y -= portSize.y;
            
            break;
        
        case SOUTH:
            dummy.setProperty(InternalProperties.IN_LAYER_CONSTRAINT, InLayerConstraint.BOTTOM);
            dummy.getSize().x = portSize.x;
            if (portBorderOffset < 0) {
                dummy.getSize().y = -portBorderOffset;
            }
            dummyPort.setSide(PortSide.NORTH);
            if (!explicitAnchor) {
                anchor.y = 0;
            }
            break;
        
        default:
            // Should never happen!
            assert false : finalExternalPortSide;
        }
        
        // Finally apply the anchor by setting the dummy port position accordingly. Also, remember the anchor on the
        // dummy itself since the hierarchical port processors depend on that
        dummyPort.getPosition().set(anchor);
        dummy.setProperty(LayeredOptions.PORT_ANCHOR, anchor);
        
        if (portConstraints.isOrderFixed()) {
            // The order of ports is fixed in some way, so what we will have to do is to remember information about it
            double informationAboutIt = 0;
            
            // If only the order is fixed _and_ the port has an explicit index set on it, remember that
            if (portConstraints == PortConstraints.FIXED_ORDER
                    && propertyHolder.hasProperty(LayeredOptions.PORT_INDEX)) {
                
                // We will have to be careful: on the SOUTH and WEST sides, the index is in reverse to what we would
                // later expect in the code, so we'll use the index * -1 there
                switch (finalExternalPortSide) {
                case NORTH:
                case EAST:
                    informationAboutIt = propertyHolder.getProperty(LayeredOptions.PORT_INDEX).doubleValue();
                    break;
                    
                case SOUTH:
                case WEST:
                    informationAboutIt = -1 * propertyHolder.getProperty(LayeredOptions.PORT_INDEX).doubleValue();
                    break;
                }
                
            } else {
                // Otherwise, we will just go with the position itself
                switch (finalExternalPortSide) {
                case WEST:
                case EAST:
                    informationAboutIt = portPosition.y;
                    if (portConstraints.isRatioFixed()) {
                        informationAboutIt /= portNodeSize.y;
                    }
                    
                    break;
                    
                case NORTH:
                case SOUTH:
                    informationAboutIt = portPosition.x;
                    if (portConstraints.isRatioFixed()) {
                        informationAboutIt /= portNodeSize.x;
                    }
                    
                    break;
                }
            }
            
            dummy.setProperty(InternalProperties.PORT_RATIO_OR_POSITION, informationAboutIt);
        }
        
        // Set the port side of the dummy
        dummy.setProperty(InternalProperties.EXT_PORT_SIDE, finalExternalPortSide);
        
        return dummy;
    }
    
    /**
     * Calculates the position of the external port's top left corner from the position of the
     * given dummy node that represents the port. The position is relative to the graph node's
     * top left corner.
     * 
     * @param graph the graph for which ports shall be placed
     * @param portDummy the dummy node representing the external port.
     * @param portWidth the external port's width.
     * @param portHeight the external port's height.
     * @return the external port's position.
     */
    public static KVector getExternalPortPosition(final LGraph graph, final LNode portDummy,
            final double portWidth, final double portHeight) {
        
        KVector portPosition = new KVector(portDummy.getPosition());
        portPosition.x += portDummy.getSize().x / 2.0;
        portPosition.y += portDummy.getSize().y / 2.0;
        double portOffset = portDummy.getProperty(LayeredOptions.PORT_BORDER_OFFSET);
        
        // Get some properties of the graph
        KVector graphSize = graph.getSize();
        LPadding padding = graph.getPadding();
        KVector graphOffset = graph.getOffset();
        
        // The exact coordinates depend on the port's side... (often enough, these calculations will
        // give the same results as the port coordinates already computed, but depending on the kind of
        // connected components processing, the computed coordinates might be wrong now)
        switch (portDummy.getProperty(InternalProperties.EXT_PORT_SIDE)) {
        case NORTH:
            portPosition.x += padding.left + graphOffset.x - (portWidth / 2.0);
            portPosition.y = -portHeight - portOffset;
            portDummy.getPosition().y = -(padding.top + portOffset + graphOffset.y);
            break;
        
        case EAST:
            portPosition.x = graphSize.x + padding.left + padding.right + portOffset;
            portPosition.y += padding.top + graphOffset.y - (portHeight / 2.0);
            portDummy.getPosition().x = graphSize.x + padding.right + portOffset - graphOffset.x;
            break;
        
        case SOUTH:
            portPosition.x += padding.left + graphOffset.x - (portWidth / 2.0);
            portPosition.y = graphSize.y + padding.top + padding.bottom + portOffset;
            portDummy.getPosition().y = graphSize.y + padding.bottom + portOffset - graphOffset.y;
            break;
        
        case WEST:
            portPosition.x = -portWidth - portOffset;
            portPosition.y += padding.top + graphOffset.y - (portHeight / 2.0);
            portDummy.getPosition().x = -(padding.left + portOffset + graphOffset.x);
            break;
        }
        
        return portPosition;
    }
    
    
    ///////////////////////////////////////////////////////////////////////////////
    // Compound Graphs
    
    /**
     * Determines whether the given child node is a descendant of the parent node.
     * 
     * @param child a child node
     * @param parent a parent node
     * @return true if {@code child} is a direct or indirect child of {@code parent}
     */
    public static boolean isDescendant(final LNode child, final LNode parent) {
        LNode current = child;
        LNode next = current.getGraph().getParentNode();
        while (next != null) {
            current = next;
            if (current == parent) {
                return true;
            }
            next = current.getGraph().getParentNode();
        }
        return false;
    }
    
    /**
     * Converts the given point from the coordinate system of {@code oldGraph} to that of
     * {@code newGraph}. Padding and graph offset are included in the calculation.
     * If the old and new graph are identical, no calculations are made.
     * 
     * @param point a relative point
     * @param oldGraph the graph to which the point is relative to
     * @param newGraph the graph to which the point is made relative to after applying this method
     */
    public static void changeCoordSystem(final KVector point, final LGraph oldGraph,
            final LGraph newGraph) {

        if (oldGraph == newGraph) {
            // nothing has to be done
            return;
        }

        // transform to absolute coordinates
        LGraph graph = oldGraph;
        LNode node;
        do {
            point.add(graph.getOffset());
            node = graph.getParentNode();
            if (node != null) {
                LPadding padding = graph.getPadding();
                point.add(padding.left, padding.top);
                point.add(node.getPosition());
                graph = node.getGraph();
            }
        } while (node != null);
        
        // transform to relative coordinates (to newGraph)
        graph = newGraph;
        do {
            point.sub(graph.getOffset());
            node = graph.getParentNode();
            if (node != null) {
                LPadding padding = graph.getPadding();
                point.sub(padding.left, padding.top);
                point.sub(node.getPosition());
                graph = node.getGraph();
            }
        } while (node != null);
    }
    
    
    ///////////////////////////////////////////////////////////////////////////////
    // Other Stuff
    
    /**
     * {@code LGraph}-local implementation of
     * {@link IndividualSpacings#getIndividualOrInherited(org.eclipse.elk.graph.ElkNode, IProperty)}.
     * 
     * @param node
     *            the node whose property value to return. If the property is not set as part of an
     *            {@link IndividualSpacings} object, we use the property value set on the graph the node is part of.
     * @param property
     *            the property whose value to return.
     * @return the property value.
     */
    public static <T> T getIndividualOrInherited(final LNode node, final IProperty<T> property) {
        T result = null;
        
        if (node.hasProperty(CoreOptions.SPACING_INDIVIDUAL)) {
            IPropertyHolder individualSpacings = node.getProperty(CoreOptions.SPACING_INDIVIDUAL);
            if (individualSpacings.hasProperty(property)) {
                result = individualSpacings.getProperty(property);
            }
        }
        
        // use the common value
        if (result == null && node.getGraph() != null) {
            result = node.getGraph().getProperty(property);
        }
        
        return result;
    }
    
    /**
     * Determine the layout direction for the given graph. If the direction option is undefined,
     * a suitable direction is chosen depending on the aspect ratio.
     * 
     * @param graph a layered graph
     * @return the layout direction to apply for the graph
     */
    public static Direction getDirection(final LGraph graph) {
        Direction direction = graph.getProperty(LayeredOptions.DIRECTION);
        if (direction == Direction.UNDEFINED) {
            double aspectRatio = graph.getProperty(LayeredOptions.ASPECT_RATIO);
            if (aspectRatio >= 1) {
                // Default due to default value of the ASPECT_RATIO option for the algorithm
                return Direction.RIGHT;
            } else {
                return Direction.DOWN;
            }
        }
        return direction;
    }

    /**
     * Determines the minimal node order of an unlayered graph.
     * 
     * @return model order or Integer.MAX_VALUE if no node with a model order could be found.
     */
    public static int getMinimalModelOrder(final LGraph graph) {
        int order = Integer.MAX_VALUE;
        for (LNode node : graph.getLayerlessNodes()) {
            if (node.hasProperty(InternalProperties.MODEL_ORDER)) {
                order = Math.min(order, node.getProperty(InternalProperties.MODEL_ORDER));
            }
        }
        return order;
    }
    
}
