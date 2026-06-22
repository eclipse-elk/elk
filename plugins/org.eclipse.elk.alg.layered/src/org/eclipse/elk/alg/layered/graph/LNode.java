/*******************************************************************************
 * Copyright (c) 2010, 2019 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.alg.layered.graph;

import java.util.Collections;
import java.util.EnumMap;
import java.util.List;

import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.layered.options.PortType;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.options.PortSide;
import org.eclipse.elk.core.util.Pair;

import java.util.ArrayList;
import java.util.function.Predicate;

/**
 * A node in a layered graph.
 */
public final class LNode extends LShape {
    
    /**
     * Definition of node types used in the layered approach.
     */
    public static enum NodeType {
        
        /** a normal node is created from a node of the original graph. */
        NORMAL,
        /** a dummy node created to split a long edge. */
        LONG_EDGE,
        /** a node representing an external port. */
        EXTERNAL_PORT,
        /** a dummy node created to cope with ports at the northern or southern side. */
        NORTH_SOUTH_PORT,
        /** a dummy node to represent a mid-label on an edge. */
        LABEL,
        /** a dummy node representing a breaking point used to 'wrap' graphs. */
        BREAKING_POINT,
        /** a dummy node serving as a placeholder to reserve space when 'unzipping' graphs. 
         * this is used when there are no edges. */
        PLACEHOLDER,
        /** a placeholder node that can't be shifted when 'unzipping' graphs. this is used in front of extra edges. */
        NONSHIFTING_PLACEHOLDER;
        
        /**
         * Return the color used when writing debug output graphs. The colors are given as strings of
         * the form "#RGB", where each component is given as a two-digit hexadecimal value.
         * 
         * @return the color string
         */
        public String getColor() {
            switch (this) {
            case EXTERNAL_PORT: return "#cc99cc";
            case LONG_EDGE: return "#eaed00";
            case NORTH_SOUTH_PORT: return "#0034de";
            case LABEL: return "#75c3c3";
            case BREAKING_POINT: return "#eeeeff"; 
            default: return "#eeeeee";
            }
        }

    }
    
    /** the serial version UID. */
    private static final long serialVersionUID = -4272570519129722541L;
    
    /** the containing graph. */
    private LGraph graph;
    /** the containing layer. */
    private Layer layer;
    /** the node's node type. */
    private NodeType type = NodeType.NORMAL;
    /** the ports of the node. */
    private final List<LPort> ports = new ArrayList<>(6);
    /** this node's labels. */
    private final List<LLabel> labels = new ArrayList<>(2);
    /** the nested graph, or @{code null}. */
    private LGraph nestedGraph;
    /** the margin area around this node. */
    private final LMargin margin = new LMargin();
    /** the padding inside this node, usually reserved for port and label placement. */
    private final LPadding padding = new LPadding();
    /** Stores beginning and end indices of ports belonging to a side. */
    private EnumMap<PortSide, Pair<Integer, Integer>> portSideIndices;
    /** Stores whether port sides have been cached. */
    private boolean portSidesCached = false;

    /**
     * Creates a node.
     * 
     * @param graph the graph for which the node is created 
     */
    public LNode(final LGraph graph) {
        this.graph = graph;
    }

    /**
     * Returns the layer that owns this node.
     * 
     * @return the owning layer
     */
    public Layer getLayer() {
        return layer;
    }

    /**
     * Sets the owning layer and adds itself to the end of the layer's list of nodes. If the node
     * was previously in another layer, it is removed from that layer's list of nodes. Be careful
     * not to use this method while iterating through the nodes list of the old layer nor of the new
     * layer, since that could lead to {@link java.util.ConcurrentModificationException}s.
     * 
     * @param thelayer
     *            the owner to set
     */
    public void setLayer(final Layer thelayer) {
        if (this.layer != null) {
            this.layer.getNodes().remove(this);
        }
        
        this.layer = thelayer;
        
        if (this.layer != null) {
            this.layer.getNodes().add(this);
        }
    }
    
    /**
     * Returns the graph that contains this node.
     * 
     * @return the containing graph
     */
    public LGraph getGraph() {
        if (graph == null && layer != null) {
            return layer.getGraph();
        }
        return graph;
    }
    
    /**
     * Set the containing graph. This method must not be used when the layer has already been
     * assigned. The graphs' lists of nodes are <em>not</em> automatically updated.
     * 
     * @param newGraph
     *            the new containing graph
     */
    public void setGraph(final LGraph newGraph) {
        assert layer == null;
        this.graph = newGraph;
    }
    
    /**
     * Sets the containing layer and adds itself to the layer's list of nodes at the specified
     * position. If the node was previously in another layer, it is removed from that layer's list
     * of nodes. Be careful not to use this method while iterating through the nodes list of the old
     * layer nor of the new layer, since that could lead to
     * {@link java.util.ConcurrentModificationException}s.
     * 
     * @param index
     *            where the node should be inserted in the layer. Must be {@code >= 0} and
     *            {@code <= layer.getNodes().size()}.
     * @param newlayer
     *            the new layer
     */
    public void setLayer(final int index, final Layer newlayer) {
        if (newlayer != null && (index < 0 || index > newlayer.getNodes().size())) {
            throw new IllegalArgumentException("index must be >= 0 and <= layer node count");
        }
        
        if (this.layer != null) {
            this.layer.getNodes().remove(this);
        }
        
        this.layer = newlayer;
        
        if (newlayer != null) {
            newlayer.getNodes().add(index, this);
        }
    }
    
    /**
     * Returns the node's node type. Parts of the algorithm will treat nodes of different types
     * differently.
     * 
     * @return the node's node type.
     */
    public NodeType getType() {
        return type;
    }
    
    /**
     * Sets this node's node type.
     * 
     * @param type
     *            the node's new node type.
     */
    public void setType(final NodeType type) {
        this.type = type;
    }

    /**
     * Returns the list of ports of this node. Note that all edges are connected to specific ports,
     * even if the original diagram does not have any ports. Before the crossing minimization phase
     * has passed, the port order in this list is arbitrary. After crossing minimization the order
     * of ports corresponds to the clockwise order in which they are drawn, starting with the north
     * side. Hence the order is
     * <ul>
     *   <li>north ports from left to right,</li>
     *   <li>east ports from top to bottom,</li>
     *   <li>south ports from right to left,</li>
     *   <li>west port from bottom to top.</li>
     * </ul>
     * 
     * @return the ports of this node
     */
    public List<LPort> getPorts() {
        return ports;
    }
    
    /**
     * Returns an iterable for all ports of given type.
     * 
     * @param portType
     *            a port type
     * @return an iterable for the ports of given type
     */
    public Iterable<LPort> getPorts(final PortType portType) {
        switch (portType) {
        case INPUT:
            return LGraphUtil.filter(ports, LPort.INPUT_PREDICATE);
        case OUTPUT:
            return LGraphUtil.filter(ports, LPort.OUTPUT_PREDICATE);
        default:
            return Collections.emptyList();
        }
    }
    
    /**
     * Returns an iterable for all ports of given side.
     * 
     * @param side a port side
     * @return an iterable for the ports of given side
     */
    public Iterable<LPort> getPorts(final PortSide side) {
        switch (side) {
        case NORTH:
            return LGraphUtil.filter(ports, LPort.NORTH_PREDICATE);
        case EAST:
            return LGraphUtil.filter(ports, LPort.EAST_PREDICATE);
        case SOUTH:
            return LGraphUtil.filter(ports, LPort.SOUTH_PREDICATE);
        case WEST:
            return LGraphUtil.filter(ports, LPort.WEST_PREDICATE);
        default:
            return Collections.emptyList();
        }
    }
    
    /**
     * Returns a sublist view for all ports of given side. WARNING: Use this only after port sides are fixed! This is
     * currently the case after running the {@link org.eclipse.elk.alg.layered.intermediate.PortListSorter}.
     * Non-structural changes to this list are reflected in the original list. A structural modification is any
     * operation that adds or deletes one or more elements; merely setting the value of an element is not a structural
     * modification. Sublist indices can be cached using {@link LNode#cachePortSides()}.
     * 
     * @param side
     *            a port side
     * @return an iterable for the ports of given side
     */
    public List<LPort> getPortSideView(final PortSide side) {
        if (!portSidesCached) {
            // If not explicitly cached, this will be repeated each time. However, this has the same complexity as
            // filtering by side.
            findPortIndices();
        }
        Pair<Integer, Integer> indices = portSideIndices.get(side);
        if (indices == null) {
            return Collections.emptyList();
        } else {
            // We must create a new sublist each time, because the order of the ports on one side can change.
            return ports.subList(indices.getFirst(), indices.getSecond());
        }
    }

    /**
     * Returns an iterable for all ports of a given type and side.
     * 
     * @param portType a port type.
     * @param side a port side.
     * @return an iterable for the ports of the given type and side.
     */
    public Iterable<LPort> getPorts(final PortType portType, final PortSide side) {
        Predicate<LPort> typePredicate = null;
        switch (portType) {
        case INPUT:
            typePredicate = LPort.INPUT_PREDICATE;
            break;
        case OUTPUT:
            typePredicate = LPort.OUTPUT_PREDICATE;
            break;
        }
        
        Predicate<LPort> sidePredicate = null;
        switch (side) {
        case NORTH:
            sidePredicate = LPort.NORTH_PREDICATE;
            break;
        case EAST:
            sidePredicate = LPort.EAST_PREDICATE;
            break;
        case SOUTH:
            sidePredicate = LPort.SOUTH_PREDICATE;
            break;
        case WEST:
            sidePredicate = LPort.WEST_PREDICATE;
            break;
        }
        
        if (typePredicate != null && sidePredicate != null) {
            return LGraphUtil.filter(ports, typePredicate.and(sidePredicate));
        } else {
            return Collections.emptyList();
        }
    }
    
    /**
     * Returns an iterable for all inomcing edges.
     * 
     * @return an iterable for all incoming edges.
     */
    public Iterable<LEdge> getIncomingEdges() {
        List<Iterable<LEdge>> iterables = new ArrayList<>();
        for (LPort port : ports) {
            iterables.add(port.getIncomingEdges());
        }
        
        return LGraphUtil.concat(iterables);
    }
    
    /**
     * Returns an iterable for all outgoing edges.
     * 
     * @return an iterable for all outgoing edges.
     */
    public Iterable<LEdge> getOutgoingEdges() {
        List<Iterable<LEdge>> iterables = new ArrayList<>();
        for (LPort port : ports) {
            iterables.add(port.getOutgoingEdges());
        }
        
        return LGraphUtil.concat(iterables);
    }
    
    /**
     * Returns an iterable for all connected edges, both incoming and outgoing.
     * 
     * @return an iterable for all connected edges.
     */
    public Iterable<LEdge> getConnectedEdges() {
        List<Iterable<LEdge>> iterables = new ArrayList<>();
        for (LPort port : ports) {
            iterables.add(port.getConnectedEdges());
        }
        
        return LGraphUtil.concat(iterables);
    }
    
    /**
     * Returns this node's labels.
     * 
     * @return this node's labels.
     */
    public List<LLabel> getLabels() {
        return labels;
    }
    
    /**
     * Returns a graph that is nested in this node, if any.
     * 
     * @return the nested graph, or @{code null}
     */
    public LGraph getNestedGraph() {
        return nestedGraph;
    }
    
    /**
     * Sets the nested graph.
     * 
     * @param nestedGraph the new nested graph
     */
    public void setNestedGraph(final LGraph nestedGraph) {
        this.nestedGraph = nestedGraph;
    }
    
    /**
     * Returns the node's margin. The margin is the space around the node that is to be reserved
     * for ports and labels.
     * 
     * <p>The margin is not automatically updated. Rather, the margin has to be calculated once
     * the port and label positions are fixed. Usually this is right before the node placement
     * starts.</p>
     *  
     * @return the node's margin. May be modified.
     */
    public LMargin getMargin() {
        return margin;
    }
    
    /**
     * Returns the node's padding. The padding describes the area inside the node that is used by
     * ports, port labels, and node labels.
     * 
     * <p>The padding is not automatically updated. Rather, the paddings have to be calculated
     * once the port and label positions are fixed. Usually this is right before node placement
     * starts.</p>
     * 
     * @return the node's padding. May be modified.
     */
    public LPadding getPadding() {
        return padding;
    }
    
    /**
     * Returns the index of the node in the containing layer's list of nodes.
     * Note that this method has linear running time in the number of nodes,
     * so use it with caution.
     * 
     * @return the index of this node, or -1 if the node has no owner
     */
    public int getIndex() {
        if (layer == null) {
            return -1;
        } else {
            return layer.getNodes().indexOf(this);
        }
    }
    
    /**
     * Converts the position of this node from coordinates relative to the parent node's border to
     * coordinates relative to that node's content area. The content area is the parent node border
     * minus padding minus offset.
     * 
     * @param horizontal
     *            if {@code true}, the x coordinate will be translated.
     * @param vertical
     *            if {@code true}, the y coordinate will be translated.
     * @throws IllegalStateException
     *             if the node is not assigned to a layer in a layered graph.
     */
    public void borderToContentAreaCoordinates(final boolean horizontal, final boolean vertical) {
        LGraph thegraph = getGraph();
        
        LPadding graphPadding = thegraph.getPadding();
        KVector offset = thegraph.getOffset();
        KVector pos = getPosition();
        
        if (horizontal) {
            pos.x = pos.x - graphPadding.left - offset.x;
        }
        
        if (vertical) {
            pos.y = pos.y - graphPadding.top - offset.y;
        }
    }

    /**
     * Returns the position of this node's interactive reference point. This position depends on the
     * graph's {@link LayeredOptions#INTERACTIVE_REFERENCE_POINT} property. It determines on which
     * basis node positions are compared with each other in interactive layout phases.
     * 
     * @return the node's anchor point position.
     */
    public KVector getInteractiveReferencePoint() {
        switch (getGraph().getProperty(LayeredOptions.INTERACTIVE_REFERENCE_POINT)) {
        case CENTER:
            KVector nodePos = getPosition();
            KVector nodeSize = getSize();
            
            return new KVector(nodePos.x + nodeSize.x / 2.0, nodePos.y + nodeSize.y / 2.0);
        
        case TOP_LEFT:
            return new KVector(getPosition());
        
        default:
            // This shouldn't happen.
            return null;
        }
    }

    /**
     * After port sides are fixed (after running the {@link org.eclipse.elk.alg.layered.intermediate.PortListSorter} the
     * index region of ports for each side will stay constant. These are cached from the port list sorter using this
     * method and a sublist view created.
     */
    public void cachePortSides() {
        portSidesCached = true;
        findPortIndices();
    }

    private void findPortIndices() {
        portSideIndices = new EnumMap<>(PortSide.class);
        int firstIndexForCurrentSide = 0;
        PortSide currentSide = PortSide.NORTH;
        int currentIndex = 0;
        for (; currentIndex < ports.size(); currentIndex++) {
            LPort port = ports.get(currentIndex);
            if (port.getSide() != currentSide) {
                if (firstIndexForCurrentSide != currentIndex) {
                    portSideIndices.put(currentSide,
                            Pair.of(firstIndexForCurrentSide, currentIndex));
                }
                currentSide = port.getSide();
                firstIndexForCurrentSide = currentIndex;
            }
        }
        portSideIndices.put(currentSide, Pair.of(firstIndexForCurrentSide, currentIndex));
    }

    @Override
    public String getDesignation() {
        if (!labels.isEmpty() && !(labels.get(0).getText() == null || labels.get(0).getText().isEmpty())) {
            return labels.get(0).getText();
        }
        String id = super.getDesignation();
        if (id != null) {
            return id;
        }
        return Integer.toString(getIndex());
    }

    
    /**
     * Checks if the labels represented by the given this label dummy are to be placed inline.
     * 
     * @return True, if the node is a label and also an inline label.
     */
    public boolean isInlineEdgeLabel() {        
        return this.getType() == NodeType.LABEL && this.getProperty(InternalProperties.REPRESENTED_LABELS).stream()
                .allMatch(label -> label.getProperty(LayeredOptions.EDGE_LABELS_INLINE));
    }

    @Override
    public String toString() {
        StringBuilder result = new StringBuilder();
        result.append("n");
        if (type != NodeType.NORMAL) {
            result.append("(").append(type.toString().toLowerCase()).append(")");
        }
        result.append("_").append(getDesignation());
        return result.toString();
    }

}
