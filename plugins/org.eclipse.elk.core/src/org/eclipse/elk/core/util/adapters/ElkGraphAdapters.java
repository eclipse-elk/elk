/*******************************************************************************
 * Copyright (c) 2014, 2020 Kiel University and others.
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
package org.eclipse.elk.core.util.adapters;

import java.util.Comparator;
import java.util.List;

import org.eclipse.elk.core.math.ElkMargin;
import org.eclipse.elk.core.math.ElkPadding;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.options.LabelSide;
import org.eclipse.elk.core.options.PortSide;
import org.eclipse.elk.core.util.adapters.GraphAdapters.EdgeAdapter;
import org.eclipse.elk.core.util.adapters.GraphAdapters.GraphAdapter;
import org.eclipse.elk.core.util.adapters.GraphAdapters.GraphElementAdapter;
import org.eclipse.elk.core.util.adapters.GraphAdapters.LabelAdapter;
import org.eclipse.elk.core.util.adapters.GraphAdapters.NodeAdapter;
import org.eclipse.elk.core.util.adapters.GraphAdapters.PortAdapter;
import org.eclipse.elk.graph.ElkConnectableShape;
import org.eclipse.elk.graph.ElkEdge;
import org.eclipse.elk.graph.ElkGraphElement;
import org.eclipse.elk.graph.ElkLabel;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.ElkPort;
import org.eclipse.elk.graph.ElkShape;
import org.eclipse.elk.graph.properties.IProperty;
import org.eclipse.elk.graph.properties.Property;
import org.eclipse.elk.graph.util.ElkGraphUtil;
import org.eclipse.emf.common.util.ECollections;

import com.google.common.collect.Lists;

/**
 * Contains implementations of the {@link GraphAdapters} interfaces for the ElkGraph. To obtain an
 * adapter for a full ElkGraph, simply call {@link #adapt(ElkNode)}. To obtain an adapter only for a
 * single node, call {@link #adaptSingleNode(ElkNode)}.
 * 
 * @author uru
 */
public final class ElkGraphAdapters {

    private ElkGraphAdapters() {
        throw new IllegalStateException("Private constructor instantiation! Bad!");
    }

    /**
     * Creates the necessary adapters for the ElkGraph rooted at the given node.
     * 
     * @param graph
     *            the graph that should be wrapped in an adapter
     * @return an {@link ElkGraphAdapter} for the passed graph.
     */
    public static ElkGraphAdapter adapt(final ElkNode graph) {
        return new ElkGraphAdapter(graph);
    }

    /**
     * Creates a single node adapter for the given node.
     * 
     * @param node
     *            the node that should be wrapped in an adapter
     * @return an {@link ElkNodeAdapter} for the passed node.
     */
    public static ElkNodeAdapter adaptSingleNode(final ElkNode node) {
        return new ElkNodeAdapter(node.getParent() == null ? null : adapt(node.getParent()), node);
    }


    /**
     * Creates a single port adapter for the given port.
     * 
     * @param port
     *            the port that should be wrapped in an adapter
     * @return an {@link ElkPortAdapter} for the passed port.
     */
    public static ElkPortAdapter adaptSinglePort(final ElkPort port) {
        return new ElkPortAdapter(port);
    }
    
    /**
     * Implements basic adpater functionality for {@link ElkGraphElement}s.
     * 
     * @param <T>
     *            the type of the underlying graph element.
     */
    private abstract static class AbstractElkGraphElementAdapter<T extends ElkShape> implements
            GraphElementAdapter<T> {
        
        private static final IProperty<Double> OFFSET_PROXY = new Property<>(CoreOptions.PORT_BORDER_OFFSET, 0.0);
        
        // let the elements be accessed by extending classes
        // CHECKSTYLEOFF VisibilityModifier
        /** The wrapped element. */
        protected T element;
        // CHECKSTYLEON VisibilityModifier
        /**
         * Internally used versatile data field. Can be used for arbitrary information.
         * No assumptions about its value or validity should be made.
         */
        private int id;
        
        
        /**
         * Creates a new adapter for the given graph element.
         * 
         * @param element
         *            the element to be wrapped in this adapter.
         */
        protected AbstractElkGraphElementAdapter(final T element) {
            this.element = element;
        }


        @SuppressWarnings("unchecked")
        public <P> P getProperty(final IProperty<P> prop) {
            // the node spacing implementation requires a default value for the offset property
            if (prop.equals(CoreOptions.PORT_BORDER_OFFSET)) {
                return (P) element.getProperty(OFFSET_PROXY);
            }
            
            return element.getProperty(prop);
        }
        
        @Override
        public <P> boolean hasProperty(final IProperty<P> prop) {
            return element.hasProperty(prop);
        }

        @Override
        public KVector getPosition() {
            return new KVector(element.getX(), element.getY());
        }

        @Override
        public KVector getSize() {
            return new KVector(element.getWidth(), element.getHeight());
        }

        @Override
        public void setSize(final KVector size) {
            element.setWidth(size.x);
            element.setHeight(size.y);
        }

        @Override
        public void setPosition(final KVector pos) {
            element.setX(pos.x);
            element.setY(pos.y);
        }

        public ElkPadding getPadding() {
            ElkPadding elkPadding = element.getProperty(CoreOptions.PADDING);
            if (elkPadding == null) {
                return new ElkPadding();
            } else {
                return new ElkPadding(elkPadding);
            }
        }

        public void setPadding(final ElkPadding padding) {
            element.setProperty(CoreOptions.PADDING, new ElkPadding(padding));
        }

        public ElkMargin getMargin() {
            ElkMargin margins = element.getProperty(CoreOptions.MARGINS);
            if (margins == null) {
                margins = new ElkMargin();
            }
            return margins;
        }

        public void setMargin(final ElkMargin margin) {
            // analog to the padding case, we copy the margins object here
            ElkMargin newMargin = new ElkMargin(margin); 
            element.setProperty(CoreOptions.MARGINS, newMargin);
        }

        @Override
        public int getVolatileId() {
            return id;
        }

        @Override
        public void setVolatileId(final int volatileId) {
            this.id = volatileId;
        }
    }

    /**
     * Adapter for ElkGraphs rooted at a given node.
     */
    public static final class ElkGraphAdapter extends AbstractElkGraphElementAdapter<ElkNode> implements
            GraphAdapter<ElkNode> {
        
        /** cached list of child node adapters. */
        private List<NodeAdapter<?>> childNodes = null;
        
        /**
         * Creates a new adapter for the ElkGraph rooted at the given node.
         * 
         * @param node root of the ElkGraph to be adapted.
         */
        private ElkGraphAdapter(final ElkNode node) {
            super(node);
        }

        @Override
        public Iterable<NodeAdapter<?>> getNodes() {
            if (childNodes == null) {
                childNodes = Lists.newArrayListWithExpectedSize(element.getChildren().size());
                for (ElkNode n : element.getChildren()) {
                    childNodes.add(new ElkNodeAdapter(this, n));
                }
            }
            return childNodes;
        }
    }

    /**
     * Adapter for {@link ElkNode}s.
     */
    public static final class ElkNodeAdapter extends AbstractElkGraphElementAdapter<ElkNode> implements
            NodeAdapter<ElkNode> {
        
        /** The graph adapter that created this node adapter. */
        private ElkGraphAdapter parentGraphAdapter = null;
        /** Cached list of label adapters. */
        private List<LabelAdapter<?>> labelAdapters = null;
        /** Cached list of port adapters. */
        private List<PortAdapter<?>> portAdapters = null;
        /** Cached list of edge adapters for incoming edges. */
        private List<EdgeAdapter<?>> incomingEdgeAdapters = null;
        /** Cached list of edge adapters for outgoing edges. */
        private List<EdgeAdapter<?>> outgoingEdgeAdapters = null;
        
        
        /**
         * Creates a new adapter for the given node.
         * 
         * @param parent
         *            the graph adapter that created this thing.
         * @param node
         *            the node to adapt.
         */
        private ElkNodeAdapter(final ElkGraphAdapter parent, final ElkNode node) {
            super(node);
            parentGraphAdapter = parent;
        }
        

        @Override
        public GraphAdapter<?> getGraph() {
            return parentGraphAdapter;
        }

        @Override
        public List<LabelAdapter<?>> getLabels() {
            if (labelAdapters == null) {
                labelAdapters = Lists.newArrayListWithExpectedSize(element.getLabels().size());
                for (ElkLabel l : element.getLabels()) {
                    labelAdapters.add(new ElkLabelAdapter(l));
                }
            }
            return labelAdapters;
        }

        @Override
        public List<PortAdapter<?>> getPorts() {
            if (portAdapters == null) {
                portAdapters = Lists.newArrayListWithExpectedSize(element.getPorts().size());
                for (ElkPort p : element.getPorts()) {
                    portAdapters.add(new ElkPortAdapter(p));
                }
            }
            return portAdapters;
        }

        @Override
        public Iterable<EdgeAdapter<?>> getIncomingEdges() {
            if (incomingEdgeAdapters == null) {
                incomingEdgeAdapters = Lists.newArrayList();
                for (ElkEdge e : ElkGraphUtil.allIncomingEdges(element)) {
                    incomingEdgeAdapters.add(new ElkEdgeAdapter(e));
                }
            }
            return incomingEdgeAdapters;
        }

        @Override
        public Iterable<EdgeAdapter<?>> getOutgoingEdges() {
            if (outgoingEdgeAdapters == null) {
                outgoingEdgeAdapters = Lists.newArrayList();
                for (ElkEdge e : ElkGraphUtil.allOutgoingEdges(element)) {
                    outgoingEdgeAdapters.add(new ElkEdgeAdapter(e));
                }
            }
            return outgoingEdgeAdapters;
        }

        @Override
        public void sortPortList() {
            sortPortList(DEFAULT_PORTLIST_SORTER);
        }

        @Override
        @SuppressWarnings("unchecked")
        public void sortPortList(final Comparator<?> comparator) {
            // Iterate through the nodes of all layers
            if (element.getProperty(CoreOptions.PORT_CONSTRAINTS).isOrderFixed()) {
                ECollections.sort(element.getPorts(), (Comparator<ElkPort>) comparator);
            }
        }

        @Override
        public boolean isCompoundNode() {
            return !element.getChildren().isEmpty() || element.getProperty(CoreOptions.INSIDE_SELF_LOOPS_ACTIVATE);
        }
    }

    /**
     * Adapter for {@link ElkLabel}s.
     */
    private static final class ElkLabelAdapter extends AbstractElkGraphElementAdapter<ElkLabel> implements
            LabelAdapter<ElkLabel> {

        /**
         * Creates a new adapter for the given label.
         * 
         * @param label
         *            the label to adapt.
         */
        private ElkLabelAdapter(final ElkLabel label) {
            super(label);
        }
        
        
        @Override
        public LabelSide getSide() {
            return element.getProperty(LabelSide.LABEL_SIDE);
        }
        
        @Override
        public String getText() {
            return element.getText();
        }
        
    }

    /**
     * Adapter for {@link ElkPort}s.
     */
    private static final class ElkPortAdapter extends AbstractElkGraphElementAdapter<ElkPort>
            implements PortAdapter<ElkPort> {
        
        /** Cached list of label adapters. */
        private List<LabelAdapter<?>> labelAdapters = null;
        /** Cached list of edge adapters for incoming edges. */
        private List<EdgeAdapter<?>> incomingEdgeAdapters = null;
        /** Cached list of edge adapters for outgoing edges. */
        private List<EdgeAdapter<?>> outgoingEdgeAdapters = null;

        
        /**
         * Creates a new adapter for the given port.
         * 
         * @param port the port to adapt.
         */
        private ElkPortAdapter(final ElkPort port) {
            super(port);
        }

        
        public PortSide getSide() {
            return element.getProperty(CoreOptions.PORT_SIDE);
        }

        public List<LabelAdapter<?>> getLabels() {
            if (labelAdapters == null) {
                labelAdapters = Lists.newArrayListWithExpectedSize(element.getLabels().size());
                for (ElkLabel l : element.getLabels()) {
                    labelAdapters.add(new ElkLabelAdapter(l));
                }
            }
            return labelAdapters;
        }

        public Iterable<EdgeAdapter<?>> getIncomingEdges() {
            if (incomingEdgeAdapters == null) {
                incomingEdgeAdapters = Lists.newArrayListWithCapacity(element.getIncomingEdges().size());
                for (ElkEdge e : element.getIncomingEdges()) {
                    incomingEdgeAdapters.add(new ElkEdgeAdapter(e));
                }
            }
            return incomingEdgeAdapters;
        }

        public Iterable<EdgeAdapter<?>> getOutgoingEdges() {
            if (outgoingEdgeAdapters == null) {
                outgoingEdgeAdapters = Lists.newArrayListWithCapacity(element.getOutgoingEdges().size());
                for (ElkEdge e : element.getOutgoingEdges()) {
                    outgoingEdgeAdapters.add(new ElkEdgeAdapter(e));
                }
            }
            return outgoingEdgeAdapters;
        }

        public boolean hasCompoundConnections() {
            ElkNode node = element.getParent();
            
            for (ElkEdge edge : element.getOutgoingEdges()) {
                for (ElkConnectableShape target : edge.getTargets()) {
                    if (ElkGraphUtil.isDescendant(ElkGraphUtil.connectableShapeToNode(target), node)) {
                        return true;
                        
                    } else if (ElkGraphUtil.connectableShapeToNode(target) == node
                            && edge.getProperty(CoreOptions.INSIDE_SELF_LOOPS_YO)) {
                        
                        // Inside self loops are treated as compound connections, too
                        return true;
                    }
                }
            }
            
            for (ElkEdge edge : element.getIncomingEdges()) {
                for (ElkConnectableShape source : edge.getSources()) {
                    if (ElkGraphUtil.isDescendant(ElkGraphUtil.connectableShapeToNode(source), node)) {
                        return true;
                    }
                }
            }
            
            return false;
        }
    }

    /**
     * Adapter for {@link ElkEdge}s.
     */
    private static final class ElkEdgeAdapter implements EdgeAdapter<ElkEdge> {
        
        /** The wrapped edge. */
        private ElkEdge element;
        /** Cached list of label adapters. */
        private List<LabelAdapter<?>> labelAdapters = null;

        
        /**
         * Creates a new adapter for the given edge.
         * 
         * @param edge
         *            the edge to adapt.
         */
        private ElkEdgeAdapter(final ElkEdge edge) {
            this.element = edge;
        }

        public Iterable<LabelAdapter<?>> getLabels() {
            if (labelAdapters == null) {
                labelAdapters = Lists.newArrayListWithExpectedSize(element.getLabels().size());
                for (ElkLabel l : element.getLabels()) {
                    labelAdapters.add(new ElkLabelAdapter(l));
                }
            }
            return labelAdapters;
        }
    }
    
    /**
     * The default comparator for ports. Ports are sorted by side (north, east, south, west) in
     * clockwise order, beginning at the top left corner.
     */
    public static final PortComparator DEFAULT_PORTLIST_SORTER = new PortComparator();
    
    /**
     * A comparator for ports. Ports are sorted by side (north, east, south, west) in clockwise order,
     * beginning at the top left corner.
     */
    public static class PortComparator implements Comparator<ElkPort> {
        
        @Override
        public int compare(final ElkPort port1, final ElkPort port2) {
            int ordinalDifference =
                    port1.getProperty(CoreOptions.PORT_SIDE).ordinal()
                            - port2.getProperty(CoreOptions.PORT_SIDE).ordinal();

            // Sort by side first
            if (ordinalDifference != 0) {
                return ordinalDifference;
            }

            // In case of equal sides, sort by port index property
            Integer index1 = port1.getProperty(CoreOptions.PORT_INDEX);
            Integer index2 = port2.getProperty(CoreOptions.PORT_INDEX);
            if (index1 != null && index2 != null) {
                int indexDifference = index1 - index2;
                if (indexDifference != 0) {
                    return indexDifference;
                }
            }

            // In case of equal index, sort by position
            switch (port1.getProperty(CoreOptions.PORT_SIDE)) {
            case NORTH:
                // Compare x coordinates
                return Double.compare(port1.getX(), port2.getX());

            case EAST:
                // Compare y coordinates
                return Double.compare(port1.getY(), port2.getY());

            case SOUTH:
                // Compare x coordinates in reversed order
                return Double.compare(port2.getX(), port1.getX());

            case WEST:
                // Compare y coordinates in reversed order
                return Double.compare(port2.getY(), port1.getY());

            default:
                // Port sides should not be undefined
                throw new IllegalStateException("Port side is undefined");
            }
        }
        
    }
}
